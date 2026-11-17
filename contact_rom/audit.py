"""Reduced contact solve and per-query audit metrics (E1/E2/E3)."""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.linalg import cho_factor, cho_solve

from rom_engine.complementarity import lcp_active_set
from rom_engine.contact_dual_basis import ContactROM

from .bases import Snapshots, pod_primal_basis
from .benchmarks import ContactProblem
from .config import Config
from .dual_bases import BasisContext, make_dual_basis
from .estimator import error_bounds, eta_adm_v2, estimator_terms
from .query_context import QueryContext, build_query_context
from .repair import rank_repair
from .timing import best_of, time_eta_check


def reduced_kkt_residual(rom: ContactROM, rom_res, f_vec: np.ndarray, c_suite: np.ndarray) -> float:
    """Infinity-norm reduced KKT residual ``||R_red||_inf`` in reduced coordinates.

    Combines stationarity, reduced gap and multiplier sign violations and
    reduced complementarity. ``c_suite = -c`` (suite convention, with
    ``B_suite = -B``), so the reduced gap is ``B_r u_r - c_r >= 0``.
    """
    c_r = rom.U_lam.T @ c_suite
    gap_r = rom.B_r @ rom_res.u_r - c_r
    stationarity = rom.A_r @ rom_res.u_r - rom.U.T @ f_vec - rom.B_r.T @ rom_res.lam_r
    parts = [
        np.max(np.abs(stationarity)),
        max(0.0, -gap_r.min()),
        max(0.0, -rom_res.lam_r.min()),
        np.max(np.abs(rom_res.lam_r * gap_r)),
    ]
    return float(max(parts))


def energy_error(K_mat: sp.spmatrix, u_exact: np.ndarray, u_lift: np.ndarray) -> float:
    """Relative energy-norm error ``||u - u~||_K / ||u||_K``."""
    diff = u_exact - u_lift
    return float(np.sqrt(diff @ (K_mat @ diff)) / np.sqrt(u_exact @ (K_mat @ u_exact)))


@dataclass
class ReducedSolution:
    """Lifted reduced contact solution (suite sign convention for ``lam_r``)."""

    u: np.ndarray
    lam: np.ndarray
    u_r: np.ndarray
    lam_r: np.ndarray
    ok: bool
    iterations: int
    seconds: float


@dataclass
class ReducedFactors:
    """Query-independent reduced operators for constant ``K`` (offline stage of the reduced model).

    ``chol`` is the Cholesky factor of ``A_r``, ``Ainv_Bt = A_r^-1 B_r^T`` and ``A_lcp`` the symmetrised reduced
    Schur complement ``B_r A_r^-1 B_r^T``. With constant ``K`` none of them depends on the query, so they are
    computed once, exactly as the full-order affine model precomputes ``K^-1`` and ``B K^-1 B^T``.
    """

    chol: tuple
    Ainv_Bt: np.ndarray
    A_lcp: np.ndarray


def reduced_factors(rom: ContactROM) -> ReducedFactors:
    """Factorise ``A_r`` and assemble the reduced Schur complement once (valid only for constant ``K``)."""
    chol = cho_factor(0.5 * (rom.A_r + rom.A_r.T))
    Ainv_Bt = cho_solve(chol, rom.B_r.T)
    A_lcp = rom.B_r @ Ainv_Bt
    return ReducedFactors(chol, Ainv_Bt, 0.5 * (A_lcp + A_lcp.T))


def solve_reduced_timed(rom: ContactROM, f_vec: np.ndarray, c_suite: np.ndarray,
                        factors: ReducedFactors | None = None) -> ReducedSolution:
    """Reduced contact solve with iteration count and wall time.

    Same algebra as ``ContactROM.solve`` / ``solve_reduced_contact``, but calls
    ``lcp_active_set`` directly because the suite wrapper hides the iteration count.
    With ``factors=None`` the factorisation and the reduced Schur complement are rebuilt inside the timed
    region (parameter-varying ``K``); with precomputed ``factors`` (constant ``K``) only the query-dependent
    triangular solves and the active-set solve are timed.
    """
    start = time.perf_counter()
    fac = reduced_factors(rom) if factors is None else factors
    b_red, c_red = rom.U.T @ f_vec, rom.U_lam.T @ c_suite
    Ainv_b, Ainv_Bt = cho_solve(fac.chol, b_red), fac.Ainv_Bt
    lam_r, ok, iters = lcp_active_set(fac.A_lcp, rom.B_r @ Ainv_b - c_red)
    u_r = Ainv_b + Ainv_Bt @ lam_r
    seconds = time.perf_counter() - start
    return ReducedSolution(rom.U @ u_r, rom.U_lam @ lam_r, u_r, lam_r, bool(ok), int(iters), seconds)


def evaluate_query(ctx: QueryContext, rom: ContactROM, C_lambda: float,
                   factors: ReducedFactors | None = None) -> dict:
    """All per-query metrics (paper-convention gap ``c - B u~``) for one reduced model.

    ``t_online`` uses the precomputed ``factors`` when given (constant ``K``); ``t_online_refact`` always
    includes the reduced factorisation (what a parameter-varying ``K`` costs per query).
    """
    res = solve_reduced_timed(rom, ctx.f_vec, -ctx.c_vec, factors)
    t_refact = res.seconds if factors is None else solve_reduced_timed(rom, ctx.f_vec, -ctx.c_vec).seconds
    terms = estimator_terms(ctx.K_mat, ctx.B_mat, ctx.f_vec, ctx.c_vec, res.u, res.lam, ctx.K_lu, ctx.B_norm)
    bounds = error_bounds(terms, float(np.linalg.norm(ctx.fom.lam_vec)), C_lambda)
    diff = ctx.fom.u_vec - res.u
    return {
        "ok": float(res.ok),
        "res_red": reduced_kkt_residual(rom, res, ctx.f_vec, -ctx.c_vec),
        "lam_min": float(res.lam.min()),
        "gap_viol": float(max(0.0, -(ctx.c_vec - ctx.B_mat @ res.u).min())),
        "err_energy": energy_error(ctx.K_mat, ctx.fom.u_vec, res.u),
        "err_abs": float(np.sqrt(diff @ (ctx.K_mat @ diff))),
        "err_lam": float(np.linalg.norm(res.lam - ctx.fom.lam_vec) / (np.linalg.norm(ctx.fom.lam_vec) + 1e-300)),
        "k_rom": float(res.iterations),
        "m_r_act": float(np.sum(res.lam_r > 1e-12)),
        "t_online": res.seconds,
        "t_online_refact": t_refact,
        **{k: terms[k] for k in ("eta_eq", "eta_sign", "eta_gap", "eta_comp", "eta_adm")},
        "eta_adm_v2": eta_adm_v2(terms, C_lambda),
        **bounds,
    }


def rom_for_query(problem: ContactProblem, cfg: Config, rom: ContactROM, ctx: QueryContext) -> ContactROM:
    """The reduced model to use for this query: rebuilt from ``K(mu)`` when K varies, else the shared one."""
    if not problem.k_varies:
        return rom
    return ContactROM(ctx.K_mat, -ctx.B_mat, rom.U, rom.U_lam, method=cfg.lcp_method)


def build_reduced_model(problem: ContactProblem, cfg: Config, snaps: Snapshots, V_pod: np.ndarray,
                        kind: str, r_dual: int, repair: bool) -> tuple:
    """Build ``(rom, meta)`` for one (basis kind, r_lambda): dual basis, optional rank repair, timing."""
    start = time.perf_counter()
    K_mat, B_mat = problem.K_mat(snaps.mu_set[0]), problem.B_mat(snaps.mu_set[0])
    ctx = BasisContext(V_pod, B_mat, problem.elem_dofs)
    U_mat = make_dual_basis(kind, snaps.Lam_snap, r_dual, ctx)
    V_mat, r_rep, cond_after = V_pod, 0, float("nan")
    if repair:
        fixed = rank_repair(V_pod, U_mat, K_mat, B_mat, cfg.cond_max)
        V_mat, r_rep, cond_after = fixed.V_mat, fixed.r_rep, fixed.cond_after
    rom = ContactROM(K_mat, -B_mat, V_mat, U_mat, method=cfg.lcp_method)
    t_offline = time.perf_counter() - start
    meta = {"r_dual_eff": float(U_mat.shape[1]), "r_rep": float(r_rep), "cond_after": cond_after,
            "t_offline": t_offline,
            "t_proj": best_of(lambda: ContactROM(K_mat, -B_mat, V_mat, U_mat, method=cfg.lcp_method), cfg.e4_repeats),
            "t_eta": time_eta_check(B_mat, V_mat, U_mat, cfg.e4_repeats)}
    return rom, meta


def run_basis_audit_sets(problem: ContactProblem, cfg: Config, snaps: Snapshots, mu_sets: dict,
                         kinds: tuple, r_dual_list: tuple, repair: bool) -> dict:
    """Audit every (basis kind, r_lambda) on several query sets: ``{set_name: {(kind, r): entry}}``.

    Bases (and rank repair) are built once per (kind, r_lambda). Queries are the outer loop so each
    full-order solve and factorisation is done once and shared by all reduced models.
    """
    V_pod = pod_primal_basis(snaps.U_snap, cfg.r_primal)
    C_lambda = float(np.linalg.norm(snaps.Lam_snap, axis=0).max())
    models = {(kind, r_dual): build_reduced_model(problem, cfg, snaps, V_pod, kind, r_dual, repair)
              for kind in kinds for r_dual in r_dual_list}
    # constant K: the reduced factorisation is offline; varying K: rebuilt per query inside evaluate_query
    factors = {key: None if problem.k_varies else reduced_factors(model[0]) for key, model in models.items()}
    audits: dict = {}
    for name, mus in mu_sets.items():
        rows = {key: [] for key in models}
        for mu in mus:
            ctx = build_query_context(problem, mu)
            for key, (rom, _) in models.items():
                rows[key].append(evaluate_query(ctx, rom_for_query(problem, cfg, rom, ctx), C_lambda, factors[key]))
        audits[name] = {key: finalize_entry(rows[key], models[key][1]) for key in models}
    return audits


def finalize_entry(rows: list, meta: dict) -> dict:
    """Stack per-query rows into arrays and attach the scalar build metadata."""
    entry = {key: np.array([row[key] for row in rows]) for key in rows[0]}
    entry.update({key: np.full(len(rows), val) for key, val in meta.items()})
    return entry


def run_basis_audit(problem: ContactProblem, cfg: Config, snaps: Snapshots, mu_test: np.ndarray,
                    kinds: tuple, r_dual_list: tuple, repair: bool) -> dict:
    """Single-set convenience wrapper around :func:`run_basis_audit_sets`."""
    return run_basis_audit_sets(problem, cfg, snaps, {"test": mu_test}, kinds, r_dual_list, repair)["test"]


def run_silent_failure_audit(problem: ContactProblem, cfg: Config, snaps: Snapshots, mu_test: np.ndarray) -> dict:
    """E1: signed-SVD vs NMF without rank repair (reproduces the silent failure)."""
    return run_basis_audit(problem, cfg, snaps, mu_test, ("svd", "nmf"), cfg.r_dual_list, repair=False)
