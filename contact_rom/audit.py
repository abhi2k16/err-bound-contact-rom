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
from .estimator import error_bounds, eta_adm_v2, estimator_terms, trace_norm_B
from .fom import FomSolution, solve_fom
from .repair import rank_repair


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


def solve_reduced_timed(rom: ContactROM, f_vec: np.ndarray, c_suite: np.ndarray) -> ReducedSolution:
    """Reduced contact solve with iteration count and wall time.

    Same algebra as ``ContactROM.solve`` / ``solve_reduced_contact``, but calls
    ``lcp_active_set`` directly because the suite wrapper hides the iteration count.
    """
    start = time.perf_counter()
    chol = cho_factor(0.5 * (rom.A_r + rom.A_r.T))
    b_red, c_red = rom.U.T @ f_vec, rom.U_lam.T @ c_suite
    Ainv_b, Ainv_Bt = cho_solve(chol, b_red), cho_solve(chol, rom.B_r.T)
    A_lcp = rom.B_r @ Ainv_Bt
    lam_r, ok, iters = lcp_active_set(0.5 * (A_lcp + A_lcp.T), rom.B_r @ Ainv_b - c_red)
    u_r = Ainv_b + Ainv_Bt @ lam_r
    seconds = time.perf_counter() - start
    return ReducedSolution(rom.U @ u_r, rom.U_lam @ lam_r, u_r, lam_r, bool(ok), int(iters), seconds)


def evaluate_query(problem: ContactProblem, rom: ContactROM, mu: np.ndarray, fom: FomSolution,
                   C_lambda: float) -> dict:
    """All per-query metrics (paper-convention gap ``c - B u~``)."""
    K_mat, B_mat = problem.K_mat(mu), problem.B_mat(mu)
    f_vec, c_vec = problem.f_vec(mu), problem.c_vec(mu)
    res = solve_reduced_timed(rom, f_vec, -c_vec)
    K_lu = spla.splu(K_mat.tocsc())
    terms = estimator_terms(K_mat, B_mat, f_vec, c_vec, res.u, res.lam, K_lu, trace_norm_B(K_lu, B_mat))
    bounds = error_bounds(terms, float(np.linalg.norm(fom.lam_vec)), C_lambda)
    diff = fom.u_vec - res.u
    return {
        "ok": float(res.ok),
        "res_red": reduced_kkt_residual(rom, res, f_vec, -c_vec),
        "lam_min": float(res.lam.min()),
        "gap_viol": float(max(0.0, -(c_vec - B_mat @ res.u).min())),
        "err_energy": energy_error(K_mat, fom.u_vec, res.u),
        "err_abs": float(np.sqrt(diff @ (K_mat @ diff))),
        "err_lam": float(np.linalg.norm(res.lam - fom.lam_vec) / (np.linalg.norm(fom.lam_vec) + 1e-300)),
        "k_rom": float(res.iterations),
        "t_online": res.seconds,
        **{k: terms[k] for k in ("eta_eq", "eta_sign", "eta_gap", "eta_comp", "eta_adm")},
        "eta_adm_v2": eta_adm_v2(terms, C_lambda),
        **bounds,
    }


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
    meta = {"r_dual_eff": float(U_mat.shape[1]), "r_rep": float(r_rep), "cond_after": cond_after,
            "t_offline": time.perf_counter() - start}
    return rom, meta


def run_basis_audit_sets(problem: ContactProblem, cfg: Config, snaps: Snapshots, mu_sets: dict,
                         kinds: tuple, r_dual_list: tuple, repair: bool) -> dict:
    """Audit every (basis kind, r_lambda) on several query sets: ``{set_name: {(kind, r): entry}}``.

    Bases (and rank repair) are built once per (kind, r_lambda) and reused for every set.
    """
    V_pod = pod_primal_basis(snaps.U_snap, cfg.r_primal)
    C_lambda = float(np.linalg.norm(snaps.Lam_snap, axis=0).max())
    fom_sets = {name: [solve_fom(problem, mu, cfg.lcp_method) for mu in mus] for name, mus in mu_sets.items()}
    audits: dict = {name: {} for name in mu_sets}
    for kind in kinds:
        for r_dual in r_dual_list:
            rom, meta = build_reduced_model(problem, cfg, snaps, V_pod, kind, r_dual, repair)
            for name, mus in mu_sets.items():
                rows = [evaluate_query(problem, rom, mu, fom, C_lambda) for mu, fom in zip(mus, fom_sets[name])]
                entry = {key: np.array([row[key] for row in rows]) for key in rows[0]}
                entry.update({key: np.full(len(rows), val) for key, val in meta.items()})
                audits[name][(kind, r_dual)] = entry
    return audits


def run_basis_audit(problem: ContactProblem, cfg: Config, snaps: Snapshots, mu_test: np.ndarray,
                    kinds: tuple, r_dual_list: tuple, repair: bool) -> dict:
    """Single-set convenience wrapper around :func:`run_basis_audit_sets`."""
    return run_basis_audit_sets(problem, cfg, snaps, {"test": mu_test}, kinds, r_dual_list, repair)["test"]


def run_silent_failure_audit(problem: ContactProblem, cfg: Config, snaps: Snapshots, mu_test: np.ndarray) -> dict:
    """E1: signed-SVD vs NMF without rank repair (reproduces the silent failure)."""
    return run_basis_audit(problem, cfg, snaps, mu_test, ("svd", "nmf"), cfg.r_dual_list, repair=False)
