"""E4: break-even cost model validation (predicted vs measured speedup, N sweep, crossover)."""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from .audit import build_reduced_model, evaluate_query, solve_reduced_timed
from .bases import Snapshots, build_snapshots, pod_primal_basis
from .benchmark_3d import build_b3
from .benchmarks import ContactProblem, build_b2
from .config import Config
from .cost import ProblemSizes, ReducedSizes, find_break_even, problem_sizes, speedup_pred
from .timing import FomTimings, time_fom_queries

REGIMES = {"aff": False, "var": True}   # name -> K parameter-varying?


def predicted_speedups(cfg: Config, sizes: ProblemSizes, red: ReducedSizes, k_fom: float, m_act: float) -> dict:
    """Eq. (16)/(17) speedups for the constant-K ("aff") and K-varying ("var") regimes."""
    return {name: speedup_pred(sizes, red, k_fom, m_act, cfg.c_fact, cfg.c_solve, cfg.c_proj, varies)
            for name, varies in REGIMES.items()}


def measured_speedups(fom_times: dict, t_rom_online: float, t_proj: float, t_eta: float) -> dict:
    """Wall-clock speedups; the K-varying ROM additionally pays the online projection ``t_proj``."""
    t_rom_aff = t_rom_online + t_eta
    return {"aff": fom_times["aff"] / t_rom_aff, "var": fom_times["full"] / (t_rom_aff + t_proj)}


def entry_row(cfg: Config, sizes: ProblemSizes, fom: FomTimings, entry: dict) -> dict:
    """E4 row for one (basis, r_lambda) audit entry."""
    red = ReducedSizes(int(cfg.r_primal + entry["r_rep"][0]), int(entry["r_dual_eff"][0]),
                       float(np.median(entry["k_rom"])), float(np.median(entry["m_r_act"])))
    k_fom, m_act = float(np.median(fom.k_fom)), float(np.median(fom.m_act))
    pred = predicted_speedups(cfg, sizes, red, k_fom, m_act)
    meas = measured_speedups({"aff": float(np.median(fom.t_affine)), "full": float(np.median(fom.t_full))},
                             float(np.median(entry["t_online"])), float(entry["t_proj"][0]), float(entry["t_eta"][0]))
    return {"r_u": red.r_u, "r_lambda": red.r_lambda, "k_rom": red.k_rom, "m_r": red.m_r_act,
            "S_pred_aff": pred["aff"], "S_meas_aff": meas["aff"], "S_pred_var": pred["var"], "S_meas_var": meas["var"]}


def run_e4_entries(cfg: Config, problem: ContactProblem, snaps: Snapshots, mu_test: np.ndarray, audit: dict) -> tuple:
    """Rows for every audit entry plus the FOM timing summary (``k_fom``, ``m_act`` medians)."""
    fom = time_fom_queries(problem, mu_test, cfg.e4_repeats)
    sizes = problem_sizes(problem.K_mat(snaps.mu_set[0]), problem.B_mat(snaps.mu_set[0]))
    rows = {key: entry_row(cfg, sizes, fom, entry) for key, entry in audit.items()}
    return rows, fom, sizes


def sweep_problem(cfg: Config, n_grid: int) -> ContactProblem:
    """Membrane with ``n_grid`` interior nodes; load scaled by ``(18 / n_grid)^2`` to keep a partial contact patch."""
    scale = (18.0 / n_grid) ** 2
    lo, hi = cfg.b2.load_range
    cfg_b2 = replace(cfg.b2, n_grid=n_grid, n_sub=max(4, round(2 * n_grid / 3)), load_range=(lo * scale, hi * scale))
    return build_b2(cfg_b2)


def sweep_one_size(cfg: Config, n_grid: int, builder=sweep_problem) -> dict | None:
    """Measured and predicted speedups for one mesh size; ``None`` if no contact occurs in the snapshots."""
    problem = builder(cfg, n_grid)
    rng = np.random.default_rng(cfg.seed + n_grid)
    snaps = build_snapshots(problem, problem.sample_mu(rng, cfg.e4_n_train), cfg.lcp_method)
    if snaps.Lam_snap.max() <= 0.0:
        return None
    mu_test = problem.sample_mu(rng, cfg.e4_n_test)
    V_pod = pod_primal_basis(snaps.U_snap, cfg.r_primal)
    rom, meta = build_reduced_model(problem, cfg, snaps, V_pod, cfg.e4_kind, cfg.e4_r_dual, repair=True)
    t_online, k_rom, m_r = [], [], []
    for mu in mu_test:
        f_vec, c_s = problem.f_vec(mu), -problem.c_vec(mu)
        runs = [solve_reduced_timed(rom, f_vec, c_s) for _ in range(cfg.e4_repeats)]
        t_online.append(min(run.seconds for run in runs))
        k_rom.append(runs[-1].iterations)
        m_r.append(int(np.sum(runs[-1].lam_r > 1e-12)))
    fom = time_fom_queries(problem, mu_test, cfg.e4_repeats)
    sizes = problem_sizes(problem.K_mat(mu_test[0]), problem.B_mat(mu_test[0]))
    red = ReducedSizes(int(cfg.r_primal + meta["r_rep"]), int(meta["r_dual_eff"]), float(np.median(k_rom)),
                       float(np.median(m_r)))
    k_fom, m_act = float(np.median(fom.k_fom)), float(np.median(fom.m_act))
    pred = predicted_speedups(cfg, sizes, red, k_fom, m_act)
    meas = measured_speedups({"aff": float(np.median(fom.t_affine)), "full": float(np.median(fom.t_full))},
                             float(np.median(t_online)), meta["t_proj"], meta["t_eta"])
    return {"N": sizes.n_dof, "n_c": sizes.n_contact, "b_w": sizes.bandwidth, "m_act": m_act, "r_u": red.r_u,
            "S_pred_aff": pred["aff"], "S_meas_aff": meas["aff"], "S_pred_var": pred["var"], "S_meas_var": meas["var"]}


def sweep_problem_b3(cfg: Config, n_xy: int) -> ContactProblem:
    """B3 with ``n_xy`` elements per side and layers per body scaled as ``0.3 n_xy`` (at least 2)."""
    nz = max(2, round(0.3 * n_xy))
    return build_b3(replace(cfg.b3, n_xy=n_xy, nz_up=nz, nz_low=nz))


def sweep_problem_b3_code(cfg: Config, code: int) -> ContactProblem:
    """B3 with ``code = 10 * n_xy + nz`` (elements per side, layers per body); used for very small meshes."""
    n_xy, nz = divmod(code, 10)
    return build_b3(replace(cfg.b3, n_xy=n_xy, nz_up=nz, nz_low=nz))


def run_sweep(cfg: Config, grids: tuple | None = None, builder=sweep_problem) -> list:
    """N sweep over ``grids`` (default ``cfg.e4_grids``); sizes without contact are skipped."""
    rows = [sweep_one_size(cfg, n_grid, builder) for n_grid in (cfg.e4_grids if grids is None else grids)]
    return [row for row in rows if row is not None]


def crossovers(rows: list) -> dict:
    """Break-even ``N*`` (S = 1) predicted and measured, per regime; NaN if no crossing in the swept range."""
    n_vals = np.array([row["N"] for row in rows], float)
    out = {}
    for regime in REGIMES:
        pred = find_break_even(n_vals, np.array([row[f"S_pred_{regime}"] for row in rows]))
        meas = find_break_even(n_vals, np.array([row[f"S_meas_{regime}"] for row in rows]))
        ratio = pred / meas if np.isfinite(pred) and np.isfinite(meas) else float("nan")
        out[regime] = {"N_pred": pred, "N_meas": meas, "ratio": ratio}
    return out


def criterion4_verdict(cross: dict) -> dict:
    """Criterion 4: predicted break-even within a factor 2 of measured (``None`` when undefined)."""
    verdict = {}
    for regime, values in cross.items():
        ratio = values["ratio"]
        verdict[regime] = None if not np.isfinite(ratio) else bool(0.5 <= ratio <= 2.0)
    return verdict
