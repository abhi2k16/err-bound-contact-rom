"""Repeated wall-clock timings with variance for the break-even sweeps (B2 membrane, B3 two-body).

The manuscript's E4 sweep reports, per mesh size, the median over queries of best-of-5 timings. Here every raw
repeat is kept (``t[query, repeat]``) so the measurement noise and the query-to-query spread can be quantified:

* coefficient of variation across repeats (timing noise), and the spread across queries,
* percentile bootstrap intervals (queries and repeats resampled) for the measured speedup at each size,
* a bootstrap interval for the measured break-even size ``N*`` (resampling every size independently).

Point estimates use the same estimator as the manuscript (minimum over repeats, median over queries), so they agree
with the E4 tables up to the larger number of repeats. Timings must be run alone on the machine, one process at a time.
"""

from __future__ import annotations

import json
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from rom_engine.contact_dual_basis import ContactROM

from .audit import build_reduced_model, reduced_factors, solve_reduced_timed
from .bases import build_snapshots, pod_primal_basis
from .break_even_experiment import (REGIMES, predicted_speedups, sweep_problem, sweep_problem_b3,
                                    sweep_problem_b3_code)
from .config import Config
from .cost import ReducedSizes, find_break_even, problem_sizes
from .timing import assemble_factors, fom_solve_with_factor


def repeat_times(func, repeats: int) -> np.ndarray:
    """Wall time of each of ``repeats`` calls (seconds)."""
    out = np.empty(repeats)
    for i in range(repeats):
        start = time.perf_counter()
        func()
        out[i] = time.perf_counter() - start
    return out


def measure_size(cfg: Config, problem, seed_offset: int, repeats: int) -> dict | None:
    """Raw repeated timings for one mesh size (NMF dual basis, ``cfg.e4_r_dual``), or ``None`` without contact."""
    rng = np.random.default_rng(cfg.seed + seed_offset)
    snaps = build_snapshots(problem, problem.sample_mu(rng, cfg.e4_n_train), cfg.lcp_method)
    if snaps.Lam_snap.max() <= 0.0:
        return None
    mu_test = problem.sample_mu(rng, cfg.e4_n_test)
    V_pod = pod_primal_basis(snaps.U_snap, cfg.r_primal)
    rom, meta = build_reduced_model(problem, cfg, snaps, V_pod, cfg.e4_kind, cfg.e4_r_dual, repair=True)
    K0, B0 = problem.K_mat(mu_test[0]), problem.B_mat(mu_test[0])
    offline = assemble_factors(K0, -B0)
    t_rom, t_rom_aff, t_aff, t_full, k_rom, m_r, k_fom, m_act = [], [], [], [], [], [], [], []
    factors = reduced_factors(rom)
    for mu in mu_test:
        f_vec, c_s = problem.f_vec(mu), -problem.c_vec(mu)
        runs = [solve_reduced_timed(rom, f_vec, c_s) for _ in range(repeats)]
        t_rom.append([r.seconds for r in runs])
        t_rom_aff.append([solve_reduced_timed(rom, f_vec, c_s, factors).seconds for _ in range(repeats)])
        k_rom.append(runs[-1].iterations)
        m_r.append(int(np.sum(runs[-1].lam_r > 1e-12)))
        t_aff.append(repeat_times(lambda: fom_solve_with_factor(*offline, -B0, f_vec, c_s), repeats))
        t_full.append(repeat_times(lambda: fom_solve_with_factor(*assemble_factors(K0, -B0), -B0, f_vec, c_s), repeats))
        _, lam, iters = fom_solve_with_factor(*offline, -B0, f_vec, c_s)
        k_fom.append(iters)
        m_act.append(int(np.sum(lam > 1e-12)))
    V_mat, U_mat = rom.U, rom.U_lam
    t_proj = repeat_times(lambda: ContactROM(K0, -B0, V_mat, U_mat, method=cfg.lcp_method), repeats)
    BV = B0 @ V_mat
    u_r, lam_r, c1 = np.ones(V_mat.shape[1]), np.ones(U_mat.shape[1]), np.ones(B0.shape[0])

    def check() -> float:
        lam = U_mat @ lam_r
        gap = BV @ u_r - c1
        lam_pos = np.maximum(lam, 0.0)
        return (np.linalg.norm(np.minimum(lam, 0.0)) + np.linalg.norm(np.maximum(gap, 0.0))
                + np.sqrt(abs(lam_pos @ gap)))

    t_eta = repeat_times(check, repeats)
    sizes = problem_sizes(K0, B0)
    red = ReducedSizes(int(cfg.r_primal + meta["r_rep"]), int(meta["r_dual_eff"]), float(np.median(k_rom)),
                       float(np.median(m_r)))
    pred = predicted_speedups(cfg, sizes, red, float(np.median(k_fom)), float(np.median(m_act)))
    return {"N": sizes.n_dof, "n_c": sizes.n_contact, "r_u": red.r_u, "S_pred_aff": pred["aff"], "S_pred_var": pred["var"],
            "t_rom": np.array(t_rom).tolist(), "t_rom_aff": np.array(t_rom_aff).tolist(), "t_aff": np.array(t_aff).tolist(), "t_full": np.array(t_full).tolist(),
            "t_proj": t_proj.tolist(), "t_eta": t_eta.tolist()}


def run(cfg: Config, bench: str, repeats: int = 10) -> Path:
    """Time every size of the benchmark's sweep with ``repeats`` raw repeats and save to JSON."""
    if bench == "b2":
        jobs = [(g, sweep_problem(cfg, g), g) for g in cfg.e4_grids]
    else:
        jobs = [(g, sweep_problem_b3(cfg, g), g) for g in cfg.e4_grids_b3]
        jobs += [(c, sweep_problem_b3_code(cfg, c), c) for c in cfg.e4_codes_b3_small]
    rows = []
    for key, problem, offset in jobs:
        print(f"timing study {bench}: size parameter {key}", flush=True)
        row = measure_size(cfg, problem, offset, repeats)
        if row is not None:
            rows.append(row)
    rows.sort(key=lambda r: r["N"])
    path = Path(cfg.out_dir) / "timing_study" / f"timing_{bench}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"bench": bench, "repeats": repeats, "n_queries": cfg.e4_n_test, "rows": rows}),
                    encoding="utf-8")
    return path


def load(cfg: Config, bench: str) -> dict:
    return json.loads((Path(cfg.out_dir) / "timing_study" / f"timing_{bench}.json").read_text("utf-8"))


def speedup_from(t_aff, t_full, t_rom, t_proj, t_eta, t_rom_aff=None) -> dict:
    """Manuscript estimator: median over queries of the best-of-repeats time; ROM adds ``t_eta`` (and ``t_proj``).

    ``t_rom_aff`` is the reduced solve with precomputed factors (constant K); ``t_rom`` re-factorises (varying K).
    Older result files without ``t_rom_aff`` fall back to ``t_rom``.
    """
    fa, ff, rr = (np.median(np.min(a, axis=-1), axis=-1) for a in (t_aff, t_full, t_rom))
    ra = rr if t_rom_aff is None else np.median(np.min(t_rom_aff, axis=-1), axis=-1)
    eta = np.min(t_eta, axis=-1)
    return {"aff": fa / (ra + eta), "var": ff / (rr + eta + np.min(t_proj, axis=-1))}


def analyse(data: dict, n_boot: int = 1000, seed: int = 0) -> dict:
    """Per-size noise statistics and speedup intervals, plus bootstrap intervals for ``N*``."""
    rows = data["rows"]
    rng = np.random.default_rng(seed)
    n_vals = np.array([r["N"] for r in rows], float)
    arr = [{k: np.array(r.get(k, r["t_rom"]) if k == "t_rom_aff" else r[k])
            for k in ("t_aff", "t_full", "t_rom", "t_rom_aff", "t_proj", "t_eta")} for r in rows]
    per_size, boots = [], {reg: np.empty((n_boot, len(rows))) for reg in REGIMES}
    for i, (r, a) in enumerate(zip(rows, arr)):
        nq, nr = a["t_rom"].shape
        point = speedup_from(a["t_aff"], a["t_full"], a["t_rom"], a["t_proj"], a["t_eta"], a["t_rom_aff"])
        for b in range(n_boot):
            qi = rng.integers(0, nq, nq)
            ri = rng.integers(0, nr, (nq, nr))
            take = lambda t: np.take_along_axis(t[qi], ri, axis=1)
            rp = rng.integers(0, nr, nr)
            res = speedup_from(take(a["t_aff"]), take(a["t_full"]), take(a["t_rom"]), a["t_proj"][rp], a["t_eta"][rp],
                               take(a["t_rom_aff"]))
            for reg in REGIMES:
                boots[reg][b, i] = res[reg]
        cv = lambda t: float(np.median(np.std(t, axis=1) / np.mean(t, axis=1)))
        per_size.append({"N": r["N"], "cv_rom": cv(a["t_rom"]), "cv_aff": cv(a["t_aff"]), "cv_full": cv(a["t_full"]),
                         "spread_q_rom": float(np.std(np.min(a["t_rom"], axis=1)) / np.mean(np.min(a["t_rom"], axis=1))),
                         **{f"S_meas_{reg}": float(point[reg]) for reg in REGIMES},
                         **{f"S_pred_{reg}": r[f"S_pred_{reg}"] for reg in REGIMES},
                         **{f"S_ci_{reg}": [float(x) for x in np.percentile(boots[reg][:, i], [2.5, 97.5])] for reg in REGIMES}})
    crossings = {}
    for reg in REGIMES:
        point = find_break_even(n_vals, np.array([p[f"S_meas_{reg}"] for p in per_size]))
        pred = find_break_even(n_vals, np.array([p[f"S_pred_{reg}"] for p in per_size]))
        vals = np.array([find_break_even(n_vals, boots[reg][b]) for b in range(n_boot)])
        finite = vals[np.isfinite(vals)]
        crossings[reg] = {"N_meas": point, "N_pred": pred, "n_defined": int(finite.size), "n_boot": n_boot,
                          "N_meas_ci": [float(x) for x in np.percentile(finite, [2.5, 97.5])] if finite.size else [np.nan] * 2}
    return {"per_size": per_size, "crossings": crossings}
