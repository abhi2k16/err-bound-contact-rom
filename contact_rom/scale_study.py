"""Out-of-range and mesh-size studies on B1 and B2 (constant stiffness, five bases, ``r_lambda = 16``).

*Extrapolation.* Reduced models are built on the training box of each benchmark and queried inside the box
(validation and test) and outside it: loads 1.3-2.0 times the training maximum (B1, B2), a steeper obstacle (B1) and an
indenter centre moved outside the training box (B2, which also moves the contact patch). For each set we report the
error, the silent-failure share, the validity of the two forms of the Theorem 1 bound (exact multiplier norm and offline
``C_lambda``) and the performance of the interface screen with the hybrid threshold calibrated on the in-range
validation queries.

*Mesh size.* The same reduced-model pipeline on several mesh sizes of B1 (chain length) and B2 (membrane grid), with
the measured speedup against the full-order solver.
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from .audit import run_basis_audit_sets
from .bases import build_snapshots, pod_primal_basis
from .benchmarks import build_b1, build_b2
from .break_even_experiment import sweep_problem
from .calibration_study import threshold
from .config import B1Config, Config
from .dual_bases import DUAL_KINDS
from .estimator import tau_adm_from_snapshots
from .metrics import silent_failure_mask

R_DUAL = 16
SEEDS = (42, 43, 44, 45, 46)


def extrap_sets(cfg: Config, bench: str, rng: np.random.Generator, n: int) -> dict:
    """Out-of-range parameter sets ``{name: (n, n_p) array}``."""
    if bench == "b1":
        b = cfg.b1
        lo, hi = b.load_range[1] * 1.3, b.load_range[1] * 2.0
        load = np.column_stack([lo + rng.random(n) * (hi - lo), b.curvature_range[0] + rng.random(n) * (b.curvature_range[1] - b.curvature_range[0])])
        c_lo, c_hi = b.curvature_range[1] * 1.25, b.curvature_range[1] * 1.7
        geom = np.column_stack([b.load_range[0] + rng.random(n) * (b.load_range[1] - b.load_range[0]), c_lo + rng.random(n) * (c_hi - c_lo)])
        return {"load": load, "geometry": geom}
    b = cfg.b4 if bench == "b4" else cfg.b2
    lo, hi = b.load_range[1] * 1.3, b.load_range[1] * 2.0
    centre = b.center_range[0] + rng.random((n, 2)) * (b.center_range[1] - b.center_range[0])
    load = np.column_stack([lo + rng.random(n) * (hi - lo), centre])
    # centre outside the training box in each coordinate: [0.20, 0.30] or [0.70, 0.80]
    side = rng.integers(0, 2, (n, 2))
    out = np.where(side == 0, 0.20 + 0.10 * rng.random((n, 2)), 0.70 + 0.10 * rng.random((n, 2)))
    shift = np.column_stack([b.load_range[0] + rng.random(n) * (b.load_range[1] - b.load_range[0]), out])
    return {"load": load, "centre shift": shift}


def pooled_set(cfg: Config, entries: list, val_entries: list, taus: list) -> dict:
    """Pool one (kind, set) over seeds: error, silent share, bound validity, screen recall and false alarms."""
    err = np.concatenate([e["err_energy"] for e in entries])
    silent = np.concatenate([silent_failure_mask(e, cfg) for e in entries])
    good = err <= cfg.err_silent
    flagged = np.concatenate([e["eta_adm"] > t for e, t in zip(entries, taus)])
    return {"err_med": float(np.median(err)), "share_above_5pct": float(np.mean(err > cfg.err_silent)),
            "silent_share": float(silent.mean()),
            "bound1_valid": float(np.mean(np.concatenate([e["bound1"] >= e["err_abs"] * (1 - 1e-8) for e in entries]))),
            "bound_Clam_valid": float(np.mean(np.concatenate([e["bound_Clam"] >= e["err_abs"] for e in entries]))),
            "n": int(err.size), "n_silent": int(silent.sum()), "missed": int((silent & ~flagged).sum()),
            "false_alarm": float(flagged[good].mean()) if good.any() else float("nan"),
            "flagged_share": float(flagged.mean())}


def run_extrapolation(cfg: Config, builders: dict, bench: str) -> Path:
    per: dict = {}
    for seed in SEEDS:
        run_cfg = replace(cfg, seed=seed)
        problem = builders[bench](run_cfg)
        rng = np.random.default_rng(seed + 1)
        snaps = build_snapshots(problem, problem.sample_mu(rng, run_cfg.n_train), run_cfg.lcp_method)
        sets = {"test": problem.sample_mu(rng, run_cfg.n_test), "val": problem.sample_mu(rng, run_cfg.n_val)}
        sets.update(extrap_sets(run_cfg, bench, rng, 40))
        audits = run_basis_audit_sets(problem, run_cfg, snaps, sets, DUAL_KINDS, (R_DUAL,), repair=True)
        tau_rel = tau_adm_from_snapshots(snaps.U_snap, problem.K_mat(snaps.mu_set[0]), run_cfg.err_silent)
        for name in sets:
            if name == "val":
                continue
            for (kind, _), entry in audits[name].items():
                tau = threshold("hybrid", audits["val"][(kind, R_DUAL)], "eta_adm", run_cfg, tau_rel)
                slot = per.setdefault((kind, name), {"e": [], "v": [], "t": []})
                slot["e"].append(entry); slot["v"].append(audits["val"][(kind, R_DUAL)]); slot["t"].append(tau)
    rows = {f"{kind}|{name}": pooled_set(cfg, s["e"], s["v"], s["t"]) for (kind, name), s in per.items()}
    path = Path(cfg.out_dir) / "scale_study" / f"extrapolation_{bench}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"bench": bench, "seeds": list(SEEDS), "rows": rows}, indent=1), encoding="utf-8")
    return path


def mesh_problem(cfg: Config, bench: str, size: int):
    if bench == "b1":
        return build_b1(replace(cfg.b1, n_dof=size, n_contact=size // 2))
    return sweep_problem(cfg, size)


MESHES = {"b1": (30, 60, 120, 240, 480), "b2": (12, 18, 24, 32, 40)}


def run_mesh(cfg: Config, bench: str) -> Path:
    from .timing import time_fom_queries
    rows: dict = {}
    for size in MESHES[bench]:
        pooled: dict = {}
        for seed in SEEDS[:3]:
            run_cfg = replace(cfg, seed=seed, e4_repeats=5)
            problem = mesh_problem(run_cfg, bench, size)
            rng = np.random.default_rng(seed + 1)
            snaps = build_snapshots(problem, problem.sample_mu(rng, run_cfg.n_train), run_cfg.lcp_method)
            if snaps.Lam_snap.max() <= 0.0:
                continue
            mu_test = problem.sample_mu(rng, run_cfg.n_test)
            audit = run_basis_audit_sets(problem, run_cfg, snaps, {"test": mu_test}, ("svd", "nmf", "gis"), (R_DUAL,), repair=True)["test"]
            t_fom = float(np.median(time_fom_queries(problem, mu_test, 5).t_affine))
            for (kind, _), e in audit.items():
                slot = pooled.setdefault(kind, {"err": [], "res": [], "t_rom": [], "t_fom": [], "bound": [], "abs": [], "eta": []})
                slot["err"].append(e["err_energy"]); slot["res"].append(e["res_red"]); slot["t_rom"].append(e["t_online"] + float(e["t_eta"][0]))
                slot["t_fom"].append(t_fom); slot["bound"].append(e["bound1"]); slot["abs"].append(e["err_abs"])
        n_dof = problem.n_dof
        for kind, s in pooled.items():
            err, res = np.concatenate(s["err"]), np.concatenate(s["res"])
            rows[f"{kind}|{n_dof}"] = {"N": n_dof, "n_c": problem.n_contact, "err_med": float(np.median(err)),
                                       "share_above_5pct": float(np.mean(err > cfg.err_silent)),
                                       "silent_share": float(np.mean((res <= cfg.tol_kkt_red) & (err > cfg.err_silent))),
                                       "bound1_valid": float(np.mean(np.concatenate(s["bound"]) >= np.concatenate(s["abs"]) * (1 - 1e-8))),
                                       "S_aff": float(np.median(s["t_fom"]) / np.median(np.concatenate(s["t_rom"])))}
        print(f"mesh {bench} size {size}: N={n_dof}", flush=True)
    path = Path(cfg.out_dir) / "scale_study" / f"mesh_{bench}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"bench": bench, "rows": rows}, indent=1), encoding="utf-8")
    return path
