"""Threshold-calibration study (B1, B2): detection metrics with uncertainty, validation-size sweep, threshold rules.

Per seed, one large pool of validation queries and one test set are audited once (all five dual bases, rank repair).
Smaller validation sets are random subsets of the pool, so the validation-size sweep costs no extra solves.

Positives are silent failures (reduced KKT converged but relative energy error above ``cfg.err_silent``);
negatives are accurate queries (error at or below ``cfg.err_silent``); queries in between are ignored,
as in the main tables.

Threshold rules (applied per basis, dual dimension and indicator):
  ``min09``  : ``tau = 0.9 * min`` of the indicator over validation silent failures (the manuscript's rule);
               undefined without a validation silent failure.
  ``rel``    : ``tau = err_silent * mean ||u||_K`` over the training snapshots. Needs no validation data; meaningful
               for the energy-norm-scaled indicator ``eta_adm_v2``.
  ``hybrid`` : ``min(min09, rel)`` when the validation set holds a silent failure, otherwise ``rel``. Always defined.
"""

from __future__ import annotations

import pickle
from dataclasses import replace
from math import sqrt
from pathlib import Path

import numpy as np

from .audit import run_basis_audit_sets
from .config import Config
from .dual_bases import DUAL_KINDS
from .estimator import tau_adm_from_snapshots
from .metrics import silent_failure_mask

INDICATORS = ("eta_adm", "eta_adm_v2")
RULES = ("min09", "rel", "hybrid")
KEEP = ("err_energy", "err_abs", "res_red", "eta_adm", "eta_adm_v2", "eta_sign", "eta_gap", "eta_comp", "bound1",
        "bound_Clam", "r_dual_eff", "r_rep")


def _slim(audit: dict) -> dict:
    return {key: {k: np.asarray(v) for k, v in entry.items() if k in KEEP} for key, entry in audit.items()}


def generate(cfg: Config, builders: dict, prepare, bench: str, seeds: tuple, n_val_pool: int, n_test: int) -> Path:
    """Audit the validation pool and test set for every seed; save slim arrays to ``results/calibration``."""
    out = Path(cfg.out_dir) / "calibration"
    out.mkdir(parents=True, exist_ok=True)
    for seed in seeds:
        path = out / f"audits_{bench}_seed{seed}.pkl"
        if path.exists():
            print(f"calibration: {bench} seed={seed} cached", flush=True)
            continue
        run_cfg = replace(cfg, seed=seed, n_val=n_val_pool, n_test=n_test)
        problem = builders[bench](run_cfg)
        print(f"calibration: {bench} seed={seed} (pool {n_val_pool}, test {n_test})", flush=True)
        snaps, mu_test, mu_val = prepare(run_cfg, problem)
        audits = run_basis_audit_sets(problem, run_cfg, snaps, {"test": mu_test, "val": mu_val}, DUAL_KINDS,
                                      run_cfg.r_dual_list_e3, repair=True)
        tau_rel = tau_adm_from_snapshots(snaps.U_snap, problem.K_mat(snaps.mu_set[0]), run_cfg.err_silent)
        payload = {"val": _slim(audits["val"]), "test": _slim(audits["test"]), "tau_rel": tau_rel,
                   "C_lambda": float(np.linalg.norm(snaps.Lam_snap, axis=0).max()), "seed": seed, "bench": bench}
        with open(path, "wb") as fh:
            pickle.dump(payload, fh)
    return out


def wilson(k: int, n: int, z: float = 1.96) -> tuple:
    """Wilson score interval for a binomial proportion (NaN when ``n = 0``)."""
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return centre - half, centre + half


def threshold(rule: str, val: dict, key: str, cfg: Config, tau_rel: float) -> float:
    silent = silent_failure_mask(val, cfg)
    tau_min = float(cfg.tau_margin * val[key][silent].min()) if silent.any() else float("nan")
    if rule == "min09":
        return tau_min
    if rule == "rel":
        return tau_rel
    return tau_rel if np.isnan(tau_min) else min(tau_min, tau_rel)


def confusion(test: dict, key: str, tau: float, cfg: Config) -> dict:
    """Counts of the 2x2 table on a test entry (NaN threshold -> not evaluable, all zero, ``evaluable`` False)."""
    pos = silent_failure_mask(test, cfg)
    neg = test["err_energy"] <= cfg.err_silent
    if np.isnan(tau):
        return {"tp": 0, "fn": 0, "fp": 0, "tn": 0, "evaluable": False, "n_pos": int(pos.sum())}
    flag = test[key] > tau
    return {"tp": int((pos & flag).sum()), "fn": int((pos & ~flag).sum()), "fp": int((neg & flag).sum()),
            "tn": int((neg & ~flag).sum()), "evaluable": True, "n_pos": int(pos.sum())}


def analyse(cfg: Config, bench: str, seeds: tuple, n_vals: tuple, draws: int = 20, rng_seed: int = 0) -> dict:
    """``{(kind, indicator, rule, n_val): stats}``: pooled confusion counts over seeds, entries and subset draws.

    A pool entry counts as *uncalibrated* when its rule gives no threshold (no silent failure in the subset).
    ``seed_recall`` lists the per-seed pooled recall (over entries and draws) for the spread across seeds.
    """
    rng = np.random.default_rng(rng_seed)
    folder = Path(cfg.out_dir) / "calibration"
    data = {s: pickle.load(open(folder / f"audits_{bench}_seed{s}.pkl", "rb")) for s in seeds}
    acc: dict = {}
    for s, payload in data.items():
        for key_entry, val_full in payload["val"].items():
            kind = key_entry[0]
            test = payload["test"][key_entry]
            n_pool = len(val_full["err_energy"])
            for n_val in n_vals:
                n_draw = 1 if n_val >= n_pool else draws
                for _ in range(n_draw):
                    idx = np.arange(n_pool) if n_val >= n_pool else rng.choice(n_pool, n_val, replace=False)
                    val = {k: v[idx] for k, v in val_full.items() if np.ndim(v) == 1 and len(v) == n_pool}
                    for ind in INDICATORS:
                        for rule in RULES:
                            tau = threshold(rule, val, ind, cfg, payload["tau_rel"])
                            c = confusion(test, ind, tau, cfg)
                            slot = acc.setdefault((kind, ind, rule, n_val), {"tp": 0, "fn": 0, "fp": 0, "tn": 0, "uncal": 0,
                                                                              "entries": 0, "pos_uncal": 0, "per_seed": {}})
                            slot["entries"] += 1
                            if not c["evaluable"]:
                                slot["uncal"] += 1
                                slot["pos_uncal"] += c["n_pos"]
                                continue
                            for f in ("tp", "fn", "fp", "tn"):
                                slot[f] += c[f]
                            ps = slot["per_seed"].setdefault(s, [0, 0])
                            ps[0] += c["tp"]
                            ps[1] += c["tp"] + c["fn"]
    stats = {}
    for key, a in acc.items():
        tp, fn, fp, tn = a["tp"], a["fn"], a["fp"], a["tn"]
        rec, spec = (tp / (tp + fn) if tp + fn else float("nan")), (tn / (tn + fp) if tn + fp else float("nan"))
        stats[key] = {
            "recall": rec, "recall_ci": wilson(tp, tp + fn), "specificity": spec, "spec_ci": wilson(tn, tn + fp),
            "fpr": 1 - spec if tn + fp else float("nan"), "fnr": 1 - rec if tp + fn else float("nan"),
            "precision": tp / (tp + fp) if tp + fp else float("nan"), "tp": tp, "fn": fn, "fp": fp, "tn": tn,
            "uncal_share": a["uncal"] / a["entries"], "pos_in_uncal": a["pos_uncal"],
            "seed_recall": {s: (v[0] / v[1] if v[1] else float("nan")) for s, v in a["per_seed"].items()},
        }
    return stats


def save_analysis(cfg: Config, bench: str, stats: dict) -> Path:
    path = Path(cfg.out_dir) / "calibration" / f"analysis_{bench}.pkl"
    with open(path, "wb") as fh:
        pickle.dump(stats, fh)
    return path


def print_analysis(stats: dict, bench: str, n_vals: tuple) -> None:
    kinds = sorted({k[0] for k in stats})
    for ind in INDICATORS:
        print(f"\n[{bench}] indicator {ind}; rows: basis, rule; cells per validation size: "
              f"recall [95% CI] / FPR / uncalibrated share")
        for kind in kinds:
            for rule in RULES:
                cells = []
                for n in n_vals:
                    s = stats.get((kind, ind, rule, n))
                    if s is None:
                        continue
                    lo, hi = s["recall_ci"]
                    cells.append(f"n={n}: {s['recall']:.3f} [{lo:.2f},{hi:.2f}] / {s['fpr']:.2f} / {s['uncal_share']:.2f}")
                print(f"  {kind:<4} {rule:<7} " + " | ".join(cells))
