"""Dual-basis comparison at matched effective dimension, matched online cost and matched accuracy (B1, B2).

The five bases are audited on identical snapshots over a wide range of nominal dual dimensions (with primal rank
repair). Per entry we record the *effective* dual dimension (rows kept by RID, columns of ``U`` otherwise), the
enriched primal dimension ``r_u = r_primal + r_rep``, the accuracy, the Eq. (17) online FLOP count (constant ``K``,
so the projection term is zero) and the measured reduced-solve time. Entries from different bases are then compared
on three common axes instead of the nominal ``r_lambda``:

* matched effective dual dimension (nearest entry within a tolerance, otherwise reported as unavailable),
* matched accuracy (cheapest entry whose median error is at most a target),
* matched cost (most accurate entry within an online FLOP budget).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .audit import run_basis_audit_sets
from .config import Config
from .cost import ReducedSizes, flops_rom, problem_sizes
from .dual_bases import DUAL_KINDS

DIMS = {"b1": (4, 8, 12, 16, 20, 24, 30), "b2": (4, 8, 12, 16, 24, 32, 48, 60)}


def summarize_entry(cfg: Config, sizes, entry: dict, nominal: int) -> dict:
    k_rom, m_act = float(np.median(entry["k_rom"])), float(np.median(entry["m_r_act"]))
    r_u, r_l = int(cfg.r_primal + entry["r_rep"][0]), int(entry["r_dual_eff"][0])
    err = entry["err_energy"]
    return {"nominal": nominal, "r_eff": r_l, "r_u": r_u, "r_rep": int(entry["r_rep"][0]),
            "err_med": float(np.median(err)), "err_p90": float(np.quantile(err, 0.9)),
            "err_max": float(err.max()), "frac_above_5pct": float(np.mean(err > cfg.err_silent)),
            "t_online": float(np.median(entry["t_online"])), "k_rom": k_rom, "m_act": m_act,
            "flops": flops_rom(sizes, ReducedSizes(r_u, r_l, k_rom, m_act), cfg.c_proj, False),
            "t_offline": float(entry["t_offline"][0])}


def run(cfg: Config, problem, prepare, bench: str, n_test: int = 100) -> Path:
    from dataclasses import replace
    run_cfg = replace(cfg, n_test=n_test)
    snaps, mu_test, _ = prepare(run_cfg, problem)
    audit = run_basis_audit_sets(problem, run_cfg, snaps, {"test": mu_test}, DUAL_KINDS, DIMS[bench], repair=True)["test"]
    sizes = problem_sizes(problem.K_mat(snaps.mu_set[0]), problem.B_mat(snaps.mu_set[0]))
    rows = {kind: [summarize_entry(run_cfg, sizes, audit[(kind, r)], r) for r in DIMS[bench]] for kind in DUAL_KINDS}
    path = Path(cfg.out_dir) / "basis_comparison" / f"basis_comparison_{bench}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"bench": bench, "n_test": n_test, "n_contact": sizes.n_contact, "rows": rows}, indent=1),
                    encoding="utf-8")
    return path


def load(cfg: Config, bench: str) -> dict:
    return json.loads((Path(cfg.out_dir) / "basis_comparison" / f"basis_comparison_{bench}.json").read_text("utf-8"))


def matched_dimension(rows: dict, target: int, tol: float = 0.25) -> dict:
    """Per basis: entry whose effective dual dimension is nearest to ``target`` if within ``tol * target``, else None."""
    out = {}
    for kind, entries in rows.items():
        best = min(entries, key=lambda e: abs(e["r_eff"] - target))
        out[kind] = best if abs(best["r_eff"] - target) <= tol * target else None
    return out


def matched_accuracy(rows: dict, err_target: float, key: str = "flops") -> dict:
    """Per basis: cheapest entry (by ``key``) with median error <= ``err_target``; None if none reaches it."""
    out = {}
    for kind, entries in rows.items():
        ok = [e for e in entries if e["err_med"] <= err_target]
        out[kind] = min(ok, key=lambda e: e[key]) if ok else None
    return out


def matched_cost(rows: dict, budget: float, key: str = "flops") -> dict:
    """Per basis: most accurate entry with cost <= ``budget``; None if even the cheapest exceeds it."""
    out = {}
    for kind, entries in rows.items():
        ok = [e for e in entries if e[key] <= budget]
        out[kind] = min(ok, key=lambda e: e["err_med"]) if ok else None
    return out
