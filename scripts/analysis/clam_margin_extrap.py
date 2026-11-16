"""Out-of-range test of a safety margin on the offline multiplier-norm constant ``C_lambda``.

Repeats the out-of-range audits of ``contact_rom.scale_study`` (B1, B2; five seeds; r_lambda = 16; rank repair) and, for
each query set and dual basis, evaluates the online bound with the constant scaled by ``m``:
``bound_m^2 = bound_Clam^2 + 2 (m - 1) C eta_gap``. A query needs ``m >= ||lambda|| / C`` for the bound to hold; this
required factor is recovered from the exact-norm bound, ``||lambda|| = C + (bound1^2 - bound_Clam^2) / (2 eta_gap)``.

Run: ``python -m scripts.analysis.clam_margin_extrap [b1|b2]``. Writes ``results/analysis/clam_margin_{bench}.json``.
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

from contact_rom.audit import run_basis_audit_sets
from contact_rom.bases import build_snapshots
from contact_rom.config import Config
from contact_rom.dual_bases import DUAL_KINDS
from contact_rom.run import BENCHMARKS
from contact_rom.scale_study import R_DUAL, SEEDS, extrap_sets

MARGINS = (1.0, 1.1, 1.25, 1.5, 2.0, 3.0)
REL_TOL = 1e-8


def main(bench: str) -> None:
    cfg = Config()
    acc: dict = {}
    for seed in SEEDS:
        run_cfg = replace(cfg, seed=seed)
        problem = BENCHMARKS[bench](run_cfg)
        rng = np.random.default_rng(seed + 1)
        snaps = build_snapshots(problem, problem.sample_mu(rng, run_cfg.n_train), run_cfg.lcp_method)
        sets = {"test": problem.sample_mu(rng, run_cfg.n_test), "val": problem.sample_mu(rng, run_cfg.n_val)}
        sets.update(extrap_sets(run_cfg, bench, rng, 40))
        audits = run_basis_audit_sets(problem, run_cfg, snaps, sets, DUAL_KINDS, (R_DUAL,), repair=True)
        C = float(np.linalg.norm(snaps.Lam_snap, axis=0).max())
        print(f"{bench} seed {seed}: C_lambda = {C:.4g}", flush=True)
        for name in sets:
            if name == "val":
                continue
            for (kind, _), e in audits[name].items():
                slot = acc.setdefault((kind, name), {"n": 0, "need": [], "valid": {str(m): 0 for m in MARGINS},
                                                     "gap0": 0, "theta": {str(m): [] for m in MARGINS},
                                                     "strict_invalid": 0, "short": []})
                err = e["err_abs"]
                n = err.size
                slot["n"] += n
                gap = e["eta_gap"]
                with np.errstate(divide="ignore", invalid="ignore"):
                    need = np.where(gap > 0, 1.0 + (e["bound1"] ** 2 - e["bound_Clam"] ** 2) / (2.0 * C * gap), 1.0)
                slot["need"].append(need)
                slot["gap0"] += int((gap <= 0).sum())
                strict_bad = e["bound_Clam"] < err
                slot["strict_invalid"] += int(strict_bad.sum())
                slot["short"].append((1.0 - e["bound_Clam"] / np.maximum(err, 1e-300))[strict_bad])
                scope = e["err_energy"] > run_cfg.err_silent
                for m in MARGINS:
                    b = np.sqrt(np.maximum(e["bound_Clam"] ** 2 + 2.0 * (m - 1.0) * C * gap, 0.0))
                    slot["valid"][str(m)] += int((b >= err * (1 - REL_TOL)).sum())
                    slot["theta"][str(m)].append((b / np.maximum(err, 1e-300))[scope])
    rows = {}
    for (kind, name), s in acc.items():
        need = np.concatenate(s["need"])
        rows[f"{kind}|{name}"] = {
            "n": s["n"], "gap_zero": s["gap0"], "strict_invalid": s["strict_invalid"],
            "shortfall_max": float(np.concatenate(s["short"]).max()) if s["strict_invalid"] else 0.0,
            "n_above_tol": int((np.concatenate(s["short"]) > REL_TOL).sum()) if s["strict_invalid"] else 0, "need_max": float(need.max()),
            "need_p99": float(np.percentile(need, 99)),
            "valid": {m: v / s["n"] for m, v in s["valid"].items()},
            "n_invalid": {m: s["n"] - v for m, v in s["valid"].items()},
            "theta_med": {m: (float(np.median(np.concatenate(t))) if sum(len(x) for x in t) else float("nan"))
                          for m, t in s["theta"].items()},
        }
    out = Path("results/analysis")
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"clam_margin_{bench}.json"
    path.write_text(json.dumps({"bench": bench, "seeds": list(SEEDS), "margins": MARGINS, "rows": rows}, indent=1), encoding="utf-8")
    print("saved", path)
    print("\nbasis|set  n  need_max  valid% at m = " + ", ".join(str(m) for m in MARGINS) + "   median theta (>5% queries) at m=1, 2")
    for key, r in rows.items():
        v = " ".join(f"{100 * r['valid'][str(m)]:6.1f}" for m in MARGINS)
        print(f"{key:<22} {r['n']:>4}  {r['need_max']:6.2f}  {v}   {r['theta_med']['1.0']:.2f}, {r['theta_med']['2.0']:.2f}")
        if r["strict_invalid"]:
            print(f"    strict m=1 violations: {r['strict_invalid']}/{r['n']}, of which {r['n_above_tol']} exceed 1e-8 relative; "
                  f"largest relative shortfall {r['shortfall_max']:.3g}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "b1")
