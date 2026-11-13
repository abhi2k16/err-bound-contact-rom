"""Step 7: acceptance criteria (Section 5.3) scored from the E3 audits and the saved E4 sweep."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .audit import run_basis_audit_sets
from .calibration import calibrate_all, pooled_summary
from .break_even_experiment import crossovers, criterion4_verdict
from .config import Config
from .dual_bases import DUAL_KINDS
from .metrics import assert_theorem1, criterion2_verdict, effectivity_stats


def score_benchmark(cfg: Config, problem, prepare) -> dict:
    """Criteria 1-3 on one benchmark: E3 audits on identical snapshots, ``tau`` calibrated on validation."""
    snaps, mu_test, mu_val = prepare(cfg, problem)
    audits = run_basis_audit_sets(problem, cfg, snaps, {"test": mu_test, "val": mu_val}, DUAL_KINDS,
                                  cfg.r_dual_list_e3, repair=True)
    test = audits["test"]
    assert_theorem1(test)
    pooled = pooled_summary(calibrate_all(audits["val"], test, cfg))
    rho = {kind: min(effectivity_stats(e, cfg)["spearman"] for (k, _), e in test.items() if k == kind)
           for kind in DUAL_KINDS}
    clam = {"c2_Clam": criterion2_verdict(test, cfg, "bound_Clam"),
            "c3_min_Clam": {kind: min(effectivity_stats(e, cfg, "bound_Clam")["spearman"]
                                      for (k, _), e in test.items() if k == kind) for kind in DUAL_KINDS},
            "c3_median_Clam": {kind: float(np.median([effectivity_stats(e, cfg, "bound_Clam")["spearman"]
                                                      for (k, _), e in test.items() if k == kind])) for kind in DUAL_KINDS},
            "valid_Clam": {kind: float(np.mean([effectivity_stats(e, cfg, "bound_Clam")["valid"]
                                                for (k, _), e in test.items() if k == kind])) for kind in DUAL_KINDS}}
    return {**clam, "c1": pooled, "c2": criterion2_verdict(test, cfg), "c3_min": rho,
            "c3_median": {kind: float(np.median([effectivity_stats(e, cfg)["spearman"]
                                                 for (k, _), e in test.items() if k == kind]))
                          for kind in DUAL_KINDS}}


def load_sweep(path: Path) -> list:
    """Rows of the saved E4 N sweep."""
    data = np.load(path)
    return [{key: float(data[key][i]) for key in data.files} for i in range(len(data["N"]))]


def criterion4(sweep_path: Path) -> dict:
    """Predicted vs measured break-even N* (membrane sweep) and the factor-2 verdict."""
    cross = crossovers(load_sweep(sweep_path))
    return {"cross": cross, "verdict": criterion4_verdict(cross)}


def print_summary(results: dict, c4_by_bench: dict) -> None:
    """Console summary table per benchmark and criterion."""
    print("\nAcceptance summary (Section 5.3)")
    for bench, res in results.items():
        print(f"\n[{bench}]")
        for kind, per in res["c1"].items():
            for ind, s in per.items():
                fa = np.mean(s["fa"]) if s["fa"] else float("nan")
                tag = "PASS" if s["missed"] == 0 and s["silent"] > 0 else ("n/a" if s["silent"] == 0 else "FAIL")
                print(f"  C1 {kind:<4}{ind:<11} missed {s['missed']}/{s['silent']}  false alarm {fa:.2f}  "
                      f"uncalibrated {s['uncal']}  {tag}")
        for kind, v in res["c2"].items():
            print(f"  C2 {kind:<4} theta in [1,10] on {v['share']:.0%} of {v['n_scope']} in-scope  "
                  f"{'PASS' if v['pass'] else 'FAIL'}")
        for kind, rho in res["c3_min"].items():
            rho_c = res["c3_min_Clam"][kind]
            print(f"  C3 {kind:<4} exact-lambda bound: Spearman min {rho:.3f}, median {res['c3_median'][kind]:.3f}  "
                  f"{'PASS' if rho > 0.90 else 'FAIL'} | C_lambda bound: min {rho_c:.3f}, median "
                  f"{res['c3_median_Clam'][kind]:.3f}  {'PASS' if rho_c > 0.90 else 'FAIL'}, valid on {res['valid_Clam'][kind]:.0%}")
    for bench, c4 in c4_by_bench.items():
        for regime, v in c4["cross"].items():
            verdict = c4["verdict"][regime]
            tag = "PASS" if verdict else "FAIL" if verdict is not None else "n/a (no crossover in swept range)"
            print(f"  C4 {bench}, {regime}: N* pred {v['N_pred']:.0f}, meas {v['N_meas']:.0f}, "
                  f"ratio {v['ratio']:.2f}  {tag}")
    print("  C4 B1: fixed size, no N sweep, not evaluated")
