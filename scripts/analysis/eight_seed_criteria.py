"""Acceptance criteria 1-3 and the C_lambda margin study over eight seeds, from the saved calibration audits.

Reads ``results/calibration/audits_{bench}_seed{42..49}.pkl`` (200-query validation pool, 100 test queries per seed,
all five dual bases, r_lambda in {4, 8, 12, 16, 24}, rank repair). No full-order solve is repeated.

Criterion 1: threshold ``0.9 * min`` of the indicator over the first ``N_VAL`` validation queries, scored on the test
queries; counts are pooled over r_lambda. Criterion 2: share of queries above the 5% tolerance with theta in [1, 10].
Criterion 3: minimum (and median) Spearman correlation between the bound and the true error over r_lambda.
Margin study: bound with the offline constant scaled, ``C -> m C``, recomputed from the stored ``bound_Clam`` and
``eta_gap`` (``bound_m^2 = bound_Clam^2 + 2 (m - 1) C eta_gap``).

Run: ``python -m scripts.analysis.eight_seed_criteria``. Writes ``results/analysis/eight_seed_criteria.json``.
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from contact_rom.calibration_study import threshold
from contact_rom.config import Config
from contact_rom.metrics import silent_failure_mask

SEEDS = tuple(range(42, 50))
KINDS = ("svd", "nmf", "cpg", "gis", "rid")
N_VAL_LIST = (40, 200)
MARGINS = (1.0, 1.1, 1.25, 1.5, 2.0)
CFG = Config()


def load(bench: str) -> dict:
    return {s: pickle.load(open(f"results/calibration/audits_{bench}_seed{s}.pkl", "rb")) for s in SEEDS}


def subset(entry: dict, n: int) -> dict:
    pool = len(entry["err_energy"])
    return {k: v[:n] for k, v in entry.items() if np.ndim(v) == 1 and len(v) == pool}


def criterion1(data: dict, n_val: int, indicator: str) -> dict:
    """Per basis kind: per-seed (missed, silent, uncalibrated entries) for the 0.9 min rule."""
    out: dict = {k: {} for k in KINDS}
    for s, p in data.items():
        for kind in KINDS:
            missed = silent = uncal = 0
            for (k, r), val in p["val"].items():
                if k != kind:
                    continue
                tau = threshold("min09", subset(val, n_val), indicator, CFG, p["tau_rel"])
                test = p["test"][(k, r)]
                sil = silent_failure_mask(test, CFG)
                if np.isnan(tau):
                    uncal += 1
                    continue
                silent += int(sil.sum())
                missed += int((sil & ~(test[indicator] > tau)).sum())
            out[kind][s] = {"missed": missed, "silent": silent, "uncal": uncal}
    return out


def theta(entry: dict, key: str) -> np.ndarray:
    return entry[key] / np.maximum(entry["err_abs"], 1e-300)


def criteria23(data: dict) -> dict:
    out: dict = {k: {} for k in KINDS}
    for s, p in data.items():
        for kind in KINDS:
            entries = [e for (k, _), e in p["test"].items() if k == kind]
            rho = [float(spearmanr(e["bound1"], e["err_abs"])[0]) for e in entries]
            rho_c = [float(spearmanr(e["bound_Clam"], e["err_abs"])[0]) for e in entries]
            th = np.concatenate([theta(e, "bound1") for e in entries])
            scope = np.concatenate([e["err_energy"] > CFG.err_silent for e in entries])
            inr = (th >= 1) & (th <= 10)
            valid = float(np.mean(np.concatenate([e["bound1"] >= e["err_abs"] * (1 - 1e-8) for e in entries])))
            valid_c = float(np.mean(np.concatenate([e["bound_Clam"] >= e["err_abs"] for e in entries])))
            out[kind][s] = {"rho_min": float(np.nanmin(rho)), "rho_med": float(np.nanmedian(rho)),
                            "rho_min_C": float(np.nanmin(rho_c)), "rho_med_C": float(np.nanmedian(rho_c)),
                            "n_scope": int(scope.sum()), "c2_share": float(inr[scope].mean()) if scope.any() else float("nan"),
                            "valid_exact": valid, "valid_C": valid_c}
    return out


def margin_study(data: dict) -> dict:
    """Validity and looseness of the online bound for scaled constants, pooled over seeds and dual dimensions."""
    out: dict = {k: {} for k in KINDS}
    for kind in KINDS:
        for m in MARGINS:
            valid, th_scope, per_entry_valid = [], [], []
            for p in data.values():
                C = p["C_lambda"]
                for (k, _), e in p["test"].items():
                    if k != kind:
                        continue
                    b2 = e["bound_Clam"] ** 2 + 2.0 * (m - 1.0) * C * e["eta_gap"]
                    b = np.sqrt(np.maximum(b2, 0.0))
                    ok = b >= e["err_abs"] * (1 - 1e-8)
                    valid.append(ok)
                    per_entry_valid.append(float(ok.mean()))
                    scope = e["err_energy"] > CFG.err_silent
                    th_scope.append((b / np.maximum(e["err_abs"], 1e-300))[scope])
            v = np.concatenate(valid)
            ts = np.concatenate(th_scope) if th_scope else np.array([])
            out[kind][str(m)] = {"valid": float(v.mean()), "n_invalid": int((~v).sum()), "n": int(v.size),
                                 "worst_entry_valid": float(min(per_entry_valid)),
                                 "theta_med_scope": float(np.median(ts)) if ts.size else float("nan"),
                                 "theta_p90_scope": float(np.percentile(ts, 90)) if ts.size else float("nan")}
    return out


def main() -> None:
    result: dict = {}
    for bench in ("b1", "b2"):
        data = load(bench)
        res: dict = {"c23": criteria23(data), "margin": margin_study(data)}
        for n_val in N_VAL_LIST:
            for ind in ("eta_adm", "eta_adm_v2"):
                res[f"c1_{ind}_{n_val}"] = criterion1(data, n_val, ind)
        result[bench] = res
        print(f"\n===== {bench.upper()} (seeds {SEEDS[0]}-{SEEDS[-1]}) =====")
        for n_val in N_VAL_LIST:
            for ind in ("eta_adm", "eta_adm_v2"):
                print(f"\nCriterion 1, {ind}, {n_val} validation queries: missed/silent pooled | zero-miss seeds / seeds with silent | uncal entries")
                for kind in KINDS:
                    cells = res[f"c1_{ind}_{n_val}"][kind]
                    miss = sum(c["missed"] for c in cells.values())
                    sil = sum(c["silent"] for c in cells.values())
                    with_sil = [c for c in cells.values() if c["silent"] > 0]
                    zero = sum(c["missed"] == 0 for c in with_sil)
                    unc = sum(c["uncal"] for c in cells.values())
                    print(f"  {kind:<4} {miss}/{sil}  | {zero}/{len(with_sil)} | {unc}")
        print("\nCriterion 2 (pooled share theta in [1,10] among >5% queries) and criterion 3 (min Spearman over r, per seed)")
        for kind in KINDS:
            c = res["c23"][kind]
            n_scope = sum(v["n_scope"] for v in c.values())
            shares = [v["c2_share"] for v in c.values() if not np.isnan(v["c2_share"])]
            mins = [v["rho_min"] for v in c.values()]
            minsC = [v["rho_min_C"] for v in c.values()]
            print(f"  {kind:<4} C2 share {np.min(shares) if shares else float('nan'):.2f}-{np.max(shares) if shares else float('nan'):.2f} ({n_scope} queries) | "
                  f"C3 exact min {np.min(mins):.3f}..{np.max(mins):.3f}, pass {sum(m > 0.9 for m in mins)}/8, "
                  f"median of medians {np.median([v['rho_med'] for v in c.values()]):.3f} | "
                  f"C3 with C_lambda pass {sum(m > 0.9 for m in minsC)}/8 | valid exact {min(v['valid_exact'] for v in c.values()):.3f}, "
                  f"valid C {min(v['valid_C'] for v in c.values()):.3f}")
        print("\nMargin study (in-range test queries, pooled over 8 seeds and 5 dual dimensions):")
        print("  basis  m      valid%   invalid/n    worst-entry  theta median / p90 (queries above 5%)")
        for kind in KINDS:
            for m in MARGINS:
                r = res["margin"][kind][str(m)]
                print(f"  {kind:<5} {m:<5} {100 * r['valid']:7.2f}  {r['n_invalid']:>4}/{r['n']:<6} {100 * r['worst_entry_valid']:8.1f}%   "
                      f"{r['theta_med_scope']:.2f} / {r['theta_p90_scope']:.2f}")
    out = Path("results/analysis")
    out.mkdir(parents=True, exist_ok=True)
    (out / "eight_seed_criteria.json").write_text(json.dumps(result, indent=1, default=float), encoding="utf-8")
    print("\nsaved", out / "eight_seed_criteria.json")


if __name__ == "__main__":
    main()
