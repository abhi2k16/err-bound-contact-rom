"""Compare our simplified gIS (top rows by snapshot energy) with Algorithm 2 of Sambataro and Ehrlacher (``gis2``).

Same protocol as the calibration study (B1, B2; seeds 42-49; 200 validation and 100 test queries per seed;
r_lambda in {4, 8, 12, 16, 24}; rank repair). Detection scoring uses the manuscript rule (threshold 0.9 min over the
first 40 validation silent failures). Writes results/analysis/gis_variant.json and manuscript_main/tables/gis_variant.tex.

Run from the repository root: ``python -m scripts.analysis.gis_variant``.
"""

import json
import pickle
from dataclasses import replace
from pathlib import Path

import numpy as np

from contact_rom.audit import run_basis_audit_sets
from contact_rom.calibration_study import _slim
from contact_rom.config import Config
from contact_rom.estimator import tau_adm_from_snapshots
from contact_rom.run import BENCHMARKS, prepare
from scripts.analysis import detector_sensitivity as ds

SEEDS = tuple(range(42, 50))
KINDS = ("gis", "gis2")
OUT = Path("results/gis_variant")


def generate(bench: str) -> dict:
    cfg = Config()
    OUT.mkdir(parents=True, exist_ok=True)
    data = {}
    for seed in SEEDS:
        path = OUT / f"audits_{bench}_seed{seed}.pkl"
        if not path.exists():
            run_cfg = replace(cfg, seed=seed, n_val=200, n_test=100)
            problem = BENCHMARKS[bench](run_cfg)
            print("gis_variant", bench, seed, flush=True)
            snaps, mu_test, mu_val = prepare(run_cfg, problem)
            audits = run_basis_audit_sets(problem, run_cfg, snaps, {"test": mu_test, "val": mu_val}, KINDS,
                                          run_cfg.r_dual_list_e3, repair=True)
            payload = {"val": _slim(audits["val"]), "test": _slim(audits["test"])}
            path.write_bytes(pickle.dumps(payload))
        data[seed] = pickle.loads(path.read_bytes())
    return data


def main() -> None:
    res = {}
    for bench in ("b1", "b2"):
        data = generate(bench)
        for kind in KINDS:
            errs, eff, rep, nsf, nacc, ntest = [], [], [], 0, 0, 0
            for d in data.values():
                for key, test in d["test"].items():
                    if key[0] != kind:
                        continue
                    errs.append(np.median(test["err_energy"]))
                    eff.append(float(np.median(test["r_dual_eff"]))); rep.append(float(np.median(test["r_rep"])))
                    sf = (test["res_red"] <= ds.TOL) & (test["err_energy"] > ds.EPS)
                    nsf += int(sf.sum()); nacc += int((test["err_energy"] <= ds.EPS).sum()); ntest += len(test["err_energy"])
            c, ne, nd, ns = ds.score(data, ds.EPS, ds.TOL, (kind,), ("adm",))
            rec, fpr, npos, nneg = ds.rates(c["adm"])
            res[f"{bench}|{kind}"] = {"median_err": float(np.median(errs)), "silent_share": nsf / ntest, "accurate_share": nacc / ntest,
                                      "eff_dim_median": float(np.median(eff)), "r_rep_median": float(np.median(rep)),
                                      "recall": rec, "fpr": fpr, "n_silent_scored": npos, "thresholds": f"{nd}/{ne}"}
            print(bench, kind, res[f"{bench}|{kind}"], flush=True)
    Path("results/analysis").mkdir(parents=True, exist_ok=True)
    Path("results/analysis/gis_variant.json").write_text(json.dumps(res, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
