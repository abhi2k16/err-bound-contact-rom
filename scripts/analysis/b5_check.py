"""Criteria 1-3, detection recall and false-alarm rate on B5 (membrane with non-node-aligned, interpolated constraint rows).

Reads ``results/calibration/audits_b5_seed{42..49}.pkl`` written by ``python -m contact_rom calibration --bench b5`` and
applies the same scoring as the eight-seed check of B1 and B2. Writes ``results/analysis/b5_check.json``.
Run from the repository root: ``python -m scripts.analysis.b5_check``.
"""

import json
from pathlib import Path

import numpy as np

from scripts.analysis import detector_sensitivity as ds
from scripts.analysis import eight_seed_criteria as es


def main() -> None:
    data = es.load("b5")
    c23 = es.criteria23(data)
    out = {}
    for kind in es.KINDS:
        row = {}
        c = c23[kind]
        row["c2_share"] = float(np.nanmean([v["c2_share"] for v in c.values()]))
        row["c2_n"] = int(sum(v["n_scope"] for v in c.values()))
        row["c3_min_over_seeds"] = float(min(v["rho_min"] for v in c.values()))
        row["c3_median"] = float(np.median([v["rho_med"] for v in c.values()]))
        row["c3_seeds_pass"] = int(sum(v["rho_min"] > 0.9 for v in c.values()))
        row["valid_exact_min"] = float(min(v["valid_exact"] for v in c.values()))
        row["valid_C_min"] = float(min(v["valid_C"] for v in c.values()))
        for n_val in es.N_VAL_LIST:
            cells = es.criterion1(data, n_val, "eta_adm")[kind]
            row[f"c1_{n_val}_missed"] = int(sum(x["missed"] for x in cells.values()))
            row[f"c1_{n_val}_silent"] = int(sum(x["silent"] for x in cells.values()))
            row[f"c1_{n_val}_uncal"] = int(sum(x["uncal"] for x in cells.values()))
        cnt, ne, nd, ns = ds.score(data, ds.EPS, ds.TOL, (kind,), ("adm",))
        rec, fpr, npos, nneg = ds.rates(cnt["adm"])
        row.update({"recall": rec, "fpr": fpr, "n_silent": npos, "n_accurate": nneg, "entries_defined": f"{nd}/{ne}"})
        errs = [np.median(e["err_energy"]) for p in data.values() for (k, _), e in p["test"].items() if k == kind]
        row["median_err"] = float(np.median(errs))
        row["silent_share"] = float(np.mean([((e["res_red"] <= ds.TOL) & (e["err_energy"] > ds.EPS)).mean()
                                             for p in data.values() for (k, _), e in p["test"].items() if k == kind]))
        out[kind] = row
        print(kind, json.dumps(row))
    Path("results/analysis").mkdir(parents=True, exist_ok=True)
    Path("results/analysis/b5_check.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
