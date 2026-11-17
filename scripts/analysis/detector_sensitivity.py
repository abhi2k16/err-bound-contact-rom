"""Single-indicator baselines and sensitivity of the silent-failure definition (B1, B2; eight seeds).

Reads the audits saved by ``python -m contact_rom calibration`` (``results/calibration/audits_{bench}_seed{s}.pkl``;
five dual bases, ``r_lambda`` in {4, 8, 12, 16, 24}, 200 validation and 100 test queries per seed) and writes

* ``results/analysis/detectors.json``       recall and false-alarm rate of each indicator used alone,
* ``tables/detectors.tex``                  the corresponding table,
* ``results/analysis/sensitivity.json``    recall, false-alarm rate and counts for other tolerances,
* ``tables/sensitivity.tex``                the corresponding table.

Calibration follows the manuscript's rule: ``tau = 0.9 * min`` of the indicator over the validation silent failures
(first 40 validation queries); entries without a validation silent failure have no threshold and are left out for
every detector, so that all detectors are scored on the same entries. A silent failure has a reduced KKT residual at
or below ``tol`` and a relative energy error above ``eps``; an accurate query has an error at or below ``eps``.

Run from the repository root with ``python -m scripts.analysis.detector_sensitivity``.
"""

import json
import pickle
from pathlib import Path

import numpy as np

SEEDS = tuple(range(42, 50))
BENCHES = ("b1", "b2")
N_VAL = 40
EPS, TOL = 0.05, 1e-12
ZERO_TOL = 1e-10                    # "any violation" detector: a sign or gap violation above this value
BASES_NONNEG = ("nmf", "cpg", "gis", "rid")


def load(bench: str, seed: int) -> dict:
    with open(Path("results/calibration") / f"audits_{bench}_seed{seed}.pkl", "rb") as fh:
        return pickle.load(fh)


def indicators(e: dict) -> dict:
    return {
        "sign": e["eta_sign"], "gap": e["eta_gap"], "comp": np.sqrt(e["eta_comp"]),
        "adm": e["eta_adm"], "adm_v2": e["eta_adm_v2"],
    }


def score(bench_data: dict, eps: float, tol: float, bases: tuple, names: tuple):
    """Counts (tp, fn, fp, tn) per detector, over entries that have a validation silent failure."""
    counts = {n: np.zeros(4) for n in names + ("any_violation",)}
    n_entries = n_defined = n_silent_all = 0
    for d in bench_data.values():
        for key, val in d["val"].items():
            if key[0] not in bases:
                continue
            test = d["test"][key]
            v = {k: np.asarray(x)[:N_VAL] for k, x in val.items()}
            silent_val = (v["res_red"] <= tol) & (v["err_energy"] > eps)
            n_entries += 1
            n_silent_all += int(((test["res_red"] <= tol) & (test["err_energy"] > eps)).sum())
            if not silent_val.any():
                continue
            n_defined += 1
            pos = (test["res_red"] <= tol) & (test["err_energy"] > eps)
            neg = test["err_energy"] <= eps
            iv, it = indicators(v), indicators(test)
            for n in names:
                tau = 0.9 * iv[n][silent_val].min()
                flag = it[n] > tau
                counts[n] += [flag[pos].sum(), (~flag)[pos].sum(), flag[neg].sum(), (~flag)[neg].sum()]
            flag = (test["eta_sign"] > ZERO_TOL) | (test["eta_gap"] > ZERO_TOL)
            counts["any_violation"] += [flag[pos].sum(), (~flag)[pos].sum(), flag[neg].sum(), (~flag)[neg].sum()]
    return counts, n_entries, n_defined, n_silent_all


def f2(x: float) -> str:
    return "--" if x != x else f"{x:.2f}"


def rates(c: np.ndarray) -> tuple:
    tp, fn, fp, tn = c
    return (tp / (tp + fn) if tp + fn else float("nan"), fp / (fp + tn) if fp + tn else float("nan"), int(tp + fn), int(fp + tn))


def main() -> None:
    data = {b: {s: load(b, s) for s in SEEDS} for b in BENCHES}
    names = ("sign", "gap", "comp", "adm", "adm_v2")
    label = {"sign": r"$\eta_{\mathrm{sign}}$", "gap": r"$\eta_{\mathrm{gap}}$", "comp": r"$\sqrt{\eta_{\mathrm{comp}}}$",
             "adm": r"$\eta_{\mathrm{adm}}$ (sum)", "adm_v2": r"$\eta_{\mathrm{adm,v2}}$",
             "any_violation": r"any sign or gap violation $>10^{-10}$"}
    out = {"detectors": {}, "sensitivity": {}}

    # ---- detectors used alone --------------------------------------------------------------------------------
    rows = []
    for group, bases in (("signed SVD", ("svd",)), ("non-negative bases", BASES_NONNEG)):
        for n in names + ("any_violation",):
            cells = []
            for b in BENCHES:
                c, ne, nd, ns = score(data[b], EPS, TOL, bases, names)
                rec, fpr, npos, nneg = rates(c[n])
                out["detectors"][f"{group}|{n}|{b}"] = {"recall": rec, "fpr": fpr, "n_silent": npos, "n_accurate": nneg,
                                                         "entries": ne, "defined": nd}
                cells.append(f"{f2(rec)} & {f2(fpr)} & {npos}")
            rows.append((group, n, cells))
    lines = []
    last = None
    for group, n, cells in rows:
        if group != last and last is not None:
            lines.append(r"\midrule")
        lines.append(f"{group if group != last else ''} & {label[n]} & " + " & ".join(cells) + r" \\")
        last = group
    head = (r"\begin{tabular}{@{}llrrrrrr@{}}" "\n" r"\toprule" "\n"
            r" & & \multicolumn{3}{c}{B1} & \multicolumn{3}{c}{B2} \\" "\n"
            r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}" "\n"
            r"Bases & Detector & Recall & FPR & $n_{\mathrm{sf}}$ & Recall & FPR & $n_{\mathrm{sf}}$ \\" "\n" r"\midrule" "\n")
    Path("tables").mkdir(exist_ok=True)
    Path("tables/detectors.tex").write_text(head + "\n".join(lines) + "\n" r"\bottomrule" "\n" r"\end{tabular}" "\n", encoding="utf-8")

    # ---- sensitivity of the silent-failure definition ----------------------------------------------------------
    srows = []
    configs = [(0.01, 1e-12), (0.02, 1e-12), (0.05, 1e-12), (0.10, 1e-12), (0.20, 1e-12),
               (0.05, 1e-10), (0.05, 1e-14)]
    for eps, tol in configs:
        cells = []
        for b in BENCHES:
            c, ne, nd, ns = score(data[b], eps, tol, ("svd",) + BASES_NONNEG, ("adm",))
            rec, fpr, npos, nneg = rates(c["adm"])
            out["sensitivity"][f"{eps}|{tol}|{b}"] = {"recall": rec, "fpr": fpr, "n_silent_scored": npos, "n_silent_all": ns,
                                                       "entries": ne, "defined": nd}
            cells.append(f"{ns} & {nd}/{ne} & {f2(rec)} & {f2(fpr)}")
        srows.append(f"{100 * eps:g}\\% & $10^{{{int(round(np.log10(tol)))}}}$ & " + " & ".join(cells) + r" \\")
    head = (r"\begin{tabular}{@{}llrrrrrrrr@{}}" "\n" r"\toprule" "\n"
            r" & & \multicolumn{4}{c}{B1} & \multicolumn{4}{c}{B2} \\" "\n"
            r"\cmidrule(lr){3-6}\cmidrule(lr){7-10}" "\n"
            r"Error threshold & Residual tol. & $n_{\mathrm{sf}}$ & Thresholds & Recall & FPR & $n_{\mathrm{sf}}$ & Thresholds & Recall & FPR \\" "\n" r"\midrule" "\n")
    Path("tables/sensitivity.tex").write_text(head + "\n".join(srows) + "\n" r"\bottomrule" "\n" r"\end{tabular}" "\n", encoding="utf-8")

    Path("results/analysis").mkdir(parents=True, exist_ok=True)
    Path("results/analysis/detectors_sensitivity.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(Path("tables/detectors.tex").read_text())
    print(Path("tables/sensitivity.tex").read_text())


if __name__ == "__main__":
    main()
