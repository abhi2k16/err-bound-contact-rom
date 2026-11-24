"""Summary table of the independent 3D benchmark B4 -> ``manuscript_main/tables/b4_summary.tex``.

Needs ``results/b4_e3/e3_b4.npz`` (``python -m contact_rom e3 --bench b4``) and ``results/calibration/analysis_b4.pkl``
(``python -m contact_rom calibration --bench b4``). One row per basis: nominal dual dimension 16 (seed 42, 40 test
queries), then the pooled calibration results of eight seeds with 40 validation queries (hybrid rule, Eq. (eta-adm)).
"""
from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np

KINDS = ("svd", "nmf", "cpg", "gis", "rid")
LAB = {"svd": "SVD", "nmf": "NMF", "cpg": "CPG", "gis": "gIS", "rid": "RID"}


def main() -> None:
    d = np.load("results/b4_e3/e3_b4.npz")
    st = pickle.load(open("results/calibration/analysis_b4.pkl", "rb"))
    lines = []
    for kind in KINDS:
        g = lambda m: d[f"{kind}_16_{m}"]
        err = g("err_energy")
        share = 100 * float(np.mean(err > 0.05))
        r_eff, r_rep = int(np.median(g("r_dual_eff"))), int(np.median(g("r_rep")))
        bound_ok = 100 * float(np.mean(g("bound1") >= g("err_abs") * (1 - 1e-8)))
        a, v = st[(kind, "eta_adm", "hybrid", 40)], st[(kind, "eta_adm_v2", "hybrid", 40)]
        rec = lambda s: "n/a" if np.isnan(s["recall"]) else f"{s['recall']:.3f}"
        fpr = lambda s: "n/a" if np.isnan(s["fpr"]) else f"{s['fpr']:.2f}"
        lines.append(f"{LAB[kind]} & {r_eff} & {r_rep} & {np.median(err):.2g} & {share:.0f} & {bound_ok:.0f} & "
                     f"{rec(a)} / {fpr(a)} & {rec(v)} / {fpr(v)} \\\\")
    head = "\n".join([r"\begin{tabular}{@{}lrrrrrcc@{}}", r"\toprule",
                      r" & & & & & & \multicolumn{2}{c}{recall / false-alarm rate} \\",
                      r"\cmidrule(lr){7-8}",
                      r"Basis & $r_\lambda^{\mathrm{eff}}$ & $r_{\mathrm{rep}}$ & median $\varepsilon_E$ & $>5\%$ (\%) & bound valid (\%) & Eq.~\eqref{eq:eta-adm} & $\eta_{\mathrm{adm,v2}}$ \\",
                      r"\midrule"])
    Path("manuscript_main/tables/b4_summary.tex").write_text(head + "\n" + "\n".join(lines) + "\n" + r"\bottomrule" + "\n" + r"\end{tabular}" + "\n",
                                             encoding="utf-8")
    print(Path("manuscript_main/tables/b4_summary.tex").read_text())


def calibration_table() -> None:
    """Supplementary table: recall / false-alarm rate of the 0.9 min and hybrid rules for Eq. (eta-adm), 40 and 200 queries."""
    st = pickle.load(open("results/calibration/analysis_b4.pkl", "rb"))
    lines = []
    for n_val in (40, 200):
        for i, kind in enumerate(("nmf", "cpg", "gis", "rid")):
            a, h = st[("%s" % kind, "eta_adm", "min09", n_val)], st[(kind, "eta_adm", "hybrid", n_val)]
            rec = lambda s: "n/a" if np.isnan(s["recall"]) else f"{s['recall']:.3f}"
            fpr = lambda s: "n/a" if np.isnan(s["fpr"]) else f"{s['fpr']:.2f}"
            lines.append(f"{n_val if i == 0 else ''} & {LAB[kind]} & {rec(a)} & {a['uncal_share']:.2f} & {rec(h)} & {fpr(h)} \\\\")
        if n_val == 40:
            lines.append(r"\midrule")
    head = "\n".join([r"\begin{tabular}{@{}rlcccc@{}}", r"\toprule",
                      r" & & \multicolumn{2}{c}{$0.9\,\min$ rule} & \multicolumn{2}{c}{hybrid rule} \\",
                      r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}",
                      r"Val. queries & Basis & Recall & Uncal. & Recall & False-alarm rate \\", r"\midrule"])
    Path("manuscript_main/tables/b4_calibration.tex").write_text(head + "\n" + "\n".join(lines) + "\n" + r"\bottomrule" + "\n" + r"\end{tabular}" + "\n",
                                                 encoding="utf-8")
    print(Path("manuscript_main/tables/b4_calibration.tex").read_text())


if __name__ == "__main__":
    main()
    calibration_table()
