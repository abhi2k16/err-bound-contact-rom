"""Calibration-study figure and table from ``results/calibration/analysis_{b1,b2}.pkl``.

Run ``python -m scripts.figures.make_calibration_figure`` after ``python -m contact_rom calibration --bench b1`` (and b2).
Writes ``figures/fig_calibration.pdf`` and ``tables/calibration_table.tex`` (a ``tabular`` fragment).
Intervals are 95% percentile bootstrap intervals over the eight seeds (the independent units).
"""

from __future__ import annotations

import pickle
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts.figures import make_figures as mf   # shared style and palette

RES = Path("results") / "calibration"
KINDS = ("nmf", "cpg", "gis", "rid")        # signed SVD is detected almost always and has no accurate queries
N_VALS = (10, 20, 40, 80, 200)
IND = "eta_adm"


def load(bench: str) -> dict:
    return pickle.load(open(RES / f"analysis_{bench}.pkl", "rb"))


def fig_calibration() -> None:
    fig, axes = plt.subplots(2, len(KINDS), figsize=(7.0, 3.3), sharex=True, sharey=True)
    for row, bench in enumerate(("b1", "b2")):
        st = load(bench)
        for col, kind in enumerate(KINDS):
            ax = axes[row, col]
            for rule, colour, ls, label in (("min09", mf.BLUE, "-", r"$0.9\,\min$ rule"),
                                            ("hybrid", mf.ORANGE, "--", "hybrid rule")):
                rec = [st[(kind, IND, rule, n)]["recall"] for n in N_VALS]
                lo = [st[(kind, IND, rule, n)]["boot"]["recall"][0] for n in N_VALS]
                hi = [st[(kind, IND, rule, n)]["boot"]["recall"][1] for n in N_VALS]
                ax.plot(N_VALS, rec, color=colour, ls=ls, lw=1.2, marker="o", ms=2.5, label=label)
                ax.fill_between(N_VALS, lo, hi, color=colour, alpha=0.15, lw=0)
            ax.set_xscale("log")
            ax.set_ylim(0.6, 1.005)
            if row == 0:
                ax.set_title(mf.LABEL[kind], fontsize=8)
            if col == 0:
                ax.set_ylabel(f"{'B1' if bench == 'b1' else 'B2'}: recall")
            if row == 1:
                ax.set_xlabel("validation queries")
    axes[0, 0].legend(loc="lower right", fontsize=6.5)
    fig.tight_layout()
    mf.save(fig, "fig_calibration")


def cell(stats: dict, key: str, which: str) -> str:
    lo, hi = stats["boot"][key]
    return f"{stats[which]:.3f} [{lo:.2f}, {hi:.2f}]" if not np.isnan(stats[which]) else "n/a"


def table(n_val: int) -> str:
    lines = []
    for bench, name in (("b1", "B1"), ("b2", "B2")):
        st = load(bench)
        for i, kind in enumerate(KINDS):
            a, h = st[(kind, IND, "min09", n_val)], st[(kind, IND, "hybrid", n_val)]
            lines.append(f"{name if i == 0 else ''} & {mf.LABEL[kind]} & {cell(a, 'recall', 'recall')} & "
                         f"{a['uncal_share']:.2f} & {cell(h, 'recall', 'recall')} & {cell(h, 'fpr', 'fpr')} \\\\")
        if bench == "b1":
            lines.append(r"\midrule")
    return "\n".join(lines)


def write_table() -> None:
    Path("tables").mkdir(exist_ok=True)
    body = ("\\begin{tabular}{@{}llcccc@{}}\n\\toprule\n & & \\multicolumn{2}{c}{$0.9\\,\\min$ rule} & "
            "\\multicolumn{2}{c}{hybrid rule} \\\\\n\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n"
            "Bench. & Basis & Recall & Uncal. & Recall & False-alarm rate \\\\\n\\midrule\n")
    for n_val in (40, 200):
        out = Path("tables") / f"calibration_table_{n_val}.tex"
        out.write_text(body + table(n_val) + "\n\\bottomrule\n\\end{tabular}\n", encoding="utf-8")
        print("wrote", out)


if __name__ == "__main__":
    fig_calibration()
    write_table()
