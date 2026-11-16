"""Figure and tables for the matched dual-basis comparison from ``results/basis_comparison/*.json``.

Run ``python -m scripts.figures.make_basis_comparison`` after ``python -m contact_rom basiscomp --bench b1`` (and b2).
Writes ``figures/fig_basis_cost.pdf``, ``tables/basis_matched_dim.tex`` and ``tables/basis_matched_acc.tex``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.figures import make_figures as mf
from contact_rom import basis_comparison as bc
from contact_rom.config import Config

CFG = Config()
KINDS = ("svd", "nmf", "cpg", "gis", "rid")
STYLE = {"svd": (mf.MUTED, "x"), "nmf": (mf.BLUE, "o"), "cpg": (mf.ORANGE, "s"), "gis": (mf.AQUA, "^"), "rid": (mf.INK, "D")}
DIM_TARGETS = {"b1": (16, 24), "b2": (16, 60)}
ACC_TARGETS = (0.05, 0.01)


def fig_cost_accuracy() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7))
    for ax, bench, title in zip(axes, ("b1", "b2"), ("B1", "B2")):
        rows = bc.load(CFG, bench)["rows"]
        for kind in KINDS:
            colour, marker = STYLE[kind]
            ax.plot([e["flops"] for e in rows[kind]], [e["err_med"] for e in rows[kind]], color=colour, marker=marker,
                    ms=3, lw=1.0, label=mf.LABEL[kind])
        ax.axhline(mf.TOL_ERR, color=mf.MUTED, lw=0.6, ls=":")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("online FLOPs per query, Eq. (17)")
        ax.set_title(title, fontsize=8)
    axes[0].set_ylabel("median relative energy error")
    axes[0].legend(fontsize=6.5, ncol=2, loc="lower left")
    fig.tight_layout()
    mf.save(fig, "fig_basis_cost")


def fmt_err(e) -> str:
    return "n/a" if e is None else f"{e['err_med']:.2g}"


def write_tables() -> None:
    Path("tables").mkdir(exist_ok=True)
    # matched effective dimension: median error, with the effective dual / primal dimensions actually used
    lines = []
    for bench, name in (("b1", "B1"), ("b2", "B2")):
        rows = bc.load(CFG, bench)["rows"]
        for i, target in enumerate(DIM_TARGETS[bench]):
            matched = bc.matched_dimension(rows, target)
            cells = []
            for kind in KINDS:
                e = matched[kind]
                cells.append("n/a" if e is None else f"{e['err_med']:.2g} ({e['r_eff']}, {e['r_u']})")
            lines.append(f"{name if i == 0 else ''} & {target} & " + " & ".join(cells) + r" \\")
        if bench == "b1":
            lines.append(r"\midrule")
    head = (r"\begin{tabular}{@{}llccccc@{}}" "\n" r"\toprule" "\n"
            r"Bench. & $r_{\lambda}^{\mathrm{eff}}$ target & SVD & NMF & CPG & gIS & RID \\" "\n" r"\midrule" "\n")
    Path("tables/basis_matched_dim.tex").write_text(head + "\n".join(lines) + "\n" r"\bottomrule" "\n" r"\end{tabular}" "\n",
                                                    encoding="utf-8")
    # matched accuracy: cheapest entry reaching the target median error
    lines = []
    for bench, name in (("b1", "B1"), ("b2", "B2")):
        rows = bc.load(CFG, bench)["rows"]
        for i, tgt in enumerate(ACC_TARGETS):
            matched = bc.matched_accuracy(rows, tgt)
            cells = []
            for kind in KINDS:
                e = matched[kind]
                cells.append("not reached" if e is None else
                             f"{e['r_eff']}/{e['r_u']}, {e['flops'] / 1e3:.0f}k, {e['t_online'] * 1e6:.0f}")
            lines.append(f"{name if i == 0 else ''} & {int(tgt * 100)}\\% & " + " & ".join(cells) + r" \\")
        if bench == "b1":
            lines.append(r"\midrule")
    head = (r"\begin{tabular}{@{}llccccc@{}}" "\n" r"\toprule" "\n"
            r"Bench. & Median error & SVD & NMF & CPG & gIS & RID \\" "\n" r"\midrule" "\n")
    lines.append(r"\multicolumn{7}{@{}l}{\scriptsize Entries: effective dual / enriched primal dimension, online kFLOPs per query, solve time in $\mu$s.} \\")
    lines.append(r"\multicolumn{7}{@{}l}{\scriptsize ``Not reached'': no entry up to $r_\lambda=30$ (B1) or $60$ (B2) attains the target.} \\")
    Path("tables/basis_matched_acc.tex").write_text(head + "\n".join(lines) + "\n" r"\bottomrule" "\n" r"\end{tabular}" "\n",
                                                    encoding="utf-8")
    print("wrote tables")


if __name__ == "__main__":
    fig_cost_accuracy()
    write_tables()
