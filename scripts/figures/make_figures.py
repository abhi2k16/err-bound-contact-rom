"""Publication figures for the manuscript, generated from the saved results in ``results/``.

Run ``python -m scripts.figures.make_figures``; vector PDFs are written to ``figures/``. Palette: first three slots of the
validated categorical palette (blue, orange, aqua); multi-basis plots use small multiples, not extra hues.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RES, OUT = Path("results"), Path("figures")
BLUE, ORANGE, AQUA, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e"
KINDS = ("svd", "nmf", "cpg", "gis", "rid")
LABEL = {"svd": "SVD", "nmf": "NMF", "cpg": "CPG", "gis": "gIS", "rid": "RID"}
TOL_ERR, TOL_KKT = 0.05, 1e-12

plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "font.size": 8,
                     "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "legend.frameon": False,
                     "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.4, "axes.axisbelow": True,
                     "mathtext.fontset": "stix", "pdf.fonttype": 42})


def pooled(path: Path, kind: str, metric: str) -> np.ndarray:
    """Concatenate ``metric`` over all ``r_lambda`` entries of ``kind`` in an audit ``.npz``."""
    data = np.load(path)
    keys = sorted(k for k in data.files if k.startswith(f"{kind}_") and k.endswith(f"_{metric}")
                  and k[len(kind) + 1:-len(metric) - 1].isdigit())
    return np.concatenate([data[k] for k in keys])


def save(fig, name: str) -> None:
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    if os.environ.get("FIG_PREVIEW_DIR"):   # optional raster preview for visual checks
        fig.savefig(Path(os.environ["FIG_PREVIEW_DIR"]) / f"{name}.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / f"{name}.pdf")


def fig_silent_failure() -> None:
    """Reduced KKT residual vs true energy error: SVD sits in the silent-failure region, NMF mostly does not."""
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.7), sharey=True)
    for ax, (title, path) in zip(axes, (("B1 spring chain", RES / "b1_e1" / "e1_b1.npz"),
                                        ("B2 membrane", next((p for p in (RES / "b2_all" / "e1_b2.npz", RES / "b2_e1" / "e1_b2.npz")
                                                              if p.exists()), RES / "b2_e1" / "e1_b2.npz")))):
        ax.fill_between([1e-17, TOL_KKT], TOL_ERR, 10, color=ORANGE, alpha=0.08, linewidths=0)
        for kind, color, marker in (("svd", BLUE, "o"), ("nmf", ORANGE, "^")):
            res = np.maximum(pooled(path, kind, "res_red"), 1e-17) * 10 ** rng.uniform(-0.15, 0.15, 200)
            err = np.maximum(pooled(path, kind, "err_energy"), 1e-6)
            ax.scatter(res, err, s=9, marker=marker, color=color, alpha=0.75, linewidths=0, label=LABEL[kind])
        ax.axhline(TOL_ERR, color=MUTED, lw=0.7, ls="--")
        ax.axvline(TOL_KKT, color=MUTED, lw=0.7, ls=":")
        ax.set(xscale="log", yscale="log", xlim=(1e-17, 1e-9), ylim=(1e-4, 5), title=title)
        ax.set_xlabel(r"reduced KKT residual $\|R_{\mathrm{red}}\|_\infty$")
        ax.text(TOL_KKT * 0.6, 3.2, "silent-failure\nregion", fontsize=6.5, color=MUTED, ha="right", va="top")
        ax.text(1.5e-9, TOL_ERR * 1.15, r"5% error", fontsize=6.5, color=MUTED, ha="right")
        ax.text(TOL_KKT * 1.3, 1.2e-4, r"$10^{-12}$", fontsize=6.5, color=MUTED, ha="left")
    axes[0].set_ylabel("relative energy-norm error")
    axes[0].legend(loc="lower right", markerscale=1.5)
    save(fig, "fig_silent_failure")


def fig_bound() -> None:
    """Certified bound vs true energy error, one panel per basis (pooled over r_lambda), B1 and B2."""
    fig, axes = plt.subplots(2, 5, figsize=(6.3, 3.3), sharex="row", sharey="row")
    for row, (bench, path) in enumerate((("B1", RES / "b1_e3" / "e3_b1.npz"), ("B2", RES / "b2_e3" / "e3_b2.npz"))):
        data = {k: (np.maximum(pooled(path, k, "err_abs"), 1e-16), np.maximum(pooled(path, k, "bound1"), 1e-16))
                for k in KINDS}
        both = np.concatenate([np.concatenate(v) for v in data.values()])
        lim = [both.min() / 2, both.max() * 2]          # one common range per row
        for col, kind in enumerate(KINDS):
            ax = axes[row, col]
            err, bound = data[kind]
            ax.fill_between(lim, lim, np.array(lim) * 10, color=BLUE, alpha=0.10, linewidths=0)
            ax.plot(lim, lim, color=MUTED, lw=0.6)
            ax.scatter(err, bound, s=3, color=BLUE, alpha=0.5, linewidths=0)
            ax.set(xscale="log", yscale="log", xlim=lim, ylim=lim)
            ax.set_title(f"{LABEL[kind]}" if row == 0 else "", fontsize=8)
            ax.tick_params(labelsize=6.5)
            if col == 0:
                ax.set_ylabel(f"{bench}: bound")
            if row == 1:
                ax.set_xlabel("true error", fontsize=7)
    fig.tight_layout()
    save(fig, "fig_bound")


def load_sweep(path: Path) -> dict:
    data = np.load(path)
    return {k: data[k] for k in data.files}


def fig_breakeven() -> None:
    """Predicted vs measured speedup against problem size for both stiffness regimes, B2 and B3."""
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.0), sharey=True)
    for ax, (title, path) in zip(axes, (("B2 membrane", RES / "b2_e4" / "e4_sweep.npz"),
                                        ("B3 two-body", RES / "b3_e4" / "e4_sweep_b3.npz"))):
        d = load_sweep(path)
        for regime, color, name in (("aff", BLUE, "constant $K$"), ("var", ORANGE, r"varying $K(\mu)$")):
            ax.plot(d["N"], d[f"S_pred_{regime}"], color=color, lw=1.0, ls="--")
            ax.plot(d["N"], d[f"S_meas_{regime}"], color=color, lw=1.2, marker="o", ms=3, label=name)
        ax.axhline(1, color=MUTED, lw=0.7)
        ax.set(xscale="log", yscale="log", title=title, xlabel="full-order size $N$")
    axes[0].set_ylabel(r"speedup $S = \mathcal{C}_{FOM}/\mathcal{C}_{ROM}$")
    axes[0].legend(loc="upper left", title="solid: measured, dashed: predicted", title_fontsize=6.5, fontsize=7)
    save(fig, "fig_breakeven")


def fig_breaker() -> None:
    """E5 admissibility screening (B3): fallback share, speedup and accepted-query error, in range vs extrapolated.

    Filled marks: r_lambda = 16; open markers: r_lambda = 8. NMF and CPG come from the 100-query run (``b3_e5_n100``);
    SVD, gIS and RID, which fall back on every query, from the 100-query run ``b3_e5_n100_rest``.
    """
    d25 = np.load(RES / "b3_e5_n100_rest" / "e5_circuit_breaker.npz")
    d100 = np.load(RES / "b3_e5_n100" / "e5_circuit_breaker.npz")
    x = np.arange(len(KINDS))
    fig, axes = plt.subplots(1, 3, figsize=(6.9, 2.4))
    regimes = (("test", "in range", BLUE, -0.18), ("extrap", "extrapolated", ORANGE, 0.18))

    def get(k, r, name, key):
        return float((d100 if k in ("nmf", "cpg") else d25)[f"{k}_{r}_eta_adm_{name}_{key}"])

    for name, _, color, off in regimes:
        # (a) fallback share: bars at r=16, open marker at r=8
        v16 = [100 * get(k, 16, name, "fallback") for k in KINDS]
        v8 = [100 * get(k, 8, name, "fallback") for k in KINDS]
        axes[0].bar(x + off, v16, 0.34, color=color)
        axes[0].scatter(x + off, v8, s=16, facecolor="white", edgecolor=INK, linewidths=0.7, zorder=4)
        for xi, v in zip(x + off, v16):
            if v < 99:   # bars at 100 % are read off the axis; labels there collide
                axes[0].text(xi, v + 6, f"{v:.0f}", ha="center", va="bottom", fontsize=6, color=INK)
        # (b) speedup
        s16 = [get(k, 16, name, "speedup") for k in KINDS]
        s8 = [get(k, 8, name, "speedup") for k in KINDS]
        axes[1].scatter(x + off, s8, s=30, facecolor="white", edgecolor=color, linewidths=0.9, zorder=3)
        axes[1].scatter(x + off, s16, s=12, color=color, zorder=4)
        for xi, v in zip(x + off, s16):
            if v > 1.5:
                axes[1].text(xi, v * 1.2, f"{v:.1f}x", ha="center", va="bottom", fontsize=6.5, color=INK)
        # (c) largest true error among accepted queries (blank where nothing was accepted)
        for r, filled in ((16, True), (8, False)):
            e = np.array([100 * get(k, r, name, "max_err_accepted") for k in KINDS])
            ok = ~np.isnan(e)
            axes[2].scatter((x + off)[ok], e[ok], s=12 if filled else 30, zorder=4 if filled else 3,
                            facecolor=color if filled else "white", edgecolor=color, linewidths=0.9)

    axes[0].set(ylim=(0, 118), ylabel="queries sent to full-order solver (%)", title="(a) fallback share")
    axes[1].set(yscale="log", ylim=(0.6, 45), ylabel=r"net speedup $S_{net}$", title="(b) speedup")
    axes[1].axhline(1, color=MUTED, lw=0.7)
    axes[2].set(yscale="log", ylim=(1, 20), ylabel="max true error of accepted queries (%)", title="(c) accepted-query error")
    axes[2].axhline(100 * TOL_ERR, color=MUTED, lw=0.7, ls="--")
    axes[2].text(4.45, 100 * TOL_ERR * 1.06, "5% tolerance", fontsize=6.5, color=MUTED, ha="right", va="bottom")
    axes[2].text(4.45, 1.12, "no mark: all queries re-solved", fontsize=6, color=MUTED, ha="right", va="bottom")
    axes[2].set_yticks([1, 2, 5, 10, 20], ["1", "2", "5", "10", "20"])
    for ax in axes:
        ax.set_xticks(x, [LABEL[k] for k in KINDS], fontsize=7)
        ax.set_xlabel("dual basis", fontsize=7)
        ax.grid(axis="x", visible=False)
        ax.title.set_fontsize(8)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c, _ in regimes]
    handles += [plt.Line2D([], [], marker="o", ls="", mfc=INK, mec=INK, ms=4),
                plt.Line2D([], [], marker="o", ls="", mfc="white", mec=INK, ms=4)]
    fig.legend(handles, ["in range", "extrapolated", r"$r_\lambda=16$", r"$r_\lambda=8$"], loc="upper center", ncol=4,
               fontsize=7, bbox_to_anchor=(0.5, 1.07))
    axes[0].text(0.02, 0.97, "n = 100 queries per set", transform=axes[0].transAxes, fontsize=5.5,
                 color=MUTED, va="top")
    fig.tight_layout()
    save(fig, "fig_breaker")


def fig_components() -> None:
    """Median value of each estimator component per basis on B2 (dot plot; log axis, zeros drawn hollow at the floor)."""
    path = RES / "b2_e3" / "e3_b2.npz"
    comps = (("eta_eq", r"$\eta_{eq}$"), ("eta_sign", r"$\eta_{sign}$"), ("eta_gap", r"$\eta_{gap}$"),
             ("eta_comp", r"$\eta_{comp}$"))
    floor = 1e-14
    fig, axes = plt.subplots(1, 4, figsize=(6.3, 1.8), sharey=True)
    for ax, (key, title) in zip(axes, comps):
        for row, kind in enumerate(KINDS):
            v = float(np.median(pooled(path, kind, key)))
            zero = v <= floor
            ax.plot(max(v, floor), row, "o", ms=4.5, color=BLUE, mfc="white" if zero else BLUE, mew=1.0)
            if zero:
                ax.text(floor * 3, row, "0", fontsize=6.5, va="center", color=MUTED)
        ax.set(xscale="log", xlim=(floor / 3, 1e1), title=title)
        ax.set_yticks(range(len(KINDS)), [LABEL[k] for k in KINDS])
        ax.invert_yaxis()
        ax.tick_params(labelsize=6.5)
    fig.supxlabel("median over test queries and dual dimensions (B2)", fontsize=7, y=-0.04)
    fig.tight_layout()
    save(fig, "fig_components")


if __name__ == "__main__":
    for make in (fig_silent_failure, fig_bound, fig_breakeven, fig_breaker, fig_components):
        make()
