"""Estimator-performance figure: admissibility indicator eta_adm vs true energy-norm error.

Run ``python make_estimator_figure.py`` (needs results/b1_e3, results/b2_e3). Writes ``figures/fig_estimator.pdf``.
Silent failure = reduced KKT residual < 1e-12 and relative energy error > 5% (same definition as elsewhere).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import make_figures as mf

FLOOR = 1e-9   # indicator values of exactly zero are drawn at the axis floor


def fig_estimator() -> None:
    fig, axes = plt.subplots(2, 5, figsize=(6.6, 3.5), sharey=True)
    for row, (bench, path) in enumerate((("B1", mf.RES / "b1_e3" / "e3_b1.npz"), ("B2", mf.RES / "b2_e3" / "e3_b2.npz"))):
        for col, kind in enumerate(mf.KINDS):
            ax = axes[row, col]
            err = np.maximum(mf.pooled(path, kind, "err_energy"), 1e-6)
            eta = np.maximum(mf.pooled(path, kind, "eta_adm"), FLOOR)
            res = mf.pooled(path, kind, "res_red")
            sf = (err > mf.TOL_ERR) & (res < mf.TOL_KKT)
            ax.scatter(err[~sf], eta[~sf], s=4, color=mf.BLUE, alpha=0.6, linewidths=0, label="accurate")
            ax.scatter(err[sf], eta[sf], s=5, marker="^", color=mf.ORANGE, alpha=0.7, linewidths=0, label="silent failure")
            ax.axvline(mf.TOL_ERR, color=mf.MUTED, lw=0.6, ls="--")
            ax.set(xscale="log", yscale="log", xlim=(1e-6, 5), ylim=(FLOOR / 3, 60))
            ax.tick_params(labelsize=6.5)
            if row == 0:
                ax.set_title(mf.LABEL[kind], fontsize=8)
            if col == 0:
                ax.set_ylabel(rf"{bench}: $\eta_{{\mathrm{{adm}}}}$")
            if row == 1:
                ax.set_xlabel("true error", fontsize=7)
    axes[0, 0].legend(loc="lower right", fontsize=6, markerscale=1.5, handletextpad=0.2)
    fig.tight_layout()
    mf.save(fig, "fig_estimator")


if __name__ == "__main__":
    fig_estimator()
