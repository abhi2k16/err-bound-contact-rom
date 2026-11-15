"""Threshold sensitivity, false-negative / acceptance trade-off, rule comparison and net speedup with fallback.

Reads ``results/calibration/audits_*.pkl`` and ``analysis_*.pkl`` (``python -m contact_rom calibration``), and the
membrane timings ``results/b2_e4/e4_b2.npz``. Writes ``figures/fig_threshold_sensitivity.pdf``,
``figures/fig_tradeoff.pdf``, ``tables/threshold_rules.tex`` and ``tables/net_speedup.tex``.
"""

from __future__ import annotations

import pickle
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import make_figures as mf
from contact_rom import calibration_study as cs
from contact_rom.config import Config

CFG = Config()
SEEDS = tuple(range(42, 50))
KINDS = ("nmf", "cpg", "gis", "rid")
IND = "eta_adm"
REL = np.geomspace(1e-3, 1e2, 120)       # threshold / tau_rel grid
STYLE = {"nmf": (mf.BLUE, "o"), "cpg": (mf.ORANGE, "s"), "gis": (mf.AQUA, "^"), "rid": (mf.INK, "D")}


def load(bench: str) -> dict:
    return {s: pickle.load(open(Path("results/calibration") / f"audits_{bench}_seed{s}.pkl", "rb")) for s in SEEDS}


def sweep(data: dict, kind: str) -> dict:
    """Pooled FNR, FPR and accepted share over seeds and dual dimensions for every relative threshold."""
    fn = np.zeros(len(REL)); pos = 0
    fp = np.zeros(len(REL)); neg = 0
    acc = np.zeros(len(REL)); tot = 0
    for payload in data.values():
        for (k, _), test in payload["test"].items():
            if k != kind:
                continue
            ind = test[IND]
            silent = cs.silent_failure_mask(test, CFG)
            good = test["err_energy"] <= CFG.err_silent
            tau = REL * payload["tau_rel"]
            flag = ind[None, :] > tau[:, None]
            fn += (silent[None, :] & ~flag).sum(axis=1); pos += silent.sum()
            fp += (good[None, :] & flag).sum(axis=1); neg += good.sum()
            acc += (~flag).sum(axis=1); tot += len(ind)
    return {"fnr": fn / max(pos, 1), "fpr": fp / max(neg, 1) if neg else np.full(len(REL), np.nan), "acc": acc / tot}


def fig_sensitivity() -> None:
    fig, axes = plt.subplots(2, len(KINDS), figsize=(7.0, 3.4), sharex=True, sharey=True)
    for row, bench in enumerate(("b1", "b2")):
        data = load(bench)
        for col, kind in enumerate(KINDS):
            ax = axes[row, col]
            sw = sweep(data, kind)
            ax.plot(REL, sw["fnr"], color=mf.BLUE, lw=1.2, label="false-negative rate")
            ax.plot(REL, sw["fpr"], color=mf.ORANGE, lw=1.2, ls="--", label="false-positive rate")
            ax.plot(REL, sw["acc"], color=mf.AQUA, lw=1.2, ls=":", label="accepted by ROM")
            ax.axvline(1.0, color=mf.MUTED, lw=0.6)
            ax.set_xscale("log")
            ax.set_ylim(-0.02, 1.02)
            if row == 0:
                ax.set_title(mf.LABEL[kind], fontsize=8)
            if col == 0:
                ax.set_ylabel(f"{bench.upper()}: share")
            if row == 1:
                ax.set_xlabel(r"$\tau_{\mathrm{adm}}/\tau_{\mathrm{rel}}$")
    axes[0, 0].legend(fontsize=5.8, loc="center left")
    fig.tight_layout()
    mf.save(fig, "fig_threshold_sensitivity")


def fig_tradeoff() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6), sharey=True)
    for ax, bench in zip(axes, ("b1", "b2")):
        data = load(bench)
        for kind in KINDS:
            colour, marker = STYLE[kind]
            sw = sweep(data, kind)
            ax.plot(sw["acc"], np.maximum(sw["fnr"], 1e-4), color=colour, lw=1.2, label=mf.LABEL[kind])
        ax.set_yscale("log")
        ax.set_xlabel("share of queries accepted by the ROM")
        ax.set_title(bench.upper(), fontsize=8)
    axes[0].set_ylabel("false-negative rate")
    axes[0].legend(fontsize=6.5, loc="lower right")
    fig.tight_layout()
    mf.save(fig, "fig_tradeoff")


def rule_table(n_val: int = 40) -> None:
    rules = ("min09", "q10", "q25", "costopt", "hybrid")
    lines = []
    for bench in ("b1", "b2"):
        st = pickle.load(open(Path("results/calibration") / f"analysis_{bench}.pkl", "rb"))
        for i, kind in enumerate(KINDS):
            cells = []
            for rule in rules:
                s = st[(kind, IND, rule, n_val)]
                cells.append(f"{s['recall']:.2f} / {s['fpr']:.2f}" + (f" ({s['uncal_share']:.2f})" if rule != "hybrid" and s["uncal_share"] > 0 else ""))
            lines.append(f"{bench.upper() if i == 0 else ''} & {mf.LABEL[kind]} & " + " & ".join(cells) + r" \\")
        if bench == "b1":
            lines.append(r"\midrule")
    head = (r"\begin{tabular}{@{}llccccc@{}}" "\n" r"\toprule" "\n"
            r"Bench. & Basis & $0.9\min$ & 10th pct. & 25th pct. & cost-optimal & hybrid \\" "\n" r"\midrule" "\n")
    Path("tables/threshold_rules.tex").write_text(head + "\n".join(lines) + "\n" r"\bottomrule" "\n" r"\end{tabular}" "\n",
                                                  encoding="utf-8")


def net_speedup_table(n_val: int = 40, draws: int = 20) -> None:
    """Membrane: net speedup ``1 / (1/S + f)`` when the hybrid-rule fallback share ``f`` of queries is re-solved."""
    e4 = np.load("results/b2_e4/e4_b2.npz")
    data = load("b2")
    rng = np.random.default_rng(0)
    lines = []
    for kind in KINDS:
        flagged, total = 0, 0
        for payload in data.values():
            for (k, r), val_full in payload["val"].items():
                if k != kind:
                    continue
                test = payload["test"][(k, r)]
                n_pool = len(val_full["err_energy"])
                for _ in range(draws):
                    idx = rng.choice(n_pool, n_val, replace=False)
                    val = {a: v[idx] for a, v in val_full.items() if np.ndim(v) == 1 and len(v) == n_pool}
                    tau = cs.threshold("hybrid", val, IND, CFG, payload["tau_rel"])
                    flagged += int((test[IND] > tau).sum()); total += len(test[IND])
        f = flagged / total
        s_aff, s_var = float(e4[f"{kind}_16_S_meas_aff"]), float(e4[f"{kind}_16_S_meas_var"])
        lines.append(f"{mf.LABEL[kind]} & {f:.2f} & {s_aff:.2f} & {1 / (1 / s_aff + f):.2f} & {s_var:.1f} & {1 / (1 / s_var + f):.1f} \\\\")
    head = (r"\begin{tabular}{@{}lccccc@{}}" "\n" r"\toprule" "\n"
            r"Basis & fallback share $f$ & $S$ (const. $K$) & net & $S$ (varying $K$) & net \\" "\n" r"\midrule" "\n")
    Path("tables/net_speedup.tex").write_text(head + "\n".join(lines) + "\n" r"\bottomrule" "\n" r"\end{tabular}" "\n", encoding="utf-8")
    print(Path("tables/net_speedup.tex").read_text())


if __name__ == "__main__":
    fig_sensitivity()
    fig_tradeoff()
    rule_table()
    net_speedup_table()
