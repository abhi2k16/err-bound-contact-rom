"""Table of the repeated-timing study (``results/timing_study/analysis_{b2,b3}.json``) -> ``tables/timing_table.tex``.

Run ``python -m contact_rom timingstudy --bench b2`` and ``--bench b3`` (one at a time, nothing else running), then
the analysis (``analyse`` in ``contact_rom/timing_study.py``) and this script.
"""

import json
from pathlib import Path

import numpy as np

from contact_rom.acceptance import criterion4

E4 = {"b2": "results/b2_e4/e4_sweep.npz", "b3": "results/b3_e4/e4_sweep_b3.npz"}
NAME = {"b2": "B2", "b3": "B3"}
REG = {"aff": "constant $K$", "var": "varying $K$"}


def fmt(x: float, spec: str = ".0f") -> str:
    """Format a number, or ``n/d`` when the break-even size is not defined (no crossing in the swept range)."""
    return "n/d" if x != x else format(x, spec)


def main() -> None:
    lines = []
    for bench in ("b2", "b3"):
        a = json.load(open(f"results/timing_study/analysis_{bench}.json"))
        e4 = criterion4(Path(E4[bench]))["cross"]
        ps = a["per_size"]
        cv = {k: (np.min([p[k] for p in ps]), np.max([p[k] for p in ps])) for k in ("cv_rom", "cv_aff", "cv_full")}
        last = ps[-1]
        for i, reg in enumerate(("aff", "var")):
            c = a["crossings"][reg]
            lo, hi = c["N_meas_ci"]
            share = c["n_defined"] / c["n_boot"]
            lines.append(
                f"{NAME[bench] if i == 0 else ''} & {REG[reg]} & {c['N_pred']:.0f} & {fmt(e4[reg]['N_meas'])} & "
                f"{fmt(c['N_meas'])} [{fmt(lo)}, {fmt(hi)}]{'$^{*}$' if share < 0.9 else ''} & "
                f"{fmt(c['N_pred'] / c['N_meas'], '.2f')} & {last['S_pred_' + reg]:.0f} & "
                f"{last['S_meas_' + reg]:.0f} [{last['S_ci_' + reg][0]:.0f}, {last['S_ci_' + reg][1]:.0f}] \\\\")
        lines.append(rf"\multicolumn{{8}}{{@{{}}l}}{{\scriptsize {NAME[bench]}: CV across repeats "
                     f"{cv['cv_rom'][0]:.2f}--{cv['cv_rom'][1]:.2f} (reduced), {cv['cv_aff'][0]:.2f}--{cv['cv_aff'][1]:.2f} "
                     f"(full-order, constant $K$), {cv['cv_full'][0]:.2f}--{cv['cv_full'][1]:.2f} (full-order, varying $K$)}} \\\\")
        if bench == "b2":
            lines.append(r"\midrule")
    head = (r"\begin{tabular}{@{}llccccrr@{}}" "\n" r"\toprule" "\n"
            r" & Regime & $N^*_{\mathrm{pred}}$ & $N^*_{\mathrm{E4}}$ & $N^*_{\mathrm{meas}}$ [95\%] & "
            r"pred/meas & $S_{\mathrm{pred}}$ & $S_{\mathrm{meas}}$ [95\%] \\" "\n" r"\midrule" "\n")
    lines.append(r"\multicolumn{8}{@{}l}{\scriptsize $^{*}$ Break-even size defined in fewer than 90\% of the bootstrap resamples. n/d: no crossing in the swept range.} \\")
    Path("tables/timing_table.tex").write_text(head + "\n".join(lines) + "\n" r"\bottomrule" "\n" r"\end{tabular}" "\n",
                                               encoding="utf-8")
    print(Path("tables/timing_table.tex").read_text())


if __name__ == "__main__":
    main()
