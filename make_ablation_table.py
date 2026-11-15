"""Ablation table from ``results/ablation/ablation_{b1,b2}.json`` -> ``tables/ablation.tex``."""

import json
from pathlib import Path

LABEL = {"svd": "SVD", "nmf": "NMF", "cpg": "CPG", "gis": "gIS"}


def fmt_cond(x: float) -> str:
    return r"$\infty$" if x == float("inf") or x > 1e30 else f"{x:.0e}".replace("e+0", "e").replace("e+", "e")


def main() -> None:
    lines = []
    for bench in ("b1", "b2"):
        rows = json.load(open(f"results/ablation/ablation_{bench}.json"))["rows"]
        first = True

        def add(name, basis, r):
            nonlocal first
            fb = f"{r['fallback_share']:.2f}" if "fallback_share" in r else "--"
            miss = f"{r['accepted_above_5pct']}" if "accepted_above_5pct" in r else "--"
            lines.append(f"{bench.upper() if first else ''} & {name} & {basis} & {r['err_med']:.2g} & {100 * r['share_above_5pct']:.0f} & "
                         f"{fmt_cond(r['cond_med'])} & {r['r_u']:.0f} & {r['flops'] / 1e3:.0f} & {r['speedup']:.2f} & {fb} & {miss} \\\\")
            first = False

        add("Primal only", "--", rows["primal|-"])
        for kind in ("svd", "nmf", "cpg", "gis"):
            add("Primal--dual, no repair", LABEL[kind], rows[f"pd_norepair|{kind}"])
            add("Primal--dual, repair", LABEL[kind], rows[f"pd_repair|{kind}"])
            add("Repair + screen + fallback", LABEL[kind], rows[f"screened|{kind}"])
        if bench == "b1":
            lines.append(r"\midrule")
    head = (r"\begin{tabular}{@{}lllrrrrrrrr@{}}" "\n" r"\toprule" "\n"
            r"Bench. & Variant & Basis & Median err. & $>5\%$ (\%) & cond. & $r_u$ & kFLOPs & Speedup & Fallback & Missed \\" "\n" r"\midrule" "\n")
    Path("tables/ablation.tex").write_text(head + "\n".join(lines) + "\n" r"\bottomrule" "\n" r"\end{tabular}" "\n", encoding="utf-8")
    print(Path("tables/ablation.tex").read_text())


if __name__ == "__main__":
    main()
