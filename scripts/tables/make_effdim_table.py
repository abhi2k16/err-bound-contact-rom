"""Nominal vs effective dimension, rank repair, offline and online cost at nominal r_lambda = 16 -> ``tables/basis_effective.tex``."""

import json
from pathlib import Path

LABEL = {"svd": "SVD + repair", "nmf": "NMF", "cpg": "CPG", "gis": "gIS", "rid": "RID"}


def main() -> None:
    lines = []
    for bench in ("b1", "b2"):
        rows = json.load(open(f"results/basis_comparison/basis_comparison_{bench}.json"))["rows"]
        for i, kind in enumerate(LABEL):
            e = [x for x in rows[kind] if x["nominal"] == 16][0]
            lines.append(f"{bench.upper() if i == 0 else ''} & {LABEL[kind]} & 16 & {e['r_eff']} & {e['r_rep']} & {e['r_u']} & "
                         f"{e['t_offline'] * 1e3:.0f} & {e['flops'] / 1e3:.0f} & {e['t_online'] * 1e6:.0f} \\\\")
        if bench == "b1":
            lines.append(r"\midrule")
    head = (r"\begin{tabular}{@{}llcccccrr@{}}" "\n" r"\toprule" "\n"
            r"Bench. & Method & Nominal $r_\lambda$ & Eff. $r_\lambda$ & Repair $r_{\mathrm{rep}}$ & $r_u$ & "
            r"Offline (ms) & Online (kFLOPs) & Online ($\mu$s) \\" "\n" r"\midrule" "\n")
    Path("tables/basis_effective.tex").write_text(head + "\n".join(lines) + "\n" r"\bottomrule" "\n" r"\end{tabular}" "\n", encoding="utf-8")
    print(Path("tables/basis_effective.tex").read_text())


if __name__ == "__main__":
    main()
