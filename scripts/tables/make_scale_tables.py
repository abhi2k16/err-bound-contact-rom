"""Tables for the out-of-range and mesh-size studies -> ``manuscript_main/tables/extrapolation.tex`` and ``manuscript_main/tables/meshsize.tex``."""

import json
from pathlib import Path

KINDS = ("svd", "nmf", "cpg", "gis", "rid")
LAB = {"svd": "SVD", "nmf": "NMF", "cpg": "CPG", "gis": "gIS", "rid": "RID"}
SETS = {"b1": (("test", "in range"), ("load", r"load $\times$1.3--2"), ("geometry", "steeper obstacle")),
        "b2": (("test", "in range"), ("load", r"load $\times$1.3--2"), ("centre shift", "centre outside")),
        "b4": (("test", "in range"), ("load", r"load $\times$1.3--2"), ("centre shift", "centre outside"))}


def extrapolation_table(benches: tuple, out: str) -> None:
    lines = []
    for n, bench in enumerate(benches):
        rows = json.load(open(f"results/scale_study/extrapolation_{bench}.json"))["rows"]
        for i, (key, name) in enumerate(SETS[bench]):
            cells = []
            for kind in KINDS:
                r = rows[f"{kind}|{key}"]
                cells.append(f"{100 * r['share_above_5pct']:.0f} ({r['missed']}/{r['n_silent']})")
            clam = min(rows[f"{k}|{key}"]["bound_Clam_valid"] for k in KINDS)
            lines.append(f"{bench.upper() if i == 0 else ''} & {name} & " + " & ".join(cells) + f" & {100 * clam:.1f} \\\\")
        if n < len(benches) - 1:
            lines.append(r"\midrule")
    head = (r"\begin{tabular}{@{}llcccccc@{}}" "\n" r"\toprule" "\n"
            r"Bench. & Query set & SVD & NMF & CPG & gIS & RID & $C_\lambda$ bound valid (\%) \\" "\n" r"\midrule" "\n")
    Path(out).write_text(head + "\n".join(lines) + "\n" r"\bottomrule" "\n" r"\end{tabular}" "\n", encoding="utf-8")


def extrapolation() -> None:
    """Main-text table (B1, B2) and the supplementary table of the independent benchmark B4."""
    extrapolation_table(("b1", "b2"), "manuscript_main/tables/extrapolation.tex")
    extrapolation_table(("b4",), "manuscript_main/tables/extrapolation_b4.tex")


def meshsize() -> None:
    lines = []
    for bench in ("b1", "b2"):
        rows = json.load(open(f"results/scale_study/mesh_{bench}.json"))["rows"]
        sizes = sorted({r["N"] for r in rows.values()})
        for i, n in enumerate(sizes):
            g = lambda k, m: rows[f"{k}|{n}"][m]
            lines.append(f"{bench.upper() if i == 0 else ''} & {n} & {rows[f'nmf|{n}']['n_c']} & "
                         f"{100 * g('svd', 'share_above_5pct'):.0f} & {100 * g('nmf', 'share_above_5pct'):.0f} ({g('nmf', 'err_med'):.2g}) & "
                         f"{100 * g('gis', 'share_above_5pct'):.0f} & {g('nmf', 'S_aff'):.2f} & {g('gis', 'S_aff'):.2f} \\\\")
        if bench == "b1":
            lines.append(r"\midrule")
    head = (r"\begin{tabular}{@{}lrrcccrr@{}}" "\n" r"\toprule" "\n"
            r"Bench. & $N$ & $n_c$ & SVD $>5\%$ & NMF $>5\%$ (med.) & gIS $>5\%$ & $S$ NMF & $S$ gIS \\" "\n" r"\midrule" "\n")
    Path("manuscript_main/tables/meshsize.tex").write_text(head + "\n".join(lines) + "\n" r"\bottomrule" "\n" r"\end{tabular}" "\n", encoding="utf-8")


if __name__ == "__main__":
    extrapolation()
    meshsize()
    print(Path("manuscript_main/tables/extrapolation.tex").read_text())
