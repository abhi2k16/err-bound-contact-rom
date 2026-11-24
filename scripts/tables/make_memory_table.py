"""Storage comparison, full-order sparse factorization vs reduced dense model -> ``results/memory_study/memory_{b2,b3}.json``
and ``manuscript_main/tables/memory_table.tex``.

Run ``python -m scripts.tables.make_memory_table``. For every mesh size of the break-even sweeps (B2 membrane, B3 two-body) the
stiffness matrix is factored with the same SuperLU call as in the timing study (``scipy.sparse.linalg.splu``,
default COLAMD ordering) and the stored entries of L and U are counted. Reduced sizes ``r_u`` come from the timing
study (primal POD dimension plus repair vectors); ``r_lambda = cfg.e4_r_dual`` is the dual dimension used there.
Counts are numbers of stored double-precision values (8 bytes each); index arrays are not counted.

full-order : nnz(L) + nnz(U) + n_c^2 (the dense Schur block A = B K^-1 B^T that the contact solver uses)
reduced    : r_u^2 + r_l^2 + n_c r_l (reduced K, reduced dual block, dual basis) + n_c r_u (B V)   [online screening and solve]
reduced+lift: reduced + N r_u (primal basis, needed to reconstruct the full field)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from contact_rom.break_even_experiment import sweep_problem, sweep_problem_b3, sweep_problem_b3_code
from contact_rom.config import Config

BYTES = 8


def measure(problem, r_u: int, r_l: int) -> dict:
    mu = problem.sample_mu(np.random.default_rng(0), 1)[0]
    K, B = sp.csc_matrix(problem.K_mat(mu)), problem.B_mat(mu)
    lu = spla.splu(K)
    N, n_c = problem.n_dof, problem.n_contact
    full = lu.L.nnz + lu.U.nnz + n_c ** 2
    red = r_u ** 2 + r_l ** 2 + n_c * r_l + n_c * r_u
    return {"N": N, "n_c": n_c, "nnz_K": int(K.nnz), "nnz_LU": int(lu.L.nnz + lu.U.nnz), "full": int(full),
            "reduced": int(red), "reduced_lift": int(red + N * r_u), "r_u": r_u, "r_l": r_l}


def run(bench: str, cfg: Config) -> list:
    rows_t = {r["N"]: r["r_u"] for r in json.load(open(f"results/timing_study/timing_{bench}.json"))["rows"]}
    if bench == "b2":
        jobs = [sweep_problem(cfg, g) for g in cfg.e4_grids]
    else:
        jobs = [sweep_problem_b3(cfg, g) for g in cfg.e4_grids_b3] + [sweep_problem_b3_code(cfg, c) for c in cfg.e4_codes_b3_small]
    out = []
    for p in jobs:
        if p.n_dof in rows_t:           # sizes with contact in the training snapshots, as in the timing study
            out.append(measure(p, rows_t[p.n_dof], cfg.e4_r_dual))
    out.sort(key=lambda r: r["N"])
    return out


def mb(n: float) -> str:
    v = n * BYTES / 2 ** 20
    return f"{v:.3f}" if v < 0.1 else f"{v:.2f}" if v < 10 else f"{v:.0f}"


def main() -> None:
    cfg = Config()
    Path("results/memory_study").mkdir(parents=True, exist_ok=True)
    lines, summary = [], {}
    for bench, name in (("b2", "B2"), ("b3", "B3")):
        rows = run(bench, cfg)
        json.dump(rows, open(f"results/memory_study/memory_{bench}.json", "w"), indent=1)
        N = np.array([r["N"] for r in rows], float)
        lu = np.array([r["nnz_LU"] for r in rows], float)
        alpha = np.polyfit(np.log(N[N >= 200]), np.log(lu[N >= 200]), 1)[0]
        summary[bench] = {"alpha_LU": float(alpha)}
        print(bench, "fitted exponent of nnz(L+U) in N (N>=200):", round(alpha, 2))
        for i, r in enumerate(rows):
            if r["N"] < 100:
                continue
            lines.append(f"{name if i == [x for x in range(len(rows)) if rows[x]['N'] >= 100][0] else ''} & {r['N']} & {r['n_c']} & "
                         f"{r['nnz_LU'] / r['nnz_K']:.1f} & {mb(r['full'])} & {mb(r['reduced'])} & {mb(r['reduced_lift'])} & "
                         f"{r['full'] / r['reduced']:.1f} & {r['full'] / r['reduced_lift']:.1f} \\\\")
        if bench == "b2":
            lines.append(r"\midrule")
    json.dump(summary, open("results/memory_study/summary.json", "w"))
    rows_tex = [
        r"\begin{tabular}{@{}lrrrrrrrr@{}}", r"\toprule",
        r" & & & fill & \multicolumn{3}{c}{storage (MiB)} & \multicolumn{2}{c}{ratio} \\",
        r"\cmidrule(lr){5-7}\cmidrule(lr){8-9}",
        r"Bench. & $N$ & $n_c$ & $\mathrm{nnz}(L{+}U)/\mathrm{nnz}(K)$ & full-order & reduced & reduced + lift & online & with lift \\",
        r"\midrule", *lines, r"\bottomrule", r"\end{tabular}"]
    body = "\n".join(rows_tex) + "\n"
    Path("manuscript_main/tables").mkdir(parents=True, exist_ok=True)
    Path("manuscript_main/tables/memory_table.tex").write_text(body, encoding="utf-8")
    print(body)


if __name__ == "__main__":
    main()
