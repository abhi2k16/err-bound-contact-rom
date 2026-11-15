import io
p = r"C:\Users\abhij\Desktop\Physical_ML\proposed-research\reduced-contact-reliability\manuscript_master_v2.tex"
CR, LF = chr(13), chr(10)
raw = io.open(p, encoding="utf-8", newline="").read()
crlf = (CR + LF) in raw
s = raw.replace(CR + LF, LF)


def rep(old, new):
    global s
    assert s.count(old) == 1, (s.count(old), old[:80])
    s = s.replace(old, new, 1)


anchor = r"\paragraph{Comparison at matched effective dimension, cost and accuracy.}"
block = r"""\paragraph{Ablation of the framework components.}
To isolate what each component contributes, we compared four configurations on B1 and B2 (constant stiffness; five seeds with 40 test and 40 validation queries each, so 200 queries per row; primal dimension 20, dual dimension 16). (i)~A \emph{primal-only} ROM reduces the displacements and keeps all $n_c$ multipliers. (ii)~A \emph{primal--dual} ROM with a dual basis (signed SVD, NMF, CPG or gIS) and no primal repair. (iii)~The same with supremizer rank repair. (iv)~Configuration (iii) with the interface screen of Eq.~\eqref{eq:eta-adm}, the hybrid threshold calibrated on the validation queries, and a full-order re-solve of every flagged query; its error and speedup are those of the answers actually served. Table~\ref{tab:ablation} gives the results.

\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{3.5pt}
\caption{Ablation on B1 and B2 (200 queries per row). Median err.: median relative energy error of the served answer; $>5\%$: share of served answers above $5\%$ error; cond.: median condition number of the reduced Schur complement (infinite when singular); speedup: mean full-order time over mean time of the configuration (net of the fallback in the last variant, whose served errors are zero for the flagged queries); Fallback: share of queries re-solved at full order; Missed: accepted queries with error above $5\%$. Speedups come from one timing per query, and models that are identical (for example NMF with and without repair, where repair adds no vector) differ by up to about $10\%$.}
\label{tab:ablation}
\input{tables/ablation.tex}
\end{table}

"""
rep(anchor, block + anchor)

item_anchor = "\\item\n  \\textbf{Matched comparison (Tables"
item = r"""\item
  \textbf{Ablation (Table~\ref{tab:ablation}).} \emph{Primal-only ROM.} Reducing only the displacements keeps primal feasibility and non-negativity by construction but is neither accurate nor cheap here: the reduced Schur complement is singular ($n_c > r_u$), the median error is $0.0023$ (B1) and $0.031$ (B2) but $44\%$ (B1) and $11\%$ (B2) of the queries exceed $5\%$, with errors up to $692$ on B1 (we attribute this to the singular LCP but did not isolate the cause), and the dense $n_c$-dimensional problem makes the ROM $25$ to $100$ times slower than the full-order solver (speedup $0.04$ and $0.01$). \emph{Dual reduction} is therefore what makes the reduced model cheap, and the dual basis is what decides whether it is admissible: the signed SVD basis fails on all queries ($100\%$ above $5\%$, median error $0.93$--$0.95$) while NMF and CPG give median errors of $0.002$--$0.004$ (B1) and $0.04$--$0.05$ (B2). \emph{Rank repair} matters for solvability, not for accuracy: it adds no vector for SVD and NMF at $r_\lambda = 16$ and leaves CPG unchanged, and only for gIS does it turn a singular Schur complement into a regular one (condition number $\infty \to 1.4\times10^{3}$ on B1, $\infty \to 1.2\times10^{5}$ on B2). On B1 this also lowers the median error from $0.18$ to $0.061$ and the share above $5\%$ from $68\%$ to $52\%$ and removes the one failed LCP solve, whereas on B2 the error is unchanged ($0.19$). \emph{Screen and fallback} converts silent failures into fallbacks at a price: the signed SVD ROM is flagged on every query (fallback $1.00$), which is safe but makes it useless (net speedup $0.36$ and $0.66$); on B2 NMF the share above $5\%$ falls from $26\%$ to $3\%$ at a fallback share of $0.39$ with six of the 200 accepted answers still above $5\%$ (largest $6.1\%$), CPG falls to zero at a fallback share of $0.83$, and gIS leaves seven accepted answers above $5\%$ (largest $7.2\%$) at a fallback share of $0.90$. On B1 no screened configuration of NMF or CPG needs a fallback because their errors are below $5\%$, but the unscreened speedup is already below one there, because the problem is too small. The net speedup of the screened ROM is below one for every basis on both benchmarks ($0.24$--$0.89$ in constant-stiffness regime at $N = 324$), consistent with Table~\ref{tab:netspeed}.
"""
rep(item_anchor, item + item_anchor)
io.open(p, "w", encoding="utf-8", newline="").write(s.replace(LF, CR + LF) if crlf else s)
print("ok")
