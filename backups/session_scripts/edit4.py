import io
p = r"C:\Users\abhij\Desktop\Physical_ML\proposed-research\reduced-contact-reliability\manuscript_master_v2.tex"
CR, LF = chr(13), chr(10)
raw = io.open(p, encoding="utf-8", newline="").read()
crlf = (CR + LF) in raw
s = raw.replace(CR + LF, LF)

head = r"\begin{itemize}" + LF + r"\tightlist" + LF + r"\item" + LF
anchor = head + r"  \textbf{Calibration study (Table~\ref{tab:calibration}, Figure~\ref{fig:calibration}).}"
assert s.count(anchor) == 1

block = r"""\paragraph{Comparison at matched effective dimension, cost and accuracy.}
Table~\ref{tab:e3} compares the bases at the same nominal $r_\lambda$, but the effective sizes differ: RID keeps all contact rows of the selected elements (23 rows on B1 and 98 on B2 at nominal 16), and gIS and RID need more primal rank-repair modes. We therefore audited all five bases on the same snapshots over $r_\lambda \in \{4,\dots,30\}$ (B1) and $\{4,\dots,60\}$ (B2), 100 test queries, with rank repair, and compare them on common axes (Tables~\ref{tab:matcheddim} and \ref{tab:matchedacc}, Figure~\ref{fig:basiscost}). The online cost is the FLOP count of Eq.~\eqref{eq:crom} for constant $K$ (reduced dense solves, reduced active-set iterations and the interface check), evaluated with the measured median iteration count and active-set size, and the measured median wall time of the reduced solve. A dual dimension above the number of training snapshots (60) is not available to the snapshot-span bases; NMF and CPG were still run up to 60.

\begin{table}[htbp]
\centering
\small
\caption{Median relative energy error at matched effective dual dimension $r_\lambda^{\mathrm{eff}}$, with the effective dual and enriched primal dimensions in parentheses. The entry nearest to the target within $25\%$ is used; ``n/a'' means the basis cannot be built at that size (RID on B2 keeps at least 61 rows). Signed SVD is the only basis whose error stays near one.}
\label{tab:matcheddim}
\input{tables/basis_matched_dim.tex}
\end{table}

\begin{table}[htbp]
\centering
\small
\caption{Cheapest entry whose median error is at most the target: effective dual / enriched primal dimension, online FLOPs per query (thousands, Eq.~\eqref{eq:crom}), measured reduced-solve time in $\mu$s. ``Not reached'' means no entry up to $r_\lambda = 30$ (B1) or $60$ (B2) attains the target.}
\label{tab:matchedacc}
\input{tables/basis_matched_acc.tex}
\end{table}

\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{fig_basis_cost.pdf}
\caption{Median relative energy error against online FLOPs per query for the five bases, each point one nominal dual dimension (with rank repair). The dotted line is the $5\%$ silent-failure level. RID moves right quickly because its effective dual dimension and rank-repair modes grow with the nominal value.}
\label{fig:basiscost}
\end{figure}

"""

item = r"""\item
  \textbf{Matched comparison (Tables~\ref{tab:matcheddim}, \ref{tab:matchedacc}, Figure~\ref{fig:basiscost}).} At the same effective dual dimension the ranking in Table~\ref{tab:e3} changes in size but not in kind for NMF and CPG. At $r_\lambda^{\mathrm{eff}} = 16$ the median errors on B1 are $0.0021$ (NMF), $0.003$ (CPG), $0.051$ (gIS) and $0.26$ (RID, nominal 8), and on B2 $0.038$ (NMF), $0.054$ (CPG) and $0.24$ (gIS), with RID unavailable. The comparison reverses once the dimension is large enough: on B1 at $r_\lambda^{\mathrm{eff}} = 24$ gIS reaches $2.9\times10^{-6}$ and RID $3.8\times10^{-4}$ against $9.3\times10^{-4}$ for NMF, and on B2 at $r_\lambda^{\mathrm{eff}} \approx 60$ the four bases are within a factor of 2.5 of each other ($0.013$ to $0.033$, gIS best). Online FLOPs at equal effective dimensions are nearly equal across bases, so the cost differences in Figure~\ref{fig:basiscost} come from the effective dimensions that each basis needs: to reach a $5\%$ median error, NMF needs 12 dual modes (28k FLOPs) on B2 against 24 for CPG, 48 for gIS (740k) and 61 contact rows for RID (1183k). A $1\%$ median error on B2 is reached only by RID (82 rows and 82 primal modes, 2732k FLOPs), because NMF and CPG reach only $1.7\%$ and $3.3\%$ at the largest dimension tried and gIS reaches $1.3\%$ only there. The nominal-dimension comparison thus understated the price of RID and gIS on B2 and hid their advantage at high accuracy; the conclusion that NMF is the cheapest way to a moderate accuracy, and that RID or gIS are needed for high accuracy, holds on both benchmarks. These are results for one seed and one snapshot set per benchmark.
"""

s = s.replace(anchor, block + head + item + r"""\item
  \textbf{Calibration study (Table~\ref{tab:calibration}, Figure~\ref{fig:calibration}).}""", 1)
io.open(p, "w", encoding="utf-8", newline="").write(s.replace(LF, CR + LF) if crlf else s)
print("ok")
