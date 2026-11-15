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
block = r"""\paragraph{Threshold sensitivity, alternative rules and the cost of fallback.}
Figure~\ref{fig:thrsens} shows how the false-negative rate, the false-positive rate and the share of queries accepted by the ROM change with the threshold, written as $\tau_{\mathrm{adm}}/\tau_{\mathrm{rel}}$ with $\tau_{\mathrm{rel}}$ the relative threshold of the calibration study; counts are pooled over seeds and dual dimensions on the test sets. Figure~\ref{fig:tradeoff} plots the false-negative rate against the share of queries accepted. In both, a rate of silent failures that is not small makes the trade-off steep: on B1 NMF and CPG accept about $90\%$ of the queries before the first miss, whereas gIS and RID accept only about $40\%$ because their silent failures are numerous; on B2 a false-negative rate of $1\%$ requires sending about $65\%$ of the NMF queries, and more than $90\%$ of the gIS queries, to the full-order solver.

\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{fig_threshold_sensitivity.pdf}
\caption{False-negative rate, false-positive rate and share of queries accepted by the ROM as functions of the threshold of Eq.~\eqref{eq:eta-adm}, relative to $\tau_{\mathrm{rel}}$ (vertical line). Eight seeds, five dual dimensions pooled, test sets.}
\label{fig:thrsens}
\end{figure}

\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{fig_tradeoff.pdf}
\caption{False-negative rate against the share of queries accepted by the ROM, obtained by sweeping the threshold of Eq.~\eqref{eq:eta-adm} (eight seeds, five dual dimensions pooled). The floor of the vertical axis is $10^{-4}$.}
\label{fig:tradeoff}
\end{figure}

Table~\ref{tab:rules} compares five calibration rules on the validation queries (40 per seed). Besides the $0.9\,\min$ and hybrid rules, two percentile rules set $\tau_{\mathrm{adm}}$ to the 10th or 25th percentile of the indicator over the validation silent failures, which accepts a validation false-negative rate of about $10\%$ or $25\%$ in exchange for fewer false alarms, and a cost-optimal rule minimizes $20\,\mathrm{FN} + \mathrm{FP}$ on the validation set, that is, it prices a false acceptance at twenty needless fallbacks. All rules except the hybrid rule need at least one validation silent failure; the share of entries without one is given in parentheses.

\begin{table}[htbp]
\centering
\small
\caption{Recall / false-alarm rate on the test sets for five calibration rules with Eq.~\eqref{eq:eta-adm}, 40 validation queries, eight seeds, five dual dimensions pooled. Parentheses: share of entries for which the rule is undefined (no validation silent failure), which are excluded from recall and false-alarm rate. The percentile rules give recall close to $1 - q$ by construction.}
\label{tab:rules}
\input{tables/threshold_rules.tex}
\end{table}

The percentile rules buy a lower false-alarm rate at a recall that falls by about the percentile (recall $0.71$--$0.89$ for the 10th and 25th), so they move along the trade-off curve and do not improve it. The cost-optimal rule behaves like the $0.9\,\min$ rule, with recall $0.94$--$1.00$ and a lower false-alarm rate in some cases (B2 NMF $0.28$ against $0.31$) and a higher one in others (B2 gIS $0.62$ against $0.52$, B1 RID $0.21$ against $0.10$). A single exception is worth noting: for gIS on B2 the 10th-percentile rule keeps recall at $0.89$ and cuts the false-alarm rate from $0.52$ to $0.01$. No rule removes the trade-off; the choice of rule is a choice of operating point, and for a deployment it should be made with the cost of a missed failure against that of a fallback in view.

The false alarms and the true failures both send queries to the full-order solver, so the fallback share can erase the speedup. Table~\ref{tab:netspeed} combines the measured membrane speedups ($r_\lambda = 16$, Section~\ref{sec:basis_cost}) with the fallback share of the hybrid rule: the net speedup of a screened ROM that re-solves a fraction $f$ of the queries at full order is $1/(1/S + f)$, with the reduced solve always executed.

\begin{table}[htbp]
\centering
\small
\caption{Net speedup on the membrane with fallback. $f$: share of test queries flagged by the hybrid rule (40 validation queries, eight seeds, five dual dimensions pooled) and sent to the full-order solver; $S$: measured speedup of the unscreened ROM at $r_\lambda = 16$ (single-run sweep). Net speedup is $1/(1/S + f)$.}
\label{tab:netspeed}
\input{tables/net_speedup.tex}
\end{table}

"""
rep(anchor, block + anchor)

item_anchor = "\\item\n  \\textbf{Criterion 1 is not met in general, and its outcome depends on the seed.}"
item = r"""\item
  \textbf{Threshold sensitivity and fallback (Figures~\ref{fig:thrsens}, \ref{fig:tradeoff}, Tables~\ref{tab:rules}, \ref{tab:netspeed}).} No calibration rule removes the trade-off between missed failures and unnecessary fallbacks; percentile rules and a cost-optimal rule move the operating point along it. The cost of the screen can exceed its benefit on the membrane: with the hybrid rule the fallback share is $0.60$ (NMF), $0.83$ (CPG), $0.94$ (gIS) and $0.07$ (RID), because many B2 queries are genuine silent failures and not only false alarms. At $N = 324$ the unscreened constant-stiffness ROM is only $1.6$--$1.8$ times faster than the full-order solver for NMF, CPG and gIS, so the net speedup falls below one ($0.66$--$0.84$); with varying stiffness the unscreened speedup is about $5$ and the net speedup is $1.3$ (NMF), $1.0$ (CPG), $0.9$ (gIS) and $1.9$ (RID, whose unscreened speedup is $2.2$). On this benchmark the screened ROM pays off only with varying stiffness, clearly for RID ($1.9$) and marginally for NMF ($1.3$); a larger mesh, where the unscreened speedup is larger, would change this.
"""
rep(item_anchor, item + item_anchor)
io.open(p, "w", encoding="utf-8", newline="").write(s.replace(LF, CR + LF) if crlf else s)
print("ok")
