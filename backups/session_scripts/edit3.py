import io
p = r"C:\Users\abhij\Desktop\Physical_ML\proposed-research\reduced-contact-reliability\manuscript_master_v2.tex"
CR, LF = chr(13), chr(10)
raw = io.open(p, encoding="utf-8", newline="").read()
crlf = (CR + LF) in raw
s = raw.replace(CR + LF, LF)

head = r"\begin{itemize}" + LF + r"\tightlist" + LF + r"\item" + LF
anchor = head + r"  \textbf{Criterion 1 is not met in general, and its outcome depends on the seed.}"
assert s.count(anchor) == 1

block = r"""\paragraph{Threshold calibration study.}
To quantify how reliable the calibrated screen is, we repeated the E3 audit for eight seeds (42--49) on B1 and B2 with a pool of 200 validation queries and 100 test queries per seed. Smaller validation sets are random subsets of the pool (20 draws each), so the validation-size sweep needs no extra solves. Positives are silent failures (reduced KKT residual below $10^{-12}$, relative energy error above $5\%$) and negatives are accurate queries (error at most $5\%$), as in the main tables; counts are pooled over the five dual dimensions. Three threshold rules are compared for Eq.~\eqref{eq:eta-adm}. The \emph{$0.9\,\min$ rule} is the one used so far, $\tau_{\mathrm{adm}} = 0.9\min(\eta_{\mathrm{adm}})$ over validation silent failures; it is undefined when the validation set contains none. The \emph{relative rule} $\tau_{\mathrm{adm}} = 0.05\,\overline{\|u\|_K}$ (training mean) needs no validation data. The \emph{hybrid rule} takes the smaller of the two when the validation set contains a silent failure and the relative value otherwise, so it is always defined. Uncertainty is a $95\%$ percentile bootstrap interval resampling whole seeds, the independent units; intervals from pooled query counts would be too narrow because queries, dual dimensions and subsets are dependent.

\begin{table}[htbp]
\centering
\small
\caption{Detection of silent failures with Eq.~\eqref{eq:eta-adm}, 40 validation queries, eight seeds, five dual dimensions pooled; brackets are $95\%$ bootstrap intervals over seeds. Recall is the share of silent failures flagged; the false-alarm rate is the share of accurate queries flagged. ``Uncal.'' is the share of (seed, dual dimension) entries for which the $0.9\,\min$ rule has no threshold because the validation set holds no silent failure; the recall of that rule is computed over the remaining entries only. Signed SVD is omitted: it has no accurate queries and its recall is at least $0.99$.}
\label{tab:calibration}
\input{tables/calibration_table_40.tex}
\end{table}

\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{fig_calibration.pdf}
\caption{Recall of the calibrated screen, Eq.~\eqref{eq:eta-adm}, against the number of validation queries (eight seeds; bands are $95\%$ bootstrap intervals over seeds). The $0.9\,\min$ rule is scored on calibrated entries only, so it is optimistic when many entries are uncalibrated (Table~\ref{tab:calibration}, for example B1 NMF, 0.74).}
\label{fig:calibration}
\end{figure}

"""

item = r"""\item
  \textbf{Calibration study (Table~\ref{tab:calibration}, Figure~\ref{fig:calibration}).} Recall rises with the number of validation queries: for the $0.9\,\min$ rule it goes from $0.73$--$0.91$ at ten queries on B1 (NMF to RID) to $0.98$--$1.00$ at 200, and it stays below one for several bases even at 200 (for example, on B1 gIS misses 13 of 2452 silent failures and on B2 CPG misses 14 of 2651). With 40 queries the $0.9\,\min$ rule has no threshold for 74\% (NMF), 46\% (CPG) and 36\% (RID) of the B1 entries and for 82\% of the B2 RID entries. The hybrid rule removes this gap and gives recall $0.93$--$1.00$ on B1 at false-alarm rates of $0.02$--$0.04$. On B2 the false-alarm rate is high whatever the rule (hybrid, 40 queries: NMF $0.31$, CPG $0.53$, gIS $0.52$, RID $0.02$), so on that benchmark the screen buys recall at the price of many unnecessary fallbacks. The relative rule alone is not adequate: with Eq.~\eqref{eq:eta-adm} on B2 it misses $32$--$44\%$ of silent failures for NMF, CPG and RID, and with $\eta_{\mathrm{adm,v2}}$ it flags every accurate query for NMF and CPG there, because the factor $C_\lambda$ makes that indicator conservative. If the validation set holds no silent failure, no data-driven threshold exists and the screen has no empirical basis for that basis and dual dimension; the hybrid rule then falls back on the relative value, which is a default and not a calibrated threshold.
"""

s = s.replace(anchor, block + head + item + r"""\item
  \textbf{Criterion 1 is not met in general, and its outcome depends on the seed.}""", 1)
io.open(p, "w", encoding="utf-8", newline="").write(s.replace(LF, CR + LF) if crlf else s)
print("ok")
