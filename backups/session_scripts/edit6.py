import io
p = r"C:\Users\abhij\Desktop\Physical_ML\proposed-research\reduced-contact-reliability\manuscript_master_v2.tex"
CR, LF = chr(13), chr(10)
raw = io.open(p, encoding="utf-8", newline="").read()
crlf = (CR + LF) in raw
s = raw.replace(CR + LF, LF)

def rep(old, new):
    global s
    assert s.count(old) == 1, old[:70]
    s = s.replace(old, new, 1)

# 1. paragraph + table before the findings list (anchor: the matched-comparison bullet)
head = r"\begin{itemize}" + LF + r"\tightlist" + LF + r"\item" + LF
anchor = head + r"  \textbf{Matched comparison (Tables~\ref{tab:matcheddim}"
block = r"""\paragraph{Repeated timings and variability.}
The break-even sizes of Section~\ref{sec:basis_cost} come from single-thread wall-clock times, so we repeated the membrane (B2) and two-body (B3) size sweeps with ten timing repeats per query and ten queries per size, keeping every raw repeat. Point estimates use the estimator of the main tables (minimum over repeats, median over queries). Intervals are $95\%$ percentile bootstrap intervals that resample queries and repeats; for the break-even size $N^*$ every size is resampled independently and $N^*$ is recomputed. The measurement noise is moderate: the coefficient of variation across repeats of one query is $0.02$--$0.21$ over all sizes and solvers (Table~\ref{tab:timing}), and the spread between queries is larger, up to about $50\%$ for the reduced solve at some sizes, because the iteration counts differ.

\begin{table}[htbp]
\centering
\small
\caption{Predicted and measured break-even sizes $N^*$ and the speedup at the largest size (B2: $N = 2304$; B3: $N = 16\,653$), with $95\%$ bootstrap intervals. $N^*_{\mathrm{E4}}$ is the value of the earlier single-run sweep used in the other tables, $N^*_{\mathrm{meas}}$ the value of the ten-repeat run. pred/meas is $N^*_{\mathrm{pred}}/N^*_{\mathrm{meas}}$. $^{*}$: $N^*$ is defined in only $65\%$ (B2) and $79\%$ (B3) of the bootstrap resamples, because in the others the speedup exceeds one over the whole range or does not stay above one; the interval is conditional on that.}
\label{tab:timing}
\input{tables/timing_table.tex}
\end{table}

"""
item = r"""\item
  \textbf{Repeated timings (Table~\ref{tab:timing}).} The measured break-even sizes are reproducible within the bootstrap interval in most cases, but the interval understates the variability between runs. For constant $K$ on B2, $N^*_{\mathrm{meas}} = 143$ [98, 164] against $154$ in the single-run sweep and $119$ predicted; on B3, $722$ [644, 937] against $1006$ in the single-run sweep, a value outside the interval, and $779$ predicted. For varying $K$ the values are $32$ [28, 35] (B2) and $67$ [64, 72] (B3) against $64$ and $55$ predicted, and the B2 ratio predicted/measured is $1.97$, against $1.72$ in the single run: the agreement within a factor of two is marginal there, and the interval reaches about $2.3$. The varying-$K$ crossover on B3 is weakly determined because the measured speedup is not monotone at the smallest sizes. The speedup magnitudes are not reproduced: at the largest sizes the model predicts $1730$ against a measured $179$ [74, 265] (B2, constant $K$), $628$ against $287$ (B2, varying $K$), $730$ against $48$ (B3, constant $K$) and $1138$ against $411$ (B3, varying $K$), and at small sizes the measured speedup is above the prediction. The model locates the break-even size to within roughly a factor of two of a single noisy measurement; it does not predict the speedup.
"""
rep(anchor, block + head + item + r"""\item
  \textbf{Matched comparison (Tables~\ref{tab:matcheddim}""")

# 2. conclusion sentence
rep(r"""the latter crossover lies below \(N = 100\) and is weakly determined, B1 was not swept, and predicted speedups overshoot measured ones at large sizes.""",
    r"""the latter crossover lies below \(N = 100\) and is weakly determined, B1 was not swept, and predicted speedups overshoot measured ones at large sizes by up to an order of magnitude. A ten-repeat re-run (Table~\ref{tab:timing}) reproduces the break-even sizes within bootstrap intervals in most cases, but the varying-stiffness ratio on the membrane rises to 1.97 and the single-run B3 value lies outside the interval of the re-run, so the factor-of-two agreement should be read as approximate.""")

io.open(p, "w", encoding="utf-8", newline="").write(s.replace(LF, CR + LF) if crlf else s)
print("ok")
