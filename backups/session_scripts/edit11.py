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


anchor = r"\paragraph{Ablation of the framework components.}"
block = r"""\paragraph{Out-of-range queries and contact-patch translation.}
The calibration and comparison studies query the reduced models inside the training box. We therefore also queried the five bases ($r_\lambda = 16$, rank repair, five seeds, 40 queries per set and seed) outside it: loads $1.3$ to $2$ times the training maximum on both benchmarks, a steeper obstacle on B1 (curvature $1.25$ to $1.7$ times the training maximum), and an indenter centre outside the training box on B2, which moves the contact patch (each coordinate in $[0.20, 0.30]$ or $[0.70, 0.80]$ against $[0.35, 0.65]$ in training). The screen uses the hybrid threshold calibrated on the in-range validation queries. Table~\ref{tab:extrap} gives the share of queries above $5\%$ error, the missed silent failures among the silent failures, and the validity of the online bound with the offline constant $C_\lambda$.

\begin{table}[htbp]
\centering
\small
\caption{Out-of-range queries. Cells: share of queries above $5\%$ error in \% (silent failures missed by the screen / silent failures), five seeds pooled (200 queries per cell). The screen uses the hybrid threshold of Eq.~\eqref{eq:eta-adm} calibrated on the 40 in-range validation queries of each seed. Last column: smallest share, over the five bases, of queries on which the bound with the offline constant $C_\lambda$ is at least the true error; the bound with the exact multiplier norm held on every query in every set.}
\label{tab:extrap}
\input{tables/extrapolation.tex}
\end{table}

\paragraph{Mesh-size study.}
To see whether the findings depend on the mesh, we repeated the audit on five mesh sizes of B1 (chain length $30$ to $480$ with $n_c = N/2$) and B2 (membrane grid $12$ to $40$, $N = 144$ to $1600$, with the contact sub-grid and the load scaled as in the break-even sweep), with three seeds, 40 test queries, the dual dimension fixed at $r_\lambda = 16$ and the primal dimension at $20$ (Table~\ref{tab:mesh}). The speedup is that of the unscreened ROM against the full-order solver in the constant-stiffness regime. We did not repeat this study on B3, whose break-even sweep (Table~\ref{tab:timing}) is the only mesh study there.

\begin{table}[htbp]
\centering
\small
\caption{Mesh-size study (three seeds, 40 test queries each, $r_\lambda = 16$): share of queries above $5\%$ error for SVD, NMF (with median error in brackets) and gIS, and the measured speedup $S$ of the unscreened ROM for NMF and gIS (constant stiffness).}
\label{tab:mesh}
\input{tables/meshsize.tex}
\end{table}

"""
rep(anchor, block + anchor)

item_anchor = "\\item\n  \\textbf{Ablation (Table"
item = r"""\item
  \textbf{Out-of-range queries (Table~\ref{tab:extrap}).} Every basis degrades outside the training box, and the screen does not always follow. On B1 a load $1.3$ to $2$ times the training maximum raises the share above $5\%$ from $0\%$ to $26\%$ (NMF) and $37\%$ (CPG), to $100\%$ for gIS and to $70\%$ for RID; on B2 the in-range share of $26\%$ (NMF) and $44\%$ (CPG) goes to $99\%$ and $100\%$, whereas RID, which had no silent failure in range, reaches $45\%$. Contact-patch translation on B2 is nearly as damaging for NMF, CPG and gIS ($95$--$96\%$ above $5\%$) but not for RID ($0\%$ in the tested shifts). A steeper obstacle on B1 leaves NMF, CPG and gIS accurate ($0$--$1\%$) and affects RID ($16\%$). The screen flagged almost all of these failures where it had a calibrated threshold (missed $0$ of $200$ to $2$ of $192$ for NMF, CPG and gIS in all sets except B2 in range). The exception that matters is RID on B2: its validation set contained no silent failure, so the hybrid rule used the relative default, and the screen missed $78$ of the $90$ silent failures under load extrapolation ($20$ of $140$ on B1). This is the failure mode that Section~\ref{sec:conclusion} warns about, now shown on a benchmark: a screen calibrated only in range can be blind outside it. The bound with the exact multiplier norm held on every query, but the online bound with the offline constant $C_\lambda$ fell below the true error on $1.5\%$ of the extrapolated queries on B1 (SVD, load) and on up to $3.5\%$ on B2 (RID, centre shift), which is the limitation stated for it before and is now measured.
\item
  \textbf{Mesh size (Table~\ref{tab:mesh}).} The bound with the exact multiplier norm holds at all meshes; the only value below $100\%$ (gIS at $N = 30$ on B1) occurs where the median error is of the order of $10^{-15}$, which we attribute to rounding without having examined individual queries. At a fixed dual dimension the silent-failure share grows with the mesh whenever the dual space no longer covers the contact set: for gIS on B1 from $0\%$ at $N = 30$ (where $n_c = 15 < r_\lambda$) to $89\%$ at $N = 480$, for gIS on B2 from $69\%$ to $97\%$, and for NMF on B2 from $12\%$ to $32\%$; NMF on B1 stays below $5\%$ at every size. Accuracy at fixed $r_\lambda$ therefore does not scale with $N$, and the dual dimension must grow with the contact set, which we did not study. The speedup of NMF grows with $N$ on both benchmarks: on B1 from $0.21$ at $N = 30$ to $2.9$ at $N = 480$ (break-even at about $N = 220$ by interpolation) and on B2 from $0.81$ at $N = 144$ to $45$ at $N = 1600$ (break-even at about $N = 180$, against $154$ in the break-even sweep of Table~\ref{tab:timing}). The first break-even size on B1 is new, since B1 had not been swept.
"""
rep(item_anchor, item + item_anchor)

# limitations: remove statements that are now done
rep(r"Contact-patch translation on the membrane, an in-sample validation set and the time-dependent extension (E6) were not run;",
    r"An in-sample validation set and the time-dependent extension (E6) were not run, and out-of-range queries were tested only for one range of loads, obstacle shapes and indenter positions on B1 and B2 (not on B3 for the five bases);")
io.open(p, "w", encoding="utf-8", newline="").write(s.replace(LF, CR + LF) if crlf else s)
print("ok")
