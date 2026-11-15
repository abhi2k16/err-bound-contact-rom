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


# ---- A. Research questions and contributions (end of introduction, before the workflow figure) -----------------
fig = s.index(r"\begin{figure}[htbp]")
rq = r"""\subsection*{Research Questions and Contributions}

The study is organized around four questions.
\begin{enumerate}
\tightlist
\item \emph{What is the limitation?} Why can a reduced contact solution satisfy its reduced KKT equations to machine precision and still give an inadmissible lifted solution? (Section~\ref{sec:formulation}.)
\item \emph{What is the mathematical contribution?} What can be said, for an arbitrary lifted primal--dual pair, about the energy-norm error and its interface components, and which part of that is cheap enough to evaluate online? (Section~\ref{sec:estimator}.)
\item \emph{What is the computational contribution?} How do the choice of dual basis, the primal rank repair and the online screen interact in accuracy and cost, and when does the reduced model beat a sparse direct full-order solver? (Sections~\ref{sec:basis_cost} and~\ref{sec:experiments}.)
\item \emph{What is demonstrated?} Under which conditions are the predictions reliable and the savings real, on the three benchmarks used? (Section~\ref{sec:experiments}.)
\end{enumerate}
Table~\ref{tab:contrib} relates the answers to the existing approaches. Entries in the last column are adaptations or combinations of cited work unless stated.

\begin{table}[htbp]
\centering
\small
\caption{Existing approaches, the limitation addressed and the contribution of this paper.}
\label{tab:contrib}
\begin{tabular}{@{}p{0.22\textwidth}p{0.26\textwidth}p{0.44\textwidth}@{}}
\toprule
Existing approach & Limitation addressed & Contribution of this paper \\
\midrule
Primal contact ROM, Lagrange-multiplier reduced models \cite{balajewicz2016projection, benaceur2020reduced, leberre2022condition} & Lifted contact admissibility can be lost without any warning from the reduced solver & Documents the silent failure with signed and non-negative bases on three benchmarks, and quantifies it with the bound and its sign, gap and complementarity components. \\
Dual-basis construction \cite{balajewicz2016projection, benaceur2020reduced, sambataro2026model, fauque2018hybrid} & The five bases differ in accuracy, effective size, primal repair and cost, and are rarely compared on common snapshots & Controlled comparison of five bases on identical snapshots at matched effective dimension, cost and accuracy, with rank-repair overhead included; no new basis is proposed. \\
\emph{A posteriori} error estimation for contact \cite{hild2007residual, banz2024generic, zhang2014duality} & Full-order residual evaluation is expensive, and existing reduced-basis bounds need a bijective constraint operator & Extension of the residual decomposition to arbitrary lifted pairs under rectangular $B$ (no new proof technique); a cheap interface screen that is empirical, with a calibration study of its failure modes. \\
Acceleration of contact ROMs \cite{koller2021hyper, aresdeparga2023petrov, pawar2025efficient} & Reduced cost does not imply a net speedup over a sparse solver & Operation-count break-even model with stated assumptions (direct solvers), tested against repeated timings, and the effect of fallback on the net speedup. \\
\bottomrule
\end{tabular}
\end{table}

"""
s = s[:fig] + rq + s[fig:]

# ---- B. Roadmaps and physical introduction ----------------------------------------------------------------------
rep(r"\label{sec:formulation}" + LF,
    r"\label{sec:formulation}" + LF + LF +
    r"""Physically, a deformable body is pressed by a load against a rigid obstacle. Where the two touch, the obstacle can only push on the body: the contact pressure is non-negative, the gap between the surfaces cannot become negative (no penetration), and pressure acts only where the gap is closed (complementarity). After discretization, the displacement unknowns are the nodal displacements and the contact pressures are Lagrange multipliers attached to candidate contact points, which gives the mixed problem below. A reduced model that approximates both displacements and pressures can break any of the three conditions, and the section shows how.
""")
rep(r"\label{sec:estimator}" + LF,
    r"\label{sec:estimator}" + LF + LF +
    r"""\emph{Roadmap.} Section~\ref{norms-multiplier-splitting-and-estimator-components} defines four residual indicators of a lifted pair. Section~\ref{certified-energy-norm-upper-bound} derives a rigorous energy-norm bound from them; it needs the multiplier norm and an equilibrium residual that is a full-order quantity in general, so it is a validation tool. Section~\ref{computational-complexity-and-online-admissibility-indicator} then builds the inexpensive interface-only indicators from the three terms that do not touch the bulk; these are an empirical screen and carry no error guarantee. The two parts are kept apart on purpose.
""")
rep(r"\label{sec:basis_cost}" + LF,
    r"\label{sec:basis_cost}" + LF + LF +
    r"""\emph{Why the dual basis matters for both admissibility and cost.} The dual basis $U$ decides whether the lifted multipliers can be negative (a signed basis can, a non-negative one cannot), how many contact rows the reduced model carries, and how much the primal space must be enlarged to keep the reduced Schur complement well conditioned. The first effect governs admissibility and accuracy, the other two govern the online cost through the dense reduced algebra, so the two cannot be assessed separately.
""")
rep(r"To systematically evaluate the \emph{a posteriori} estimator (Section 3)",
    r"The section first describes the three benchmarks and their parameter ranges, then the experiments and acceptance criteria, then the results, with observations kept apart from their interpretation, which is collected in Section~\ref{sec:conclusion}. To systematically evaluate the \emph{a posteriori} estimator (Section 3)")

# ---- C. Effective dimension table in the matched paragraph --------------------------------------------------------
rep(r"""\begin{table}[htbp]
\centering
\small
\caption{Median relative energy error at matched effective dual dimension""",
    r"""Table~\ref{tab:effdim} lists, at nominal $r_\lambda = 16$, the effective dual dimension, the repair vectors, the enriched primal dimension and the offline and online costs of the five methods. The offline cost is that of the dual basis and the rank repair, measured in our implementation.

\begin{table}[htbp]
\centering
\small
\caption{Nominal and effective dimensions, primal repair, and offline and online cost of the five methods at nominal $r_\lambda = 16$ (100 test queries; online time is the median reduced-solve time).}
\label{tab:effdim}
\input{tables/basis_effective.tex}
\end{table}

\begin{table}[htbp]
\centering
\small
\caption{Median relative energy error at matched effective dual dimension""")

# ---- D. Discussion (interpretation) before the practical guidelines -----------------------------------------------------
disc = r"""\subsection{Discussion}\label{sec:discussion}

The observations are in Section~\ref{sec:experiments}. Here we interpret them. The explanations are plausible readings of the data and were not isolated by dedicated experiments.

\emph{When does each dual basis work well?} Signed SVD fails because nothing prevents negative lifted multipliers, so $\eta_{\mathrm{sign}}$ is large and the bound with it. NMF and CPG remove the sign defect by construction, which is why their errors are of the order of a few percent or less, but their columns do not enforce complementarity at individual contact points, and on the membrane the error saturates at about $2$--$3\%$ as the dimension grows. gIS keeps unit vectors at selected contact indices, so its errors are large until the selected set captures the points that are in contact; on B1 the median error falls from $0.05$ at $r_\lambda = 16$ to $3\times10^{-6}$ at $24$, a sharp transition rather than a gradual one. RID retains all contact rows of the selected elements, which gives the highest accuracy at high dimension ($10^{-3}$ or below) but at the price of an effective dual dimension of 23 (B1) and 98 (B2) at nominal 16 and of up to 78 repair vectors.

\emph{When is the framework computationally beneficial?} Our measured break-even sizes are about $N = 150$ (membrane, constant stiffness), $30$--$40$ (membrane, varying stiffness), about $700$--$1000$ (two-body, constant stiffness) and $55$--$90$ (two-body, varying stiffness), for direct full-order solvers and the code and machine used. With constant stiffness the saving only appears at larger sizes because the full-order cost is then the active-set solve alone. The savings are reduced by the fallback: on the membrane at $N = 324$ the hybrid rule sends between $7\%$ and $94\%$ of the queries to the full-order solver depending on the basis, and the net speedup of the screened ROM is below one for constant stiffness. The offline cost of the dual basis is small in absolute terms at these sizes (1--2\,ms for SVD, about 140--190\,ms for NMF), and primal repair is the dominant effect on the online cost only for RID.

\emph{Limitations.} The interface screen can miss inadmissible solutions, and some bases give high false-alarm rates; the computational advantage depends on the benchmark and the problem size; the time-dependent extension was not demonstrated; and changing contact geometries and contact-patch translation were not studied. These are developed in Section~\ref{limitations-and-future-work}.

"""
anchor = r"\hypertarget{practical-guidelines-for-contact-rom-deployment}{%"
i = s.index(anchor)
s = s[:i] + disc + s[i:]
io.open(p, "w", encoding="utf-8", newline="").write(s.replace(LF, CR + LF) if crlf else s)
print("ok")
