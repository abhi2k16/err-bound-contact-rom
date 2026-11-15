import sys
sys.path.insert(0, r"C:\Users\abhij\AppData\Local\Temp\claude\c5fb98397dc30c53\658d0914-ed8c-448a-92cd-0817b079476d\scratchpad")
from lang import Doc

d = Doc()

# ---------------------------------------------------------------- title
d.rep(r"\title{A Computable Admissibility Estimator and Break-Even Cost Model for Lagrange-Multiplier Reduced Contact Models}",
      r"\title{Lifted Error Bounds, Admissibility Screening and Break-Even Analysis for Lagrange-Multiplier Reduced-Order Contact Models}")

# ---------------------------------------------------------------- abstract
d.span(r"Projection-based reduced-order models (ROMs) for unilateral contact problems with Lagrange multipliers can satisfy",
       r"(for the latter with varying stiffness only at very small sizes).",
       r"""Projection-based reduced-order models (ROMs) for unilateral contact problems formulated with Lagrange multipliers may satisfy their projected Karush--Kuhn--Tucker (KKT) conditions to machine precision while producing lifted displacement--multiplier pairs $(\tilde{u}, \tilde{\lambda}) = (V u_r, U \lambda_r)$ that violate full-order physical admissibility. In our experiments the projected residuals were below $3 \times 10^{-14}$, whereas the relative energy-norm errors were of order one; we call such cases silent failures. Existing certified reduced-basis error bounds for variational inequalities have been developed under restrictive assumptions, including bijective constraint operators and dual feasibility enforced by construction. These assumptions do not hold for general rectangular constraint matrices $B(\mu) \in \mathbb{R}^{n_c \times N}$.

This work derives a residual-based \emph{a posteriori} energy-norm error bound for arbitrary lifted primal--dual pairs under rectangular constraint matrices. The bound consists of a bulk equilibrium residual ($\eta_{\mathrm{eq}}$) and three contact-interface contributions: multiplier sign violation ($\eta_{\mathrm{sign}}$), penetration ($\eta_{\mathrm{gap}}$) and complementarity defect ($\eta_{\mathrm{comp}}$). Because the equilibrium residual generally requires operations on the full displacement space, we also consider interface-based admissibility indicators ($\eta_{\mathrm{adm}}$ and a dimensionally consistent variant $\eta_{\mathrm{adm,v2}}$), which require $\mathcal{O}(n_c (r_u + r_\lambda))$ operations and no full-order bulk quantities. These indicators are used for online screening with a threshold calibrated on validation queries; they are not error bounds. Five dual-basis constructions (signed singular value decomposition with primal rank repair, non-negative matrix factorization, cone-projected greedy, greedy index selection and reduced integration domain restriction) are compared on identical snapshot sets, and a FLOP- and iteration-based cost model estimates the computational break-even point against sparse direct full-order contact solvers.

Numerical experiments on one-, two- and three-dimensional frictionless contact benchmarks show that the bound, evaluated with the exact multiplier norm, holds on every tested query. The interface indicators detect many but not all silent failures, and their performance depends on the dual basis, on the calibration data and, for queries outside the training range, can deteriorate. Predicted and measured break-even sizes agree within approximately a factor of two on the membrane and three-dimensional benchmarks, whereas absolute speedups are overpredicted at large sizes, and full-order fallback for rejected queries can remove the net speedup on the membrane.""")

# ---------------------------------------------------------------- introduction, para 2
d.span(r"A critical and under-appreciated failure mode in Lagrange-multiplier contact ROMs is \textbf{silent reduced-KKT convergence}.",
       r"(200 test queries per basis, dual dimensions 4--20 pooled).",
       r"""An important but insufficiently addressed failure mode of Lagrange-multiplier contact ROMs is \textbf{silent failure}: the reduced problem converges, but the reconstructed full-order solution is inadmissible. In this work, a silent failure denotes a query for which the projected KKT residual is below $10^{-12}$ while the relative energy-norm error of the lifted solution exceeds $5\%$. During online evaluation, a reduced contact solver enforces the projected KKT conditions only within the low-dimensional subspace. Consequently, these conditions may be satisfied to machine precision (residuals below $3 \times 10^{-14}$ in our experiments) while the lifted primal--dual pair $(\tilde{u}, \tilde{\lambda}) = (V u_r, U \lambda_r)$ violates full-order physical admissibility, with relative energy-norm errors of order one. When a standard proper orthogonal decomposition (POD) or singular value decomposition (SVD) is used to build the dual basis $U$, linear combinations of signed basis modes can produce negative components in the lifted multiplier vector ($\tilde{\lambda} < 0$, non-physical tensile tractions) together with interpenetration \cite{balajewicz2016projection, benaceur2020reduced}. In benchmark tests on parametrized membrane and discrete chain contact problems, a signed SVD dual basis satisfied reduced KKT tolerances below $10^{-14}$ while producing lifted multipliers as negative as $-1.6$ (chain) and $-0.83$ (membrane) and relative energy-norm errors of order one (median $0.96$ and $1.05$), whereas an NMF basis of identical dimension gave median errors of $0.4\%$ and $4\%$ (200 test queries per basis, dual dimensions 4--20 pooled).""")

# ---------------------------------------------------------------- research gaps
d.rep(r"\textbf{Positivity-preserving dual bases without lifted error bounds:}", r"\textbf{Positivity-preserving dual bases and the absence of lifted error bounds:}")
d.rep(r"nor have they been evaluated head-to-head on identical snapshot sets.", r"nor, to our knowledge, have they been compared systematically on identical snapshot sets.")
d.rep(r"\textbf{Restricted \emph{a posteriori} error estimation for reduced variational inequalities:}", r"\textbf{Limitations of existing \emph{a posteriori} error estimates:}")
d.rep(r"where reduced KKT satisfaction can mask catastrophic lifted inadmissibility.", r"where satisfaction of the projected KKT conditions can mask substantial inadmissibility of the lifted solution.")
d.rep(r"\textbf{Unquantified enrichment, conditioning, and rank-repair overhead:}", r"\textbf{Unquantified computational costs of enrichment and rank repair:}")
d.rep(r"the literature lacks an explicit break-even cost model predicting when rank-repair enrichment, active-set iteration growth, and dense reduced linear algebra eliminate the computational savings over a sparse full-order solver.",
      r"a systematic cost model that quantifies when rank-repair enrichment, active-set iteration growth and dense reduced linear algebra eliminate the computational savings over a sparse full-order solver remains limited.")

# ---------------------------------------------------------------- contributions
d.rep(r"\textbf{Computable Lifted Admissibility Estimator:}", r"\textbf{Lifted Error Bound and Admissibility Screening:}", count=2)
d.span(r"The three contact indicators evaluate in $\mathcal{O}(n_c (r_u + r_\lambda))$ operations without touching the bulk mesh, acting as an empirical online circuit breaker.",
       r"not a certificate).",
       r"""The three contact indicators can be evaluated in $\mathcal{O}(n_c (r_u + r_\lambda))$ operations without full-order bulk quantities and form the basis of an admissibility screening procedure (an online ``circuit breaker''). Three aspects are distinguished: the extension of the residual decomposition to arbitrary lifted primal--dual pairs, the resulting rigorous energy-norm error bound, and the interface-based screening indicator, which is an empirical detection mechanism and not a certified error estimator.""")
d.rep(r"With a validation-calibrated threshold it detects most, but not all, silent reduced-KKT failures in our experiments", r"With a threshold calibrated on validation queries, it detects most but not all silent failures in our experiments")
d.rep(r"evaluating lifted accuracy, estimator effectivity and, for NMF and CPG on the three-dimensional problem, out-of-range extrapolation.",
      r"evaluating lifted accuracy, bound effectivity and, for queries outside the training range, all five bases on B1 and B2 and NMF and CPG on B3.")
d.rep(r"identifying the regimes in which Lagrange-multiplier contact ROMs outperform full-order solvers;", r"identifying the regimes in which Lagrange-multiplier contact ROMs provide a computational advantage over full-order solvers;")
d.rep(r"for varying stiffness on B3 the break-even lies at very small sizes and is weakly determined, B1 could not be scored, and predicted speedup magnitudes overshoot measurements at large sizes.",
      r"for varying stiffness on B3 the break-even lies at small sizes and is weakly determined, B1 was not scored against the model, and predicted speedups exceed the measured ones at large sizes.")

# ---------------------------------------------------------------- positioning
d.span(r"The three contributions are of different kinds and should be weighed separately.",
       r"can be applied to any of them.",
       r"""The contributions of this work differ in mathematical and computational scope and are discussed separately. (i)~\emph{Extension of a known decomposition.} The splitting of the error into an equilibrium residual and sign, gap and complementarity defects is established for full-order contact \cite{hild2007residual, dabaghi2020posteriori, banz2024generic, banz2026posteriori}. We apply it to lifted reduced pairs that are neither Galerkin-orthogonal nor dually feasible, and we show that the sign defect vanishes identically for non-negative dual bases (Remark~2). (ii)~\emph{An error bound.} Theorem~\ref{thm:bound} is derived from the energy identity, the KKT conditions and standard inequalities; its contribution lies in its applicability to arbitrary lifted pairs under general rectangular constraint operators and not in a new proof technique. Relative to the duality-based bound for reduced variational inequalities \cite{zhang2014duality}, it requires neither a bijective constraint operator nor a dual approximation that is feasible by construction. Its practical value is limited by the multiplier norm it contains (Remark~1) and by the equilibrium residual, which is expensive in general. (iii)~\emph{A screening mechanism.} The interface indicators are an empirical detection tool with a calibrated threshold and not a certificate. Existing contact reduction strategies, namely positive dual bases \cite{balajewicz2016projection, benaceur2020reduced, sambataro2026model}, reduced integration domains \cite{fauque2018hybrid, leberre2022condition}, and empirical cubature and operator-inference approaches \cite{koller2021hyper, aresdeparga2023petrov, manvelyan2026inference}, aim at admissible reduced solutions but, to our knowledge, provide no per-query error check. Limitations of low-rank approximations in contact problems have been reported \cite{kollepara2022limitations}, and the authors of the greedy index selection postpone an error indicator to future work \cite{sambataro2026model}. The proposed bound and indicators are complementary to these constructions and can be applied to any of them.""")

# ---------------------------------------------------------------- research questions
d.span(r"The study is organized around four questions.", r"(Section~\ref{sec:experiments}.)" + "\n" + r"\end{enumerate}",
       r"""The study addresses four research questions.
\begin{enumerate}
\tightlist
\item[\textbf{RQ1}] Under what conditions can convergence of the projected KKT system coexist with inadmissibility of the lifted solution? (Section~\ref{sec:formulation}.)
\item[\textbf{RQ2}] Can a computable energy-norm error bound be established for arbitrary lifted primal--dual approximations under rectangular contact constraints, and which parts of it are inexpensive enough for online evaluation? (Section~\ref{sec:estimator}.)
\item[\textbf{RQ3}] How do the dual-basis construction, the primal rank repair and the online screening influence accuracy, admissibility and computational cost, and when does the reduced model outperform a sparse direct full-order solver? (Sections~\ref{sec:basis_cost} and~\ref{sec:experiments}.)
\item[\textbf{RQ4}] Under what conditions does the framework provide reliable predictions and net computational savings on the three benchmarks considered? (Section~\ref{sec:experiments}.)
\end{enumerate}""")
d.rep(r"Table~\ref{tab:contrib} relates the answers to the existing approaches. Entries in the last column are adaptations or combinations of cited work unless stated.",
      r"Table~\ref{tab:contrib} relates the contributions to existing approaches. Entries in the last column are adaptations or combinations of cited work unless stated otherwise.")

# ---------------------------------------------------------------- Section 2
d.rep(r"\section{Problem Formulation and Silent Reduced-KKT Failure}", r"\section{Problem Formulation and Silent Failure}")
d.rep(r"\subsection{Mechanisms of Silent Reduced-KKT Failure}", r"\subsection{Mechanisms of Silent Failure}")
d.span(r"Physically, a deformable body is pressed by a load against a rigid obstacle.", r"and the section shows how.",
       r"""In a frictionless unilateral contact problem, a deformable body interacts with a rigid obstacle under an applied load. The interaction is compressive: the normal contact pressure is non-negative, interpenetration is prohibited, and the pressure vanishes wherever a positive gap exists (complementarity). After discretization, the unknowns are the nodal displacements and the contact pressures, which are represented as Lagrange multipliers attached to candidate contact points; this gives the mixed problem below. A reduced model that approximates both fields may violate one or more of these conditions after reconstruction in the full-order space, and this section describes how.""")
d.rep(r"Equivalently, eliminating \(u(\mu) = K(\mu)^{-1}(f(\mu) - B(\mu)^T \lambda(\mu))\) via Schur complement condensation yields the dense full-order linear complementarity problem (LCP) for the multipliers",
      r"Eliminating the displacement, \(u(\mu) = K(\mu)^{-1}(f(\mu) - B(\mu)^T \lambda(\mu))\), from the equilibrium equation and substituting it into the contact conditions gives the following condensed (Schur-complement) linear complementarity problem (LCP) for the multiplier vector")
d.rep(r"A fundamental pitfall of the reduced formulation above is that enforcing the reduced KKT conditions to machine precision does \textbf{not} imply physical admissibility of the lifted pair",
      r"Enforcing the projected KKT conditions to machine precision does \textbf{not} imply physical admissibility of the lifted pair")
d.rep(r"the lifted solution can fail catastrophically through three coupled algebraic mechanisms:", r"the lifted solution can exhibit large errors through three coupled algebraic mechanisms:")
d.rep(r"Negative components in \(\tilde{\lambda}\) act as artificial adhesive tractions that pull the elastic body across the contact boundary to satisfy reduced equilibrium",
      r"Negative components in \(\tilde{\lambda}\) correspond to non-physical tensile tractions that pull the body across the contact boundary and can balance the projected equilibrium equation")
d.rep(r"Crucially, the reduced solver alone provides no warning that this breakdown has occurred, motivating the computable lifted admissibility estimator derived in Section~\ref{sec:estimator}.",
      r"The reduced solver itself gives no indication of this failure, which motivates the lifted error bound and admissibility indicators derived in Section~\ref{sec:estimator}.")
d.rep(r"The SVD-lifted multiplier is adhesive (negative) inside the patch", r"The SVD-lifted multiplier is negative (tensile) inside the patch")

# ---------------------------------------------------------------- simple global replacements
d.rep(r"with the online circuit breaker \(\eta_{\text{adm}} \le \tau_{\text{adm}}\)", r"with the online admissibility screening \(\eta_{\text{adm}} \le \tau_{\text{adm}}\)")
d.rep(r"the circuit breaker accepted no unsafe answer", r"the admissibility screening accepted no unsafe answer")
d.rep(r"integrating the online circuit breaker \(\eta_{\text{adm}}\)", r"integrating the online admissibility screening \(\eta_{\text{adm}}\)")
d.rep(r"Table~\ref{tab:e3} gives the head-to-head comparison of the five bases", r"Table~\ref{tab:e3} gives the comparative assessment of the five bases")
d.save()
print("ok")
