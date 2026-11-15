import sys
sys.path.insert(0, r"C:\Users\abhij\AppData\Local\Temp\claude\c5fb98397dc30c53\658d0914-ed8c-448a-92cd-0817b079476d\scratchpad")
from lang import Doc

d = Doc()

# ================================================================ Section 3
d.rep(r"\section{Computable \emph{A Posteriori} Lifted Admissibility Estimator}",
      r"\section{\emph{A Posteriori} Error Bound and Online Admissibility Indicators for Lifted Solutions}")
d.span(r"To detect silent reduced-KKT convergence without solving the full-order contact problem, we derive", r"\cite{banz2024generic, banz2026posteriori}.",
       r"""To detect silent failures without solving the full-order contact problem, we derive a residual-based \emph{a posteriori} error bound for arbitrary lifted primal--dual pairs \((\tilde{u}, \tilde{\lambda}) = (V u_r, U \lambda_r)\). Existing reduced-basis error estimates for variational inequalities require a bijective constraint operator and strictly feasible dual spaces \cite{zhang2014duality}. The present analysis applies to general rectangular constraint matrices \(B(\mu) \in \mathbb{R}^{n_c \times N}\) and makes no assumption on how \((u_r, \lambda_r)\) is computed or on whether \(U\) preserves non-negativity \cite{banz2024generic, banz2026posteriori}.""")
d.rep(r"The two parts are kept apart on purpose.",
      r"The two parts are kept apart on purpose. Terminology is used as follows. The \emph{error bound} is the mathematically derived upper bound of Section~\ref{certified-energy-norm-upper-bound}. The \emph{admissibility indicators} are the individual sign, gap and complementarity measures, and the \emph{online admissibility indicator} is their combination in Section~\ref{computational-complexity-and-online-admissibility-indicator}. The \emph{screening threshold} is the calibrated value against which the online indicator is compared, and a \emph{fallback solve} is the full-order solve performed for a rejected prediction. The terms estimator and certificate are not used for the online indicator.")
d.rep(r"into its physically admissible (compressive) and inadmissible (adhesive) parts:",
      r"into its physically admissible (compressive) part and its non-physical (negative, tensile) part:")
d.rep(r"collects the non-physical adhesive components", r"collects the non-physical negative (tensile) components")
d.rep(r"is the magnitude of the adhesive part", r"is the norm of the negative multiplier component")
d.rep(r"\(\eta_{\mathrm{sign}}\) measures the Euclidean norm of non-physical adhesive contact tractions (\(\tilde{\lambda}_i < 0\))",
      r"\(\eta_{\mathrm{sign}}\) quantifies the violation of multiplier non-negativity as the Euclidean norm of the negative component (\(\tilde{\lambda}_i < 0\))")
d.rep(r"is rigorously bounded from above by these four quantities and the multiplier norm", r"is bounded from above by these four quantities and the multiplier norm")
d.rep(r"Let \((u, \lambda) \in \mathbb{R}^N \times \mathbb{R}^{n_c}\) be the exact full-order KKT solution of the discrete contact problem, and let",
      r"Let \(K\) be symmetric positive definite, let the admissible set \(\{v : B v \le c\}\) be non-empty, let \((u, \lambda) \in \mathbb{R}^N \times \mathbb{R}^{n_c}\) be a KKT solution of the discrete contact problem~\eqref{eq:kkt}, and let")
d.rep(r"The third term is bounded immediately by", r"The third term is bounded by")
d.rep(r"Three different statements are made in this paper and they carry different guarantees.", r"Three statements are made in this paper, and they carry different guarantees.")
d.rep(r"It can be called trustworthy only to the extent of the empirical recall reported in Section~\ref{sec:experiments}, and the screening indicator is meant to trigger a fallback, not to replace the bound of item (a).",
      r"The reliability of an accepted prediction is therefore limited to the empirical recall reported in Section~\ref{sec:experiments}. The screening indicator is intended to trigger a fallback and does not replace the bound of item (a).")
d.rep(r"Whenever an indicator exceeds a user-specified tolerance \(\tau_{\mathrm{adm}}\), the reduced prediction is flagged as physically inadmissible---even if the reduced KKT residual is at machine precision---triggering either adaptive basis enrichment or fallback to the full-order solver.",
      r"When an indicator exceeds the screening threshold \(\tau_{\mathrm{adm}}\), the reduced prediction is rejected, even if the projected KKT residual is at machine precision, and either the basis is enriched or the full-order solver is invoked as a fallback.")

# ================================================================ Section 4
d.span(r"A common misconception in reduced-order modeling for contact is that", r"erode the computational advantage of the ROM:",
       r"""A reduction in the number of degrees of freedom does not necessarily translate into a reduction in the online computational cost of contact problems. Two principal factors can reduce the computational advantage of the ROM:""")
d.rep(r"To quantify the computational crossover, we model", r"To estimate the computational break-even point, we model")
d.rep(r"\emph{(Note: When \(K\) is parameter-independent, the first two terms are amortized offline \cite{pawar2025efficient}, leaving only",
      r"When \(K\) is parameter-independent, the first two terms can be evaluated offline and amortized over subsequent parameter queries \cite{pawar2025efficient}, leaving only")
d.rep(r"all three terms execute online.)}", r"all three terms are evaluated online.")
d.rep(r"the reduced contact model achieves computational break-even (\(S > 1\)) if and only if", r"the reduced model becomes computationally advantageous (\(S > 1\)) if and only if")
d.rep(r"the cubic dense cost \(r_u^3\) immediately exceeds \(2 k_{\text{FOM}} m_{\text{act}}^3\), causing \(S < 1\).",
      r"the cubic scaling of the dense reduced operations, \(r_u^3\), may exceed \(2 k_{\text{FOM}} m_{\text{act}}^3\), giving \(S < 1\).")
d.rep(r"The gIS choice by snapshot energy is the greedy choice that minimizes the neglected multiplier energy, because the columns of \(U\) are orthonormal.",
      r"The gIS strategy selects contact indices according to their accumulated multiplier energy. Because the columns of \(U\) are orthonormal, each greedy choice minimizes the multiplier energy that remains neglected.")
d.rep(r"The repair loop stops when \(\kappa \le \kappa_{\max}\), after at most \(r_\lambda\) added vectors, or when the supremizer already lies in the span of \(V\); in the last two cases the target is not reached and the algorithm returns the enriched basis with the condition number it achieved.",
      r"The enrichment terminates when the condition-number limit \(\kappa_{\max}\) is reached, when \(r_\lambda\) vectors have been added, or when the computed supremizer already lies in the span of \(V\). In the last two cases the limit is not attained, and the algorithm returns the enriched basis together with the condition number it achieved.")
d.rep(r"The model counts operations, not memory traffic or storage.", r"The model accounts for floating-point operations but not for memory traffic or storage requirements.")
d.rep(r"The $N^\alpha$ form quoted above for three dimensions was not used.",
      r"Although an $N^\alpha$ scaling for three-dimensional sparse factorization is discussed above, it was not used in the numerical break-even calculations.")
d.rep(r"Only the location of the break-even size, not the magnitude of the speedup, should be read from the model, and both are specific to this code and machine.",
      r"Accordingly, the model should be interpreted as an estimate of the break-even problem size and not as a quantitative predictor of the absolute speedup, and both are specific to this code and machine.")
d.save()
print("ok")
