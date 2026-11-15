import sys
sys.path.insert(0, r"C:\Users\abhij\AppData\Local\Temp\claude\c5fb98397dc30c53\658d0914-ed8c-448a-92cd-0817b079476d\scratchpad")
from lang import Doc

d = Doc()

# ================================================================ Section 5: introduction and campaign
d.span(r"The section first describes the three benchmarks and their parameter ranges", r"leberre2022condition}.",
       r"""This section describes the benchmark problems and their parameter ranges, outlines the experimental campaign and the evaluation criteria, and presents the numerical results. Their interpretation is discussed in Section~\ref{sec:conclusion}. The experiments evaluate the \emph{a posteriori} error bound and the admissibility indicators (Section 3), compare the five dual-basis constructions (Section 4.1) and validate the break-even cost model (Section 4.2) in a five-part suite (E1--E5) on one-, two- and three-dimensional unilateral contact problems. All reduced-order models are evaluated on identical snapshot sets \(\{(u(\mu_i), \lambda(\mu_i))\}_{i=1}^{n_s}\) and solved with a common active-set LCP solver (the projected KKT residual stayed below \(3 \times 10^{-14}\) in all experiments) \cite{balajewicz2016projection, benaceur2020reduced, leberre2022condition}.""")
d.rep(r"\subsection{Experimental Campaign (E1--E6)}", r"\subsection{Experimental Campaign (E1--E5)}")
d.rep(r"\textbf{Experiment E1: Silent Reduced-KKT Failure Audit (B1, B2)}", r"\textbf{Experiment E1: Silent Failure Audit (B1, B2)}")
d.rep(r"masks catastrophic lifted error.", r"masks substantial error in the lifted solution.")
d.rep(r"\textbf{Experiment E2: Estimator Effectivity and Extrapolation Stress Test (B1, B2, B3)}", r"\textbf{Experiment E2: Bound Effectivity and Extrapolation Stress Test (B1, B2, B3)}")
d.rep(r"Contact-patch translation on the membrane and an in-sample validation set were not run.",
      r"Out-of-range queries on B1 and B2 (loads, obstacle shape and indenter position, including contact-patch translation) are reported separately below, and an in-sample validation set was not run.")
d.rep(r"\textbf{Experiment E3: Head-to-Head Five-Basis Benchmark (B1, B2)}", r"\textbf{Experiment E3: Comparative Five-Basis Benchmark (B1, B2)}")
d.rep(r"\textbf{Experiment E5: Large-Scale 3D Contact Certification (B3)}", r"\textbf{Experiment E5: Admissibility Screening on Large-Scale 3D Contact (B3)}")
d.span("\\item\n  \\textbf{Experiment E6", r"\end{itemize}",
       r"""\end{itemize}

An extension to time-dependent contact dynamics was considered as a further experiment but was not performed in the present study; it is listed as future work.""")

# ================================================================ Section 5: criteria
d.span(r"Every test query with \(\Vert R_{\text{red}}\Vert_\infty \le 10^{-12}\) whose true relative energy-norm error exceeds \(5\%\)", r"must satisfy \(\eta_{\text{adm}} > \tau_{\text{adm}}\).",
       r"""A test query is classified as a silent failure when the projected KKT residual \(\Vert R_{\text{red}}\Vert_\infty\) is below \(10^{-12}\) while the relative energy-norm error exceeds \(5\%\). The criterion is satisfied when the calibrated indicator rejects every such query, that is, \(\eta_{\text{adm}} > \tau_{\text{adm}}\).""")
d.rep(r"the minimum over all (basis, \(r_\lambda\)) entries is the strict reading and is reported together with the median;",
      r"the minimum correlation across all basis--dimension combinations is reported as the conservative measure, together with the median;")

# ================================================================ Section 5: results wording
d.rep(r"Signed SVD is the cheapest to build, but its lifted solution is wrong: energy errors near one and adhesive multipliers.",
      r"Although the signed SVD basis has the lowest offline construction cost, it produces substantial errors in the reconstructed displacement field and negative contact multipliers, with energy-norm errors of order one.")
d.rep(r"gIS is cheap and non-negative but the least accurate of the non-negative bases (energy error 0.07 on B1 and 0.24 on B2) and needs six repair vectors.",
      r"The gIS strategy provides an inexpensive non-negative basis but has larger displacement errors than NMF and CPG at the tested dimension (energy error 0.07 on B1 and 0.24 on B2) and needs six repair vectors.")
d.rep(r"RID has the smallest energy error, but it retains additional contact rows (23 and 98 dual rows for a nominal 16), needs 3 and 78 repair vectors, has a large multiplier error on B1 (0.87), and is the slowest online on B2 (1.15\,ms). The values depend on \(r_\lambda\); only \(r_\lambda=16\) is shown.",
      r"RID achieves the smallest displacement errors in several configurations; this is accompanied by a larger effective dual dimension (23 and 98 dual rows for a nominal 16), 3 and 78 repair vectors, a large multiplier error on B1 (0.87) and the slowest online solve on B2 (1.15\,ms). The results vary with the dual dimension, and only \(r_\lambda=16\) is shown.")
d.rep(r"a cheap interface screen that is empirical", r"an inexpensive interface-based screening indicator that is empirical")
d.rep(r"these are an empirical screen and carry no error guarantee", r"these form an empirical screening indicator and carry no error guarantee")
d.rep(r"The screen uses the hybrid threshold", r"The screening indicator uses the hybrid threshold", count=2)
d.rep(r"silent failures missed by the screen /", r"silent failures missed by the screening indicator /")
d.rep(r"with the interface screen of Eq.~\eqref{eq:eta-adm}", r"with the interface-based screening indicator of Eq.~\eqref{eq:eta-adm}")
d.rep(r"Every basis degrades outside the training box, and the screen does not always follow.", r"Every basis degrades outside the training box, and the screening indicator does not always detect the failures.")
d.rep(r"The screen flagged almost all of these failures", r"The screening indicator flagged almost all of these failures")
d.rep(r"and the screen missed $78$ of the $90$ silent failures", r"and the screening indicator missed $78$ of the $90$ silent failures")
d.rep(r"This is the failure mode that Section~\ref{sec:conclusion} warns about, now shown on a benchmark: a screen calibrated only in range can be blind outside it.",
      r"This illustrates the limitation discussed in Section~\ref{sec:conclusion}: a screening threshold calibrated only in range may fail to detect failures outside it.")
d.rep(r"so on that benchmark the screen buys recall at the price of many unnecessary fallbacks.",
      r"so on that benchmark the increased detection rate is accompanied by a substantial number of unnecessary full-order fallback solves.")
d.rep(r"the screen has no empirical basis for that basis and dual dimension;", r"the screening indicator has no empirical basis for that basis and dual dimension;")
d.rep(r"The cost of the screen can exceed its benefit on the membrane:", r"The cost of screening can exceed its benefit on the membrane:")
d.rep(r"On this benchmark the screened ROM pays off only with varying stiffness,", r"On this benchmark the screened ROM provides a net speedup only with varying stiffness,")

# ================================================================ Section 6: summary
d.span(r"Projection-based reduced-order models for unilateral contact problems face a fundamental reliability hazard:", r"we addressed this gap through three interconnected contributions:",
       r"""Projection-based reduced-order models for unilateral contact problems may exhibit a significant discrepancy between convergence of the reduced system and the physical admissibility of the reconstructed solution. In particular, the reduced solver may satisfy its projected KKT conditions to machine precision (residuals below \(3 \times 10^{-14}\) in our experiments) while the lifted displacement--multiplier pair \((\tilde{u}, \tilde{\lambda}) = (V u_r, U \lambda_r)\) violates the full-order contact constraints, with relative energy-norm errors of order one \cite{balajewicz2016projection, benaceur2020reduced}.

This work makes three contributions.""")
d.rep(r"evaluates in \(\mathcal{O}(n_c (r_u + r_\lambda))\) operations without touching the bulk mesh, acting as an online circuit breaker against silent reduced-KKT convergence.",
      r"can be evaluated in \(\mathcal{O}(n_c (r_u + r_\lambda))\) operations without full-order bulk quantities and serves as an online admissibility screening indicator.")
d.rep(r"extrapolation robustness was tested only for NMF and CPG on the three-dimensional problem (experiment E5).",
      r"out-of-range robustness was tested for all five bases on B1 and B2 and for NMF and CPG on the three-dimensional problem (experiment E5).")
d.rep(r"and predicted speedups overshoot measured ones at large sizes by up to an order of magnitude.", r"and predicted speedups exceed measured ones at large sizes by up to an order of magnitude.")
d.rep(r"the interface screen detects most but not all silent failures at the price of false alarms that are high on some benchmarks,",
      r"the interface-based screening indicator detects most but not all silent failures at the price of false alarms that are high on some benchmarks,")

# ================================================================ Section 6: discussion
d.rep(r"\emph{When does each dual basis work well?} Signed SVD fails because nothing prevents negative lifted multipliers, so $\eta_{\mathrm{sign}}$ is large and the bound with it.",
      r"\paragraph{Performance of the dual-basis strategies.} The signed SVD basis does not preserve multiplier non-negativity, which allows negative components in the lifted multiplier vector and a large sign indicator $\eta_{\mathrm{sign}}$; the bound increases accordingly.")
d.rep(r"which is why their errors are of the order of a few percent or less, but their columns do not enforce complementarity at individual contact points, and on the membrane the error saturates at about $2$--$3\%$ as the dimension grows.",
      r"which is consistent with the comparatively smaller errors observed in the experiments. Their columns do not enforce complementarity at individual contact points, and on the membrane the error exhibits a plateau of approximately $2$--$3\%$ as the dimension grows.")
d.rep(r"so its errors are large until the selected set captures the points that are in contact;", r"so its errors are large until the selected contact indices adequately represent the active contact region;")
d.rep(r"\emph{When is the framework computationally beneficial?} Our measured", r"\paragraph{Computational break-even behaviour.} Our measured")
d.rep(r"With constant stiffness the saving only appears at larger sizes because the full-order cost is then the active-set solve alone.",
      r"With constant stiffness, a computational advantage is observed only beyond a sufficiently large problem size, because the full-order cost then consists of the active-set solve alone.")
d.rep(r"\emph{Limitations.} The interface screen can miss inadmissible solutions,", r"\paragraph{Limitations and applicability.} The interface-based screening indicator can miss inadmissible solutions,")

# ================================================================ Section 6: guidelines
d.rep(r"we summarize four actionable guidelines for practitioners deploying Lagrange-multiplier reduced contact models:", r"we summarize five guidelines for practitioners deploying Lagrange-multiplier reduced contact models:")
d.span(r"Always evaluate the \(\mathcal{O}(n_c (r_u + r_\lambda))\) online admissibility indicator", r"\cite{dabaghi2020posteriori}.",
       r"""Evaluate the online admissibility indicator \(\eta_{\text{adm}}\) after lifting \((u_r, \lambda_r)\) as an additional screening measure. Because the indicator is empirical and may miss inadmissible predictions, it should not be treated as a standalone guarantee of reliability. When it exceeds the calibrated screening threshold, reject the reduced prediction and invoke the full-order fallback \cite{dabaghi2020posteriori}. If such rejections occur frequently, the reduced basis may require additional enrichment \cite{leberre2022condition} or reconstruction.""")
OLD = chr(92) + "item" + chr(10) + "  " + chr(92) + "textbf{Use " + chr(92) + "(" + chr(92) + "eta_{" + chr(92) + "text{eq}}" + chr(92) + ") offline to enrich bulk primal modes during greedy training:}"
NEW5 = chr(92) + "item" + chr(10) + "  " + chr(92) + "textbf{Calibrate the screening threshold on queries that contain failures and treat it as valid only within the calibrated range:}" + chr(10) + "  The threshold is an empirical quantity. It cannot be calibrated when the validation set contains no silent failure, its recall depends on the number of validation queries, and in our out-of-range tests it missed most silent failures of one basis whose validation set had none. Report the recall and false-alarm rate together with the number of validation queries, and re-calibrate when the parameter range changes." + chr(10)
d.rep(OLD, NEW5 + OLD)

# ================================================================ Section 6: limitations
d.span(r"The experiments have a limited scope.", r"before any claim of certification could be made.",
       r"""The conclusions should be interpreted within the numerical scope of the study. The main tables use a single random seed (42) with 40 test and 40 validation queries per benchmark, and the five-seed check of Table~\ref{tab:robust} shows that the outcome of the detection criterion depends on the seed. The calibration study of Section~\ref{results-and-acceptance-outcome} (eight seeds, bootstrap intervals over seeds) quantifies this uncertainty on B1 and B2 only, and B3 was not repeated over seeds. Threshold calibration requires validation queries that contain silent failures, and several non-negative basis entries had none. The threshold was calibrated in range, and outside the training range its performance deteriorated for some bases (Table~\ref{tab:extrap}). An in-sample validation set was not evaluated.

The benchmarks and the cost model have further limitations. Three-dimensional results come from one mesh family of a single two-body geometry, and no second, independent three-dimensional configuration was run. Out-of-range queries were tested for one range of loads, obstacle shapes and indenter positions on B1 and B2, and for NMF and CPG only on B3. The dual-basis comparison at matched cost rests on one seed per benchmark, and the mesh-size study used a fixed dual dimension. The break-even model is validated on the membrane and on the three-dimensional problem; on B1 only the measured speedup was swept. The B3 varying-stiffness crossover is weakly determined at small sizes, the ten-repeat timing study shows run-to-run differences of up to about 40\% in the measured break-even and a varying-stiffness ratio on the membrane of 1.97, and the model overestimates absolute speedups at large sizes by up to an order of magnitude and assumes direct full-order solvers (Section~\ref{cost-model-assumptions}).

The effectivity and correlation results use the bound with the exact multiplier norm of the full-order solution, which is a validation quantity. The online-usable version with the offline constant \(C_\lambda\) is rigorous only while the query multiplier norm stays below \(C_\lambda\). It fell below the true error on up to \(0.5\%\) of in-range test queries for gIS and RID on B1 and on up to \(1.5\%\) (B1) and \(3.5\%\) (B2) of the out-of-range queries.

Two directions warrant further investigation. First, the experiments focus on frictionless unilateral contact; extending the constant-sign multiplier decomposition to Coulomb frictional contact \cite{niakh2024reduced} requires coupling normal traction admissibility with tangential friction-cone projection residuals. Second, integrating the online admissibility screening \(\eta_{\text{adm}}\) into an automated, on-the-fly basis enrichment loop for large-scale 3D nonlinear elastoplastic and substructured dynamic contact \cite{koller2021hyper, aresdeparga2023petrov, manvelyan2026inference} is a possible next step; it would need the bound with a valid multiplier-norm estimate and a calibrated screening threshold before any claim of certification could be made. The time-dependent extension was considered but not performed and is also left to future work.""")
d.save()
print("ok")
