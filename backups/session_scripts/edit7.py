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


def insert_before(anchor, text):
    global s
    assert s.count(anchor) == 1, (s.count(anchor), anchor[:80])
    i = s.index(anchor)
    s = s[:i] + text + s[i:]


def insert_after(anchor, text):
    global s
    assert s.count(anchor) == 1, (s.count(anchor), anchor[:80])
    i = s.index(anchor) + len(anchor)
    s = s[:i] + text + s[i:]


# ---- "certified" wording -------------------------------------------------------------------------------------
rep(r"\subsection{Certified Energy-Norm Upper Bound}", r"\subsection{Energy-Norm Upper Bound}")
rep(r"To systematically evaluate the certified \emph{a posteriori} estimator (Section 3)", r"To systematically evaluate the \emph{a posteriori} estimator (Section 3)")
rep(r"We evaluate the certified upper bound from Theorem", r"We evaluate the upper bound from Theorem")
rep(r"For the certified bound \(\eta_{\text{bound}}\)", r"For the bound \(\eta_{\text{bound}}\)")
rep(r"represents a natural next step toward certified industrial contact digital twins.",
    r"is a possible next step; it would need the bound with a valid multiplier-norm estimate and a calibrated screening threshold before any claim of certification could be made.")

# ---- Section 2.1: feasibility assumptions --------------------------------------------------------------------
insert_before(r"Because \(K(\mu)\) is SPD, a unique primal solution",
              r"""We assume that the admissible set is non-empty for every \(\mu \in \mathcal{P}\), that is, \(\{v : B(\mu) v \le c(\mu)\} \neq \emptyset\); this holds, for example, whenever \(c(\mu) \ge 0\), since then \(v = 0\) is admissible, which is the case in all our benchmarks. The set \(\mathcal{K}(\mu)\) is then a non-empty closed convex polyhedron. Because the constraints are linear, a KKT multiplier \(\lambda(\mu)\) exists without further constraint qualifications. """)

# ---- Section 3.1: sign convention ------------------------------------------------------------------------------
insert_before(r"Using \(\lambda_+\) as a feasible dual proxy",
              r"""Sign convention: the constraint is \(B u - c \le 0\) and \(\lambda \ge 0\), so a positive multiplier component is a compressive contact traction and \(\lambda_+\) is its admissible part, while \(\lambda_- \le 0\) collects the non-physical adhesive components (the reduction packages used in our code use the opposite sign for both \(B\) and \(\lambda\); we convert before evaluating any indicator). Accordingly \(\eta_{\mathrm{sign}} = \|\lambda_-\|_2\) is the magnitude of the adhesive part, \((B\tilde{u} - c)_+\) collects penetrations, and \(\lambda_+\), \(\lambda_-\) have disjoint supports. """)

# ---- Section 3.3: parameter-dependent B ------------------------------------------------------------------------
insert_after(r"---completely independent of the bulk mesh dimension \(N\).",
             r""" When \(B(\mu)\) or \(c(\mu)\) is not affine, the product \(B(\mu) V\) cannot be precomputed. It is then formed per query as \(B_{:,\Gamma}(\mu) V_{\Gamma,:}\), where \(\Gamma\) is the set of degrees of freedom touched by the contact rows, at the cost \(\mathcal{O}(\mathrm{nnz}(B)\, r_u)\) instead of \(\mathcal{O}(n_c r_u)\) for the product with a precomputed \(BV\); for node-to-segment type constraints with a bounded number of degrees of freedom per row this is still of order \(n_c r_u\) and independent of \(N\), but the constant grows with the stencil. Only the rows of \(V\) on \(\Gamma\) are needed online. In our benchmarks \(B\) and \(c\) do not depend on \(\mu\) except through \(c\), which enters as a vector.""")

# ---- Algorithm 4: failure case -------------------------------------------------------------------------------
insert_after(r"because the columns of \(U\) are orthonormal.",
             r""" The repair loop stops when \(\kappa \le \kappa_{\max}\), after at most \(r_\lambda\) added vectors, or when the supremizer already lies in the span of \(V\); in the last two cases the target is not reached and the algorithm returns the enriched basis with the condition number it achieved. The implementation records this value but does not abort. In the baseline E3 audit on B1 three entries ended above \(\kappa_{\max} = 10^{6}\) (NMF with \(r_\lambda = 24\): \(4.4\times10^{7}\); CPG with \(r_\lambda = 16\): \(6.8\times10^{7}\); CPG with \(r_\lambda = 24\): infinite, a numerically singular reduced Schur complement), and none did on B2. A reduced model with an ill-conditioned or singular Schur complement can lose accuracy in the reduced LCP solve, so such entries should not be deployed as they are: either \(r_\lambda\) should be reduced, discarding the weakest dual modes, or the primal space enlarged by other means, or the query should be sent to the full-order solver. The reported accuracy of these entries is therefore that of a model that did not meet its own conditioning target.""")

rep(r"required to maintain \(\text{cond}(B_r K_r^{-1} B_r^T) \le 10^6\) \cite{leberre2022condition};",
    r"required to bring \(\text{cond}(B_r K_r^{-1} B_r^T)\) to at most \(10^6\) (reached in all cases except three B1 entries, see the discussion of Algorithm~\ref{alg:repair}) \cite{leberre2022condition};")

# ---- Acceptance table legend ------------------------------------------------------------------------------------
i = s.index(r"4. Break-even within a factor 2")
j = s.index(r"\end{longtable}", i) + len(r"\end{longtable}")
s = s[:j] + r"""

\emph{How to read the table.} Entries of criterion 1 are missed silent failures over silent failures, pooled over the five dual dimensions and scored on the test queries with the threshold calibrated on the 40 validation queries (seed 42). ``No threshold exists'' means that, for that basis and dual dimension, the validation set contained no silent failure, so the rule $\tau_{\mathrm{adm}} = 0.9\min(\eta_{\mathrm{adm}})$ is undefined; such entries are left out of the counts. A count of 0 in the table is therefore a statement about the calibrated entries only, and the denominators differ between bases. ``Not evaluable'' means the test set itself contained no silent failure for that basis. Section~\ref{results-and-acceptance-outcome} gives the threshold calibration study that quantifies these effects with more seeds and validation queries.""" + s[j:]

# ---- Reproducibility paragraph ---------------------------------------------------------------------------------------
anchor = r"\subsection{Benchmark Test Cases}"
k = s.index(anchor)
prev = s.rfind("\n", 0, k - 1) + 1
insert_at = prev if s[prev:k].startswith(r"\hypertarget") else k
repro = r"""\paragraph{Implementation and reproducibility.}
All experiments use the Python package \texttt{contact\_rom} that accompanies the manuscript, on Python 3.14.3 with NumPy 2.5.4 and SciPy 1.18.1, the finite-element package \texttt{fea\_engine} 1.0.1 and the reduction package \texttt{rom\_engine} 0.2.0; every run logs these versions with its configuration. Timings were taken on a Windows 11 laptop with an Intel Core i7-9750H (6 cores, 12 threads) and 16\,GB of memory, with BLAS and OpenMP restricted to one thread. Full-order sparse systems are factored with SuperLU and the contact problems are solved by a dense active-set LCP method for both the full-order and the reduced model. Tolerances: full-order KKT verification $10^{-10}$, matrix symmetry $10^{-12}$, reduced KKT residual counted as converged below $10^{-12}$, relative energy error counted as failure above $5\%$, rank-repair target $\kappa_{\max} = 10^{6}$, primal dimension $r_u = 20$ before repair. Parameters are drawn independently and uniformly in the boxes of the benchmark descriptions below from one seeded NumPy generator (default seed 42; training 60, test 40, validation 40 unless stated otherwise), training first, then test, then validation, so the training set does not depend on the sizes of the later sets. B1 is a chain of 60 nodes with spring stiffness 1 and a weak ground spring $0.05$, with $\mu_1 \in [0.05, 0.3]$ and curvature $\mu_2 \in [1, 6]$; B2 is an $18\times18$ interior grid with a $12\times12$ contact sub-grid and $F_0 \in [0.05, 0.25]$, centre $(x_c, y_c) \in [0.35, 0.65]^2$; B3 is a quarter-symmetric two-body problem of two blocks of Hex8 elements ($20\times20$ elements per side and 6 layers per body in the main run, Poisson ratio $0.3$) with stiffness ratio $\in [0.5, 2]$, peak pressure $\in [0.01, 0.06]$ and initial clearance $\in [0.05, 0.3]$. The commands that regenerate every table and figure are listed in the repository \texttt{README}.

"""
s = s[:insert_at] + repro + s[insert_at:]

# ---- Introduction: positioning ---------------------------------------------------------------------------------------
fig = s.index(r"\begin{figure}[htbp]")
assert s.rfind(r"\end{itemize}", 0, fig) > s.rfind(r"\begin{itemize}", 0, fig)
pos = r"""\paragraph{Positioning of the contributions.}
The three contributions are of different kinds and should be weighed separately. (i)~\emph{Extension of a known decomposition.} The splitting of the error into an equilibrium residual and sign, gap and complementarity defects is established for full-order contact \cite{hild2007residual, dabaghi2020posteriori, banz2024generic, banz2026posteriori}; we apply it to lifted reduced pairs that are neither Galerkin-orthogonal nor dually feasible, and we show that the sign defect is exactly what separates signed from non-negative dual bases. (ii)~\emph{A bound.} Theorem~\ref{thm:bound} follows from the energy identity and the KKT conditions with elementary inequalities, and we do not claim a new proof technique. Its distinguishing property relative to the duality-based bound for reduced variational inequalities \cite{zhang2014duality} is that it needs neither a bijective constraint operator nor a feasible dual approximation by construction, and relative to reduced-basis variational inequalities \cite{haasdonk2012reduced} that it applies to rectangular $B$; its practical value is limited by the multiplier norm it contains (Remark~1) and by the equilibrium residual, which is not cheap in general. (iii)~\emph{A screening mechanism.} The interface indicators are an empirical detector with a calibrated threshold, not a certificate. Existing contact reduction strategies, namely positive dual bases \cite{balajewicz2016projection, benaceur2020reduced, sambataro2026model}, reduced integration domains \cite{fauque2018hybrid, leberre2022condition}, empirical cubature and operator-inference approaches \cite{koller2021hyper, aresdeparga2023petrov, manvelyan2026inference}, aim at admissible reduced solutions but, to our knowledge, provide no per-query error check; the authors of the greedy index selection explicitly postpone an error indicator to future work \cite{sambataro2026model}. Our estimator is complementary to these constructions and can be applied to any of them.

"""
s = s[:fig] + pos + s[fig:]

# ---- Limitations -------------------------------------------------------------------------------------------------
rep(r"rates such as false alarms or misses carry sampling uncertainty that we did not quantify with confidence intervals.",
    r"the calibration study of Section~\ref{results-and-acceptance-outcome} (eight seeds, bootstrap intervals over seeds) quantifies this uncertainty on B1 and B2 only, and B3 was not repeated over seeds.")
rep(r"repeated timing runs differ by a factor of about 1.4 in the measured break-even, and the model overestimates absolute speedups at large sizes.",
    r"the ten-repeat timing study shows run-to-run differences of up to about 40\% in the measured break-even and a varying-stiffness ratio on the membrane of 1.97, and the model overestimates absolute speedups at large sizes by up to an order of magnitude and assumes direct full-order solvers (Section~\ref{cost-model-assumptions}).")
rep(r"Three-dimensional results come from one mesh family of a single two-body geometry.",
    r"Three-dimensional results come from one mesh family of a single two-body geometry, and no second, independent three-dimensional configuration was run; the dual-basis comparison at matched cost rests on one seed per benchmark.")

# ---- Conclusion: conservative statement ----------------------------------------------------------------------------
insert_before(r"\hypertarget{practical-guidelines-for-contact-rom-deployment}{%",
              r"""The evidence supports the framework as a diagnostic and planning tool inside the tested range: the error bound holds with the exact multiplier norm on every audited query, signed dual bases fail silently, non-negative bases reduce but do not remove the hazard, and the break-even sizes of the cost model are located within roughly a factor of two for direct full-order solvers. It does not support a general claim of certified online reliability. The online-usable bound depends on an offline multiplier-norm constant, the interface screen detects most but not all silent failures at the price of false alarms that are high on some benchmarks, its threshold needs validation queries that contain failures, and the cost model predicts break-even sizes, not speedups. Results on other geometries, contact laws, solvers and hardware may differ.

""")

io.open(p, "w", encoding="utf-8", newline="").write(s.replace(LF, CR + LF) if crlf else s)
print("ok")
