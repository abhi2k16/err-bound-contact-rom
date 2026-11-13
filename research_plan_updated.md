# Research Plan: Reliability and Cost of Lagrange-Multiplier Reduced Contact Models

**Status:** UPDATED DRAFT (incorporating literature verification, estimator proof outline, FLOP-aware cost model, and 3D benchmark specification).  
**Codebase Foundation:** Built on the measured results of Chapter 14 (`packages_book/book_v3`) and existing Python assets (`ch13_assets.py`, `ch14_intro.py`, `ch14_basis.py`, `ch14_more.py`).

---

## 1. Working Title
**A Computable Admissibility Estimator and Break-Even Cost Model for Lagrange-Multiplier Reduced Contact Models**

---

## 2. Problem, Related Work, and Research Gap

### 2.1 Core Failure Mode: Silent Reduced-KKT Convergence
A reduced contact model with Lagrange multipliers can satisfy its own reduced Karush–Kuhn–Tucker (KKT) conditions to machine precision ($10^{-16}$) and still return a lifted primal–dual pair $(\tilde{u}, \tilde{\lambda}) = (V u_r, U \lambda_r)$ that is wrong by orders of magnitude. As measured in Chapter 14 of the book:
* Using a signed SVD dual basis yields lifted displacement errors of $3 \times 10^2$ to $4 \times 10^2$ and non-physical negative lifted multipliers down to $-1.4 \times 10^2$.
* Using a non-negative matrix factorization (NMF) dual basis at identical reduced dimensions yields displacement errors of $3 \times 10^{-2}$ to $9 \times 10^{-2}$ (an improvement of four orders of magnitude).

### 2.2 Verified Literature Positioning
1. **Non-Negative Dual Bases & Hyper-Reduction:**
   * Positivity-preserving dual constructions include NMF (Balajewicz et al., 2016), cone-projected greedy algorithms (Benaceur et al., 2020), greedy index selection (gIS) (Sambataro & Ehrlacher, 2026), and restriction to a reduced integration domain (RID) (Fauque et al., 2018; Le Berre et al., 2022), alongside empirical cubature and parametric contact reduction strategies (Koller et al., 2021; Ares de Parga et al., 2023; Pawar & Kulkarni, 2025; Manvelyan-Stroot et al., 2026).
   * *Gap:* None of these works provide a computable *a posteriori* error bound for the lifted primal–dual solution $(\tilde{u}, \tilde{\lambda})$, nor do they evaluate all five dual-basis strategies head-to-head on identical snapshot sets.
2. ***A Posteriori* Error Estimation for Variational Inequalities:**
   * Full-order residual-based and equilibrated *a posteriori* error estimators for unilateral contact decompose the error into equilibrium, sign-violation, and complementarity residuals without requiring Galerkin orthogonality (Coorevits et al., 2000; Hild & Nicaise, 2007; Dabaghi et al., 2020; Banz & Bertrand, 2024; Fontana & Di Pietro, 2024; Banz et al., 2026).
   * Certified duality-based reduced-basis error bounds (Zhang et al., 2014) assume a **bijective** (square, invertible) constraint operator $B$ and enforce dual feasibility by construction.
   * *Gap:* Existing ROM bounds do not cover general rectangular constraint matrices $B \in \mathbb{R}^{n_c \times N}$ arising in multi-body or non-matching contact, nor do they address cases where a small reduced KKT residual masks severe lifted inadmissibility.
3. **Conditioning & Computational Cost:**
   * Full-order saddle-point preconditioning (Franceschini et al., 2022) and projected contact conditioning (Le Berre et al., 2022) highlight the impact of dual basis precision.
   * *Gap:* The literature lacks an explicit break-even cost model quantifying when primal rank-repair enrichment, contact-row retention, active-set iteration growth, and dense reduced linear algebra eliminate the speedup over a sparse full-order solver.

---

## 3. Research Questions

* **RQ1 (Reliability):** Can a computable *a posteriori* estimator evaluated from the lifted pair $(\tilde{u}, \tilde{\lambda})$ reliably detect inadmissible reduced solutions that the reduced KKT residual accepts?
* **RQ2 (Basis Trade-Offs):** Which dual basis (NMF, cone-projected greedy, greedy index selection, signed SVD with rank repair, or RID restriction) achieves the best accuracy-to-admissibility trade-off across in-sample and out-of-range active-set regimes?
* **RQ3 (Computational Break-Even):** Under what problem dimensions, contact-zone sizes, rank-repair enrichment counts, and active-set iteration counts does the total online cost of the enriched ROM fall below that of a sparse full-order contact solver?

---

## 4. Proposed Method

### 4.1 Four-Component Lifted Admissibility Estimator
Consider the full-order discrete unilateral contact KKT system:
$$K u + B^T \lambda = f, \qquad \lambda \ge 0, \quad B u - c \le 0, \quad \lambda^T (B u - c) = 0,$$
where $K \in \mathbb{R}^{N \times N}$ is symmetric positive definite, $B \in \mathbb{R}^{n_c \times N}$ is a general rectangular constraint matrix ($n_c \le N$), and $c \in \mathbb{R}^{n_c}$ is the initial gap vector.

Given any lifted reduced-order pair $(\tilde{u}, \tilde{\lambda}) = (V u_r, U \lambda_r)$, define the positive and negative multiplier parts $\lambda_+ := \max(\tilde{\lambda}, 0) \ge 0$ and $\lambda_- := \min(\tilde{\lambda}, 0) \le 0$, and the four estimator components:
* **Multiplier sign violation:** $\eta_{\text{sign}} := \|\lambda_-\|_2$, online cost $\mathcal{O}(n_c r_\lambda)$.
* **Primal gap penetration:** $\eta_{\text{gap}} := \|(B \tilde{u} - c)_+\|_2$, online cost $\mathcal{O}(n_c r_u)$ via precomputed $B V \in \mathbb{R}^{n_c \times r_u}$.
* **Complementarity defect:** $\eta_{\text{comp}} := |\lambda_+^T (B \tilde{u} - c)|$, online cost $\mathcal{O}(n_c)$.
* **Dual-norm equilibrium residual:** $\eta_{\text{eq}} := \|K \tilde{u} + B^T \lambda_+ - f\|_{K^{-1}}$, full-order validation / offline greedy selection metric (or $\mathcal{O}((r_u + r_\lambda)^2)$ online via offline Riesz representation decomposition in affine parameter settings).

**Certified Energy-Norm Upper Bound (Hild–Nicaise Residual Split):**
Without assuming bijectivity of $B$ or reduced Galerkin orthogonality, expanding $\|u - \tilde{u}\|_K^2$ and splitting the contact cross-term $(B(u - \tilde{u}))^T (\lambda_+ - \lambda)$ yields:
$$\|u - \tilde{u}\|_K^2 \le \eta_{\text{eq}}^2 + 2\,\eta_{\text{comp}} + 2\,\|\lambda\|_2\,\eta_{\text{gap}} \le 2\,\tilde{r}_{\text{eq}}^2 + 2\,\|B\|_{K \to 2}^2\,\eta_{\text{sign}}^2 + 2\,\eta_{\text{comp}} + 2\,\|\lambda\|_2\,\eta_{\text{gap}},$$
where $\tilde{r}_{\text{eq}} := \|K \tilde{u} + B^T \tilde{\lambda} - f\|_{K^{-1}}$ is the unprojected equilibrium residual and $\|B\|_{K \to 2}$ is the full-order discrete trace norm.

**Online Admissibility Surrogate ("Circuit Breaker"):**
To monitor physical admissibility strictly in $\mathcal{O}(n_c (r_u + r_\lambda))$ operations without touching the bulk mesh, define:
$$\eta_{\text{adm}} := \eta_{\text{sign}} + \eta_{\text{gap}} + \sqrt{\eta_{\text{comp}}},$$
where $\sqrt{\eta_{\text{comp}}}$ ensures all three terms scale linearly with $\|u - \tilde{u}\|_K$.

### 4.2 FLOP- and Iteration-Aware Break-Even Cost Model
Let $r_u$ denote POD displacement modes, $r_{\text{rep}}$ denote unit-vector rank-repair modes (e.g., $14$ on the book's discrete chain), $r_\lambda$ denote reduced multipliers, and $r_{\text{tot}} = r_u + r_{\text{rep}} + r_\lambda$.
* **Dense ROM Online FLOPs:**
  $$\mathcal{C}_{\text{ROM}} = k_{\text{ROM}} \cdot \frac{2}{3} r_{\text{tot}}^3 + \mathcal{O}\big(n_c (r_u + r_{\text{rep}} + r_\lambda)\big),$$
  explicitly accounting for active-set/LCP iteration count $k_{\text{ROM}}$ and dense saddle-point factorization.
* **Sparse FOM Online FLOPs:**
  $$\mathcal{C}_{\text{FOM}} = k_{\text{FOM}} \cdot \mathcal{O}(N \cdot \text{bw}^2) \quad (\text{or } \mathcal{O}(N^\alpha) \text{ for sparse Cholesky/AMG}).$$
* **Break-Even Condition:** Predicts the critical contact-zone fraction $n_c / N$, rank-repair overhead $r_{\text{rep}}$, and iteration ratio $k_{\text{ROM}} / k_{\text{FOM}}$ where $\mathcal{C}_{\text{ROM}} \ge \mathcal{C}_{\text{FOM}}$.

---

## 5. Experiments

* **E1 (Failure Reproduction & Component Audit):** Reproduce the silent reduced-KKT failure on the 324-DOF obstacle membrane and 60-spring chain; compute all four estimator components ($\eta_{\text{sign}}, \eta_{\text{gap}}, \eta_{\text{comp}}, \eta_{\text{eq}}$) across every dual basis and rank.
* **E2 (Effectivity & Parameter Extrapolation Stress Test):** Evaluate the effectivity index $\eta / \|u - \tilde{u}\|_K$ across three regimes:
  1. In-sample training parameters;
  2. *Load-magnitude extrapolation* (increasing applied force so the contact patch grows beyond the training snapshots);
  3. *Contact-patch shift* (translating the load so unseen rows of $B$ activate).
* **E3 (Unified Five-Basis Benchmark):** Compare NMF, cone-projected greedy, greedy index selection (gIS), signed SVD with rank repair, and RID restriction on identical snapshot sets across lifted displacement/multiplier error, $\eta_{\text{adm}}$ violations, active-set iterations $k_{\text{ROM}}$, online wall-clock time, and offline construction time.
* **E4 (Cost Model Validation):** Compare predicted vs. measured FLOP and wall-clock break-even curves on the 1D chain, 2D membrane, and 3D contact benchmark.
* **E5 (3D Two-Body Contact Scale-Up):** Simulate a 3D Hertzian elastic hemisphere-on-block (or bolted flange) frictionless contact problem with $N \in [2 \times 10^4, 10^5]$ DOFs. Use **FEniCSx + PETSc TAO/SNES** (or **GetFEM**) to assemble and export exact $(K, B, c, f)$ matrices directly into the existing Python `ContactROM` pipeline.
* **E6 (Optional Dynamic Extension):** Evaluate $\eta_{\text{adm}}$ per time step on the BDF2 reduced LCP dynamics solver (`solve_contact_dynamics`).

---

## 6. Success Criteria (Set Before Running)

* **Zero False Negatives on Silent Contact Failures:** $100\%$ of cases where the reduced KKT residual is below $10^{-10}$ while relative lifted error exceeds $5\%$ are flagged by $\eta_{\text{adm}}$.
* **Estimator Effectivity:** Full four-component estimator effectivity index between $1$ and $10$ on at least $90\%$ of cases for the best-performing basis, with explicit characterization of bulk-dominated failure modes where $\eta_{\text{eq}}$ dominates.
* **Monotonicity / Rank Correlation:** Spearman rank correlation $> 0.90$ between the full estimator and the true energy-norm error across bases and reduced dimensions.
* **Break-Even Accuracy:** Predicted FLOP/wall-clock break-even threshold within a factor of $2$ of measured wall-clock crossover on the 2D and 3D benchmarks.

---

## 7. Risks and Mitigations

* **3D Solver Implementation Overhead:** Building a 3D mortar or node-to-segment solver from scratch risks schedule slippage. *Mitigation:* Use FEniCSx or GetFEM to export exact $(K, B, c, f)$ matrices into the existing Python ROM framework.
* **Offline Cost of $\eta_{\text{eq}}$:** Evaluating $K^{-1}$ requires a full-order solve in non-affine problems. *Mitigation:* Restrict $\eta_{\text{eq}}$ to offline validation and greedy basis enrichment; use $\eta_{\text{adm}}$ online, and document the offline Riesz representation extension for affine parameter dependence.
* **Overlap with Le Berre et al. (2022) on Conditioning:** *Mitigation:* Cite explicitly and demonstrate the quantitative impact of rank-repair enrichment $r_{\text{rep}}$ and active-set iteration growth $k_{\text{ROM}}$ on computational break-even.
* **Recent Preprint Follow-Ups (Sambataro & Ehrlacher, 2026; Banz et al., 2026):** *Mitigation:* Monitor arXiv prior to submission and position gIS directly inside Benchmark E3.
* **Dynamic E6 Scope Creep:** *Mitigation:* Keep E6 strictly as a brief concluding demonstration or omit if E1–E5 fill the journal page budget.

---

## 8. Timeline (8 Months)

* **Months 1–2:** Finalize estimator theorem proof; implement $\eta_{\text{adm}}$ and run E1–E2 on existing `ch13_assets.py` benchmarks.
* **Months 3–4:** Implement unified 5-basis comparison (E3) and parameter-shift extrapolation tests.
* **Months 5–6:** Set up 3D FEniCSx/GetFEM matrix export (E5) and validate break-even cost model (E4).
* **Month 7:** Run E6 (if time permits), generate publication figures/tables, and draft manuscript.
* **Month 8:** Internal review, open-source repository preparation, and submission.

---

## 9. Target Journals and Reproducibility Plan

* **Primary Venues:** *Computer Methods in Applied Mechanics and Engineering* (*CMAME*) or *International Journal for Numerical Methods in Engineering* (*IJNME*).
* **Secondary Venue:** *Computers & Structures*.
* **Data & Code Availability:** Release benchmark matrices $(K, B, c, f)$, unified Python implementations of all five dual-basis builders, and estimator verification scripts in a public repository upon submission.

---

## 10. Reuse from Existing Package
* `ch13_assets.py`: `ObstacleMembrane`, `NonlinearChainContact`, `solve_contact_dynamics`, `ContactROM`, `nonneg_dual_basis`, `svd_dual_basis`, `greedy_contact_sampling`, `nonlinearity_indicator`.
* `ch14_intro.py`, `ch14_basis.py`, `ch14_more.py`: Measured Chapter 14 failure benchmarks and baseline plots.

---

## 11. Literature Verification Log (Resolved)
* **Balajewicz et al. (2016) & Benaceur et al. (2020):** Confirmed absent of lifted multiplier *a posteriori* error bounds.
* **Banz et al. (2026, arXiv:2609.26001) & Banz & Bertrand (2024, *CAMWA*):** Confirmed abstract and two-field full-order contact error frameworks without Galerkin orthogonality; integrated into Section 2 and Section 4.
* **Sambataro & Ehrlacher (2026, arXiv:2605.04037):** Confirmed gIS dual basis construction; integrated into E3 benchmark.
* **Elsevier Precedents Added:** Coorevits et al. (2000), Dabaghi et al. (2020), Koller et al. (2021), Franceschini et al. (2022), Ares de Parga et al. (2023), Fontana & Di Pietro (2024), Banz & Bertrand (2024), and Pawar & Kulkarni (2025).

---

## 12. References

1. Ares de Parga, S., Hernández, J. A., & Rossi, R. (2023). Petrov–Galerkin hyper-reduction via the empirical cubature method for non-symmetric Jacobians. *Computer Methods in Applied Mechanics and Engineering*, 414, 116180.
2. Balajewicz, M., Amsallem, D., & Farhat, C. (2016). Projection-based model reduction for contact problems. *International Journal for Numerical Methods in Engineering*, 106(8), 644–663.
3. Banz, L., & Bertrand, F. (2024). Generic a posteriori error estimators for two-field contact problems. *Computers & Mathematics with Applications*, 169, 188–205.
4. Banz, L., Schönauer, M., & Schröder, A. (2026). A posteriori error estimates for variational inequalities of the second kind in an abstract framework. *arXiv preprint arXiv:2609.26001*.
5. Benaceur, A., Ern, A., & Ehrlacher, V. (2020). A reduced basis method for parametrized variational inequalities applied to contact mechanics. *International Journal for Numerical Methods in Engineering*, 121(6), 1170–1197.
6. Coorevits, P., Hild, P., & Pelle, J.-P. (2000). A posteriori error estimation for unilateral contact with matching and non-matching meshes. *Computer Methods in Applied Mechanics and Engineering*, 186(1), 65–83.
7. Dabaghi, J., Martin, V., & Vohralík, M. (2020). A posteriori error estimates distinguishing the error components and adaptive stopping criteria for numerical approximations of parabolic variational inequalities. *Computer Methods in Applied Mechanics and Engineering*, 367, 113105.
8. Fauque, J., Ramière, I., & Ryckelynck, D. (2018). Hybrid hyper-reduced modeling for contact mechanics problems. *International Journal for Numerical Methods in Engineering*, 115(1), 117–139.
9. Fontana, M., & Di Pietro, D. A. (2024). Guaranteed a posteriori error bounds for finite element approximations of contact problems. *Computer Methods in Applied Mechanics and Engineering*, 424, 116895.
10. Franceschini, A., Castelletto, N., White, J. A., & Tchelepi, H. A. (2022). Algebraically-preconditioned Krylov solvers for Lagrange-multiplier contact problems. *Computer Methods in Applied Mechanics and Engineering*, 392, 114632.
11. Hild, P., & Nicaise, S. (2007). Residual a posteriori error estimators for contact problems in elasticity. *ESAIM: Mathematical Modelling and Numerical Analysis*, 41(5), 897–923.
12. Kollepara, K. S., Aguado, J. V., Le Guennec, Y., Silva, L., & Borzacchiello, D. (2024). Sparse approximations for contact mechanics. *arXiv preprint arXiv:2406.11461*.
13. Koller, L., Rixen, D. J., et al. (2021). Hyper-reduction for finite element structures with nonlinear surface loads and contact. *Computer Methods in Applied Mechanics and Engineering*, 385, 114023.
14. Le Berre, S., Ramière, I., Fauque, J., & Ryckelynck, D. (2022). Condition number and clustering-based efficiency improvement of reduced-order solvers for contact problems using Lagrange multipliers. *Mathematics*, 10(9), 1495.
15. Manvelyan-Stroot, D., Filanova, Y., Pontes Duff, I., Benner, P., & Wever, U. (2026). Inference of substructured reduced-order models for dynamic contact from contact-free simulations. *Computational Mechanics*.
16. Niakh, I., Drouet, G., Ehrlacher, V., & Ern, A. (2024). A reduced basis method for frictional contact problems formulated with Nitsche's method. *SMAI Journal of Computational Mathematics*, 10, 29–54.
17. Pawar, G. S., & Kulkarni, S. S. (2025). Efficient parametric model order reduction in contact mechanics. *Computers & Structures*, 315, 107838.
18. Sambataro, G., & Ehrlacher, V. (2026). Model order reduction for parametrized variational inequalities: application to crowd motion. *arXiv preprint arXiv:2605.04037*.
19. Zhang, Z., Bader, E., & Veroy, K. (2014). A duality approach to error estimation for variational inequalities. *arXiv preprint arXiv:1410.2095*.
