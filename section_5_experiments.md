# Section 5: Experiments — All Versions

---

## v2 (current)

### 5. Numerical Experiments Protocol

To systematically evaluate the certified *a posteriori* estimator (Section 3), compare the five dual-basis constructions under controlled conditions (Section 4.1), and validate the FLOP- and iteration-aware break-even cost model (Section 4.2), we design a six-part experimental benchmark suite spanning 1D, 2D, and 3D unilateral contact problems. All reduced-order models are evaluated on identical snapshot sets $\{(u(\mu_i), \lambda(\mu_i))\}_{i=1}^{n_s}$ and solved using a unified active-set / Lemke LCP solver converged to a reduced KKT residual tolerance of $10^{-14}$ \cite{balajewicz2016projection, benaceur2020reduced, leberre2022condition}.

#### 5.1 Benchmark Test Cases
We consider three progressively larger frictionless unilateral contact problems covering distinct interface topologies and sparsity structures:

1. **Benchmark B1: 1D Discrete Spring-Chain Contact ($N = 60$, $n_c = 30$)**
   A chain of $N = 60$ linear springs with a weak ground spring, under a uniform load and constrained at every second node by a rigid parabolic obstacle $c(\mu)$ \cite{pawar2025efficient, balajewicz2016projection}. The parameter vector $\mu = (\mu_1, \mu_2) \in \mathcal{P} \subset \mathbb{R}^2$ controls the load amplitude and the obstacle curvature. Because $n_c / N = 0.50$ is large, this benchmark isolates rank-repair overhead \cite{leberre2022condition} in the constant-$K$ regime (the supremizer repair added between $0$ and $7$ vectors, depending on the basis and $r_\lambda$).
2. **Benchmark B2: 2D Parametrized Obstacle Membrane ($N = 324$, $n_c = 144$)**
   A clamped two-dimensional elastic membrane discretized on an $18 \times 18$ interior grid ($N = 324$) pressed against a rigid curved indenter \cite{benaceur2020reduced, hild2007residual, coorevits2000posteriori}. The parameter vector $\mu = (F_0, x_c, y_c) \in \mathbb{R}^3$ varies the peak load magnitude $F_0$ and the indenter center $(x_c, y_c)$, causing the active contact patch $m_{\text{act}}$ to expand, contract, and translate across the candidate interface.
3. **Benchmark B3: 3D Two-Body Contact ($N = 16\,653$, $n_c = 441$)**
   A quarter model of two stacked linear-elastic blocks, discretized with trilinear hexahedral elements on a uniform $20 \times 20$ grid with six element layers per block, with a matched node-to-node frictionless contact interface \cite{franceschini2022block, koller2021hyper, fauque2018hybrid} and a curved initial gap $\delta_0 + \kappa r^2$. A Gaussian pressure acts on the upper face, symmetry planes remove the other quadrants, the lower base is clamped, and a weak ground spring on the upper face keeps $K$ positive definite before contact. The element stiffness is taken from the Hex8 element of the `fea_engine` package and the matrices are assembled in SciPy; the full-order solve is a sparse LU factorization with an active-set LCP, cross-checked against the sparse primal--dual contact solver of the same package. Parameters $\mu = (E_1/E_2, P_{\text{load}}, \delta_0)$ vary the stiffness contrast (testing parameter-varying $K(\mu)$ \cite{koller2021hyper}), applied pressure, and initial clearance.

| Benchmark | Spatial Dim. | Full-Order $N$ | Contact Rows $n_c$ | Stiffness $K(\mu)$ | Primary Purpose |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **B1 (Spring Chain)** \cite{pawar2025efficient} | 1D | $60$ | $30$ | Constant | Rank-repair cost crossover & silent SVD failure |
| **B2 (Obstacle Membrane)** \cite{benaceur2020reduced} | 2D | $324$ | $144$ | Constant | Five-basis comparison & contact-patch translation |
| **B3 (3D Two-Body Contact)** \cite{franceschini2022block, koller2021hyper} | 3D | $16\,653$ | $441$ | Parametrized $K(\mu)$ | Large-scale wall-clock speedup & 3D effectivity |

#### 5.2 Experimental Campaign (E1--E6)

* **Experiment E1: Silent Reduced-KKT Failure Audit (B1, B2)**
  Using identical primal POD bases $V_0$ and signed-SVD and NMF dual bases $U$ of matched dimension $r_\lambda \in \{4, 8, 12, 16, 20\}$, we log the reduced KKT residual $\Vert R_{\text{red}}(u_r, \lambda_r)\Vert_\infty$ alongside the true relative energy-norm error $\Vert u - \tilde{u}\Vert_K / \Vert u\Vert_K$, the four estimator terms $(\eta_{\text{eq}}, \eta_{\text{sign}}, \eta_{\text{gap}}, \eta_{\text{comp}})$ \cite{hild2007residual, banz2024generic, banz2026posteriori}, and the online admissibility indicator $\eta_{\text{adm}}$. This experiment tests whether $\eta_{\text{adm}}$ flags every case where $\Vert R_{\text{red}}\Vert_\infty \le 10^{-12}$ masks catastrophic lifted error.

* **Experiment E2: Estimator Effectivity and Extrapolation Stress Test (B1, B2, B3)**
  We evaluate the certified upper bound from Theorem 1 and compute the **effectivity index** \cite{hild2007residual, zhang2014duality, fontana2024guaranteed}
  $$\theta_{\text{eff}}(\mu) := \frac{\eta_{\text{bound}}(\mu)}{\Vert u(\mu) - \tilde{u}(\mu)\Vert_K}, \qquad \eta_{\text{bound}}(\mu) := \left( \eta_{\text{eq}}^2 + 2\,\eta_{\text{comp}} + 2\,\Vert \lambda\Vert_2\,\eta_{\text{gap}} \right)^{1/2}, \qquad \eta_{\text{bound}}^{C}(\mu) := \left( \eta_{\text{eq}}^2 + 2\,\eta_{\text{comp}} + 2\,C_\lambda\,\eta_{\text{gap}} \right)^{1/2}.$$
  on two kinds of test sets:
  1. *In-range test queries* drawn from the training range of the parameters (B1, B2, B3);
  2. *Load-magnitude extrapolation* on B3 ($P$ between $1.3$ and $2.0$ times the training maximum), expanding the active contact set beyond the training snapshots \cite{benaceur2020reduced}.
  Contact-patch translation on the membrane and an in-sample validation set were not run.

* **Experiment E3: Head-to-Head Five-Basis Benchmark (B1, B2)**
  On identical snapshot matrices $(U_{\text{snap}}, \Lambda_{\text{snap}})$, we construct and compare all five dual-basis strategies from Section 4.1---**SVD+Repair** \cite{leberre2022condition}, **NMF** \cite{balajewicz2016projection}, **CPG** \cite{benaceur2020reduced}, **gIS** \cite{sambataro2026model, kollepara2024sparse}, and **RID** \cite{fauque2018hybrid, leberre2022condition}---across dual dimensions $r_\lambda \in \{4, 8, 12, 16, 24\}$. Each method is evaluated along five Pareto axes:
  1. Lifted energy-norm displacement error and multiplier $\ell_2$-error;
  2. Online admissibility violation $\eta_{\text{adm}}$;
  3. Primal rank-repair enrichment $r_{\text{rep}}$ required to maintain $\text{cond}(B_r K_r^{-1} B_r^T) \le 10^6$ \cite{leberre2022condition};
  4. Active-set iteration count $k_{\text{ROM}}$ and online wall-clock time;
  5. Offline basis construction wall-clock time.

* **Experiment E4: Break-Even Cost Model Validation (B1, B2, B3)**
  We measure single-thread wall-clock execution times as functions of the enriched primal dimension $r_u = r_{u,0} + r_{\text{rep}}$, active contact size $m_{\text{act}}$, and mesh size $N$ \cite{franceschini2022block, leberre2022condition, pawar2025efficient}. Measured speedup curves $S_{\text{meas}} = t_{\text{FOM}} / t_{\text{ROM}}$ are plotted against the theoretical predictions $S_{\text{pred}} = \mathcal{C}_{\text{FOM}} / \mathcal{C}_{\text{ROM}}$ from Equations (16)--(17) for the constant-$K$ and parameter-varying-$K$ cost models. The per-basis comparison covers B1 and B2, and the size sweeps cover the membrane ($N = 25$--$2304$) and B3 ($N = 24$--$16\,653$); B1 has a fixed size and is not swept. Break-even is the size from which the speedup stays above one.

* **Experiment E5: Large-Scale 3D Contact Certification (B3)**
  On the 3D two-body contact model ($N = 16\,653$, $n_c = 441$) \cite{franceschini2022block, koller2021hyper}, we deploy NMF and CPG dual bases ($r_\lambda \in \{8, 16\}$) with the online circuit breaker $\eta_{\text{adm}} \le \tau_{\text{adm}}$ \cite{dabaghi2020posteriori}, with $\tau_{\text{adm}}$ calibrated on 25 validation queries and tested on 25 in-range and 25 load-extrapolated queries, and a fallback to the full-order solve whenever $\eta_{\text{adm}} > \tau_{\text{adm}}$.

* **Experiment E6: Time-Dependent Contact Dynamics Extension (Optional, B1/B2; not carried out in this study)**
  To demonstrate applicability beyond static equilibrium, we integrate the reduced LCP with a second-order backward differentiation formula (BDF2) time-stepping scheme \cite{dabaghi2020posteriori, manvelyan2026inference}, evaluating $\eta_{\text{adm}}(t_n)$ at every discrete time step $t_n$ to detect transient impact admissibility violations.

#### 5.3 Quantitative Evaluation Metrics and Acceptance Criteria
To ensure objective assessment across all benchmarks, we predefine four quantitative validation metrics:
1. **Admissibility Detection Recall (Zero False Negatives):** Every test query with $\Vert R_{\text{red}}\Vert_\infty \le 10^{-12}$ whose true relative energy-norm error exceeds $5\%$ due to contact inadmissibility must satisfy $\eta_{\text{adm}} > \tau_{\text{adm}}$. *Scope:* $\tau_{\text{adm}}$ is calibrated per (basis, $r_\lambda$) on a held-out validation set as $0.9$ times the smallest indicator value among validation silent failures and then scored on separate test queries. The criterion is therefore an empirical, finite-sample statement, not a guarantee, and it cannot be evaluated when the validation set contains no silent failure.
2. **Effectivity Sharpness:** For the certified bound $\eta_{\text{bound}}$ (evaluated with the exact multiplier norm of the full-order solution, a validation quantity; the variant $\eta_{\text{bound}}^{C}$ with the offline constant $C_\lambda = \max_{\text{train}}\Vert\lambda\Vert_2$, usable online, is reported alongside \cite{hild2007residual, benaceur2020reduced, fontana2024guaranteed}), the effectivity index $\theta_{\text{eff}}$ lies within $[1, 10]$ on at least $90\%$ of admissible and borderline test cases, with explicit isolation of bulk-dominated regimes where $\eta_{\text{eq}}$ governs the error.
3. **Monotonic Error Correlation:** The Spearman rank correlation between $\eta_{\text{bound}}$ and the true energy-norm error $\Vert u - \tilde{u}\Vert_K$ exceeds $0.90$ across all dual bases and reduced dimensions. *Scope:* the minimum over all (basis, $r_\lambda$) entries is the strict reading and is reported together with the median; the correlation is computed over the test queries of each entry.
4. **Cost Crossover Accuracy:** The predicted break-even threshold $\mathcal{C}_{\text{ROM}} = \mathcal{C}_{\text{FOM}}$ matches the measured wall-clock break-even point within a factor of $2$ across benchmarks B1--B3.

#### 5.4 Results and Acceptance Outcome

The criteria were scored by `python -m contact_rom summary` (B1: spring chain, B2: membrane; five dual bases with rank repair, $r_\lambda \in \{4,8,12,16,24\}$, 40 test and 40 validation queries, seed 42). Theorem 1 ($\eta_{\text{bound}} \ge$ true energy error) held on every test query. Criterion 4 uses the E4 break-even sweeps on B2 and B3; E5 on B3 is reported separately. Outcomes are mixed and are reported as obtained, not tuned.

| Criterion | B1 | B2 |
|---|---|---|
| 1. Zero false negatives (calibrated $\tau_{\text{adm}}$, pooled over $r_\lambda$) | Met only by CPG (both indicators, 0/10) and RID with the consistent indicator $\eta_{\text{adm,v2}}$ (0/108). Missed: SVD 5/200 (1/200 with $\eta_{\text{adm,v2}}$), NMF 1/7, GIS 2/119, RID 3/108 with Eq. (15). Several NMF, CPG and RID entries have no validation silent failure, so no threshold exists. | Met by SVD (0/200), NMF (0/78) and GIS (0/188). CPG misses 1/115. RID has no silent failures (not evaluable). False-alarm rates on accurate queries are high: NMF 37--41%, CPG 55--58%, GIS 95--100%. |
| 2. $\theta_{\text{eff}} \in [1,10]$ on $\ge 90\%$ of queries above the 5% tolerance | Met: SVD and NMF 100% (200 and 7 queries). | Met: SVD and NMF 100% (200 and 78 queries). |
| 3. Spearman $> 0.90$ (min over $r_\lambda$; median in brackets) | Met by SVD 0.974 (0.978) and GIS 0.993 (0.996). Not met: NMF 0.515 (0.646), CPG 0.609 (0.686), RID $-0.637$ (0.995). | Met by SVD 0.992 (0.995) and GIS 0.953 (0.994). Not met: NMF 0.486 (0.682), CPG 0.283 (0.494), RID 0.811 (0.979). |
| 4. Break-even within a factor 2 | Not evaluated (fixed size, no sweep). | Membrane sweep, $N = 25$--$2304$: predicted vs measured $N^*$ is 119 vs 154 (ratio 0.78) in the affine regime and 64 vs 37 (ratio 1.72) when $K$ varies. Met in both. **B3 sweep** ($N = 24$--$16653$): constant $K$ 779 vs 1006 (ratio 0.77), met. Varying $K$: 55 vs 86 (ratio 0.64), met, but the measured speedup is not monotone below $N = 105$ and the crossover rests on one point, so this is weak evidence. |

**Seed robustness (seeds 42--46; entries are B1 / B2).** Criterion 1: seeds, out of five, with no missed silent failure (among seeds with at least one silent failure); criterion 3: seeds with minimum Spearman above 0.9.

| Basis | Crit. 3 (exact) | Crit. 3 ($C_\lambda$) | Crit. 1, 40 val., Eq. (15) | Crit. 1, 40 val., v2 | Crit. 1, 200 val., Eq. (15) | Crit. 1, 200 val., v2 |
|---|---|---|---|---|---|---|
| SVD | 5 / 5 | 5 / 5 | 3 / 4 | 2 / 4 | 5 / 5 | 5 / 5 |
| NMF | 0 / 0 | 0 / 3 | 3 / 1 | 3 / 2 | 5 / 3 | 5 / 5 |
| CPG | 0 / 0 | 0 / 0 | 3 / 3 | 4 / 3 | 5 / 4 | 5 / 5 |
| gIS | 2 / 5 | 2 / 5 | 0 / 1 | 0 / 1 | 3 / 3 | 4 / 4 |
| RID | 0 / 3 | 0 / 4 | 0 / 2 | 3 / 2 | 3 / 3 | 5 / 3 |

**Interpretation.**
- *Criterion 1 is not met in general, and its outcome depends on the seed.* Over five seeds (table below), with 40 validation queries the calibrated threshold missed at least one silent failure in most seeds for every basis (seeds without a miss, out of five: SVD 2--4, NMF 1--3, CPG 3--4, gIS 0--1, RID 0--3, depending on benchmark and indicator). Raising the validation set to 200 queries, with the 0.9 margin unchanged, removes most misses: the consistent indicator $\eta_{\text{adm,v2}}$ is then free of misses in all five seeds for SVD, NMF and CPG on both benchmarks, whereas gIS still misses in one seed on each benchmark (1/118 on B1, 3/162 on B2); Eq. (15) is free of misses in only three to five seeds. The false-alarm rate can be high (seed 42, B2: gIS 95--100%, CPG 55--58%). Many NMF, CPG and RID entries have no silent failure in the validation set, so no threshold exists for them (18--19 of 25 NMF entries on B1). Criterion 1 is therefore stated as a calibrated empirical recall with reported misses and false alarms; zero false negatives is not claimed. The main table keeps the original seed-42, 40-query protocol.
- *Criterion 3 is met reliably for SVD only.* Over five seeds the minimum Spearman correlation over dual dimensions exceeds 0.9 for SVD on both benchmarks in every seed. gIS passes on B2 in all five seeds but on B1 in only two (minimum between 0.02 and 0.99); RID passes on B2 in three seeds and on B1 in none; NMF (minimum 0.30--0.69) and CPG (0.28--0.61) never pass. The medians of gIS and RID are 0.99, so their failures come from a few dual dimensions. The bound is a valid upper bound (Theorem 1) but should not be used as a rank-preserving error predictor for NMF, CPG, RID or, on B1, gIS. The online-usable bound $\eta_{\text{bound}}^{C}$ with the offline constant $C_\lambda$ was at or above the true error on all test queries for SVD, NMF and CPG and on at least 99.5% of them for gIS and RID in every seed. Its rank correlation is higher for NMF and RID on B2 (minimum 0.87--0.93 for NMF, passing in three of five seeds; RID passing in four), but the conclusion is the same: only SVD, and gIS on B2, pass reliably. Its effectivity in the scoped range is 100% in [1, 10] except for NMF on B1 (86--100%). Restricting the criterion to the median or to selected bases would be a post hoc change and is not adopted.
- *Criterion 2 is met in its scoped form* (queries whose error exceeds the tolerance). Over all queries, effectivity grows like $\text{err}^{-1/2}$ as errors shrink and is a diagnostic only.
- *Criterion 4* is met on B2 (both regimes) and on B3 (constant $K$: 779 vs 1006, ratio 0.77; varying $K$: 55 vs 86, ratio 0.64); it is not evaluable for B1 (fixed size, no sweep). Two caveats apply. The B3 varying-stiffness break-even lies below $N = 100$, where the reduced dimension (about 20) is comparable to the model size and the sub-millisecond timings are dominated by fixed overheads; the measured speedup is not monotone there (1.24, 1.30, 0.78, 1.16 at $N = 24$, 40, 63, 105), so the measured break-even is set by a single point. Repeating the B3 sweep moved the measured constant-stiffness break-even from 703 to 1006, so timing uncertainty is of the order of the margin. The criterion concerns the break-even size only: the predicted speedup is much larger than measured at large $N$ (B3, $N = 16653$, constant $K$: 730 predicted vs 55 measured; varying $K$: 1138 vs 443), so the FLOP model should not be used to quote absolute speedups.
- *E5 (B3, default configuration, in-range vs extrapolated load).* With $\tau_{\text{adm}}$ calibrated on validation queries, the circuit breaker accepted no unsafe answer for NMF and CPG at $r_\lambda \in \{8,16\}$ (0 misses). On in-range queries about 4% fall back to the full-order solve and the speedup is about 16--18x. On extrapolated loads 96--100% fall back and the speedup is about 1x: the breaker keeps the answers safe but removes the acceleration. Only 4% of in-range queries were unsafe without the breaker, so the in-range zero-miss result rests on very few positives.

---

## v1 (original)

### 5. Numerical Experiments Protocol

To systematically evaluate the certified *a posteriori* estimator (Section 3), compare the five dual-basis constructions under controlled conditions (Section 4.1), and validate the FLOP- and iteration-aware break-even cost model (Section 4.2), we design a six-part experimental benchmark suite spanning 1D, 2D, and 3D unilateral contact problems. All reduced-order models are evaluated on identical snapshot sets $\{(u(\mu_i), \lambda(\mu_i))\}_{i=1}^{n_s}$ and solved using a unified active-set / Lemke LCP solver converged to a reduced KKT residual tolerance of $10^{-14}$.

#### 5.1 Benchmark Test Cases
We consider three progressively larger frictionless unilateral contact problems covering distinct interface topologies and sparsity structures:

1. **Benchmark B1: 1D Discrete Spring-Chain Contact ($N = 60, n_c = 30$)**
   A nonlinearly coupled or parametrized chain of $N = 60$ spring-mass nodes suspended above a rigid variable-profile obstacle $c(\mu)$ (\cite{pawar2025efficient}). The parameter vector $\mu = (\mu_1, \mu_2) \in \mathcal{P} \subset \mathbb{R}^2$ controls the spatially varying downward load amplitude and obstacle curvature. Because $n_c / N = 0.50$ is large relative to $N$, this benchmark isolates rank-repair overhead ($r_{\text{rep}} = 14$ unit vectors required for full row rank of $B_r$) in the constant-$K$ regime.
2. **Benchmark B2: 2D Parametrized Obstacle Membrane ($N = 324, n_c = 144$)**
   A clamped two-dimensional elastic membrane discretized on an $18 \times 18$ interior grid ($N = 324$) pressed against a rigid curved indenter (\cite{benaceur2020reduced}). The parameter vector $\mu = (F_0, x_c, y_c) \in \mathbb{R}^3$ varies the peak load magnitude $F_0$ and the indenter center $(x_c, y_c)$, causing the active contact patch $m_{\text{act}}$ to expand, contract, and translate across the candidate interface.
3. **Benchmark B3: 3D Two-Body Hertzian / Flange Contact ($N \in [2 \times 10^4, 10^5], n_c \in [400, 2500]$)**
   A three-dimensional linear-elastic deformable hemisphere (or bolted structural flange) pressed onto a deformable elastic block with a node-to-segment frictionless contact interface (\cite{franceschini2022block, koller2021hyper}). Full-order finite element matrices $(K(\mu), B, c(\mu), f(\mu))$ are assembled and verified using FEniCSx / GetFEM with a sparse direct Cholesky / PETSc TAO contact solver, and exported directly into the reduced-order pipeline. Parameters $\mu = (E_1/E_2, P_{\text{load}}, \delta_0)$ vary the stiffness contrast (testing parameter-varying $K(\mu)$), applied pressure, and initial clearance.

| Benchmark | Spatial Dim. | Full-Order $N$ | Contact Rows $n_c$ | Stiffness $K(\mu)$ | Primary Purpose |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **B1 (Spring Chain)** | 1D | 60 | 30 | Constant | Rank-repair cost crossover & silent SVD failure |
| **B2 (Obstacle Membrane)** | 2D | 324 | 144 | Constant | Five-basis comparison & contact-patch translation |
| **B3 (3D Two-Body Contact)** | 3D | $2 \times 10^4$--$10^5$ | 400--2500 | Parametrized $K(\mu)$ | Large-scale wall-clock speedup & 3D effectivity |

#### 5.2 Experimental Campaign (E1–E6)

* **Experiment E1: Silent Reduced-KKT Failure Audit (B1, B2)**
  Using identical primal POD bases $V_0$ and dual bases $U$ of matched dimension $r_\lambda \in \{4, 8, 12, 16, 20\}$, we log the reduced KKT residual $\|R_{\text{red}}(u_r, \lambda_r)\|_\infty$ alongside the true relative energy-norm error $\|u - \tilde{u}\|_K / \|u\|_K$, the four estimator terms $(\eta_{\text{eq}}, \eta_{\text{sign}}, \eta_{\text{gap}}, \eta_{\text{comp}})$, and the online admissibility indicator $\eta_{\text{adm}}$. This experiment verifies that $\eta_{\text{adm}}$ flags every case where $\|R_{\text{red}}\|_\infty \le 10^{-14}$ masks catastrophic lifted error.

* **Experiment E2: Estimator Effectivity and Extrapolation Stress Test (B1, B2, B3)**
  We evaluate the certified upper bound from Theorem 1 and compute the **effectivity index**
  $$\theta_{\text{eff}}(\mu) := \frac{\eta_{\text{bound}}(\mu)}{\|u(\mu) - \tilde{u}(\mu)\|_K}, \qquad \eta_{\text{bound}}(\mu) := \left( \eta_{\text{eq}}^2 + 2\,\eta_{\text{comp}} + 2\,\|\lambda_+\|_2\,\eta_{\text{gap}} \right)^{1/2},$$
  across three distinct parameter test sets:
  1. *In-sample validation set* $\mathcal{P}_{\text{val}} \subset \mathcal{P}_{\text{train}}$ (interpolatory regime);
  2. *Load-magnitude extrapolation* ($F_0 > \max_{\mu \in \mathcal{P}_{\text{train}}} F_0$), expanding the active contact set beyond training snapshots;
  3. *Contact-patch translation* (shifting $(x_c, y_c)$ so that previously inactive rows of $B$ activate online, testing whether $\eta_{\text{gap}}$ and $\eta_{\text{comp}}$ capture unseen active constraints).

* **Experiment E3: Head-to-Head Five-Basis Benchmark (B2, B3)**
  On identical snapshot matrices $(U_{\text{snap}}, \Lambda_{\text{snap}})$, we construct and compare all five dual-basis strategies from Section 4.1—**SVD+Repair**, **NMF** (\cite{balajewicz2016projection}), **CPG** (\cite{benaceur2020reduced}), **gIS** (\cite{sambataro2026model}), and **RID** (\cite{fauque2018hybrid, leberre2022condition})—across dual dimensions $r_\lambda \in \{4, 8, 12, 16, 24, 32\}$. Each method is evaluated along five Pareto axes:
  1. Lifted energy-norm displacement error and multiplier $\ell_2$-error;
  2. Online admissibility violation $\eta_{\text{adm}}$;
  3. Primal rank-repair enrichment $r_{\text{rep}}$ required to maintain $\text{cond}(B_r K_r^{-1} B_r^T) \le 10^6$;
  4. Active-set iteration count $k_{\text{ROM}}$ and online wall-clock time;
  5. Offline basis construction wall-clock time.

* **Experiment E4: Break-Even Cost Model Validation (B1, B2, B3)**
  We measure empirical floating-point operation counts and single-thread wall-clock execution times as functions of the enriched primal dimension $r_u = r_{u,0} + r_{\text{rep}}$, active contact size $m_{\text{act}}$, and mesh size $N$. Measured speedup curves $S_{\text{meas}} = t_{\text{FOM}} / t_{\text{ROM}}$ are plotted against the theoretical predictions $S_{\text{pred}} = \mathcal{C}_{\text{FOM}} / \mathcal{C}_{\text{ROM}}$ from Equations (4)–(5) in both the constant-$K$ (B1, B2) and parameter-varying $K(\mu)$ (B3) regimes.

* **Experiment E5: Large-Scale 3D Contact Certification (B3)**
  On the 3D two-body contact model ($N = 5 \times 10^4$ and $N = 10^5$), we deploy the top-performing non-negative dual bases from E3 paired with the online circuit breaker $\eta_{\text{adm}} \le \tau_{\text{adm}}$, demonstrating certified online acceleration and automatic fallback when extrapolation triggers $\eta_{\text{adm}} > \tau_{\text{adm}}$.

* **Experiment E6: Time-Dependent Contact Dynamics Extension (Optional, B1/B2)**
  To demonstrate applicability beyond static equilibrium, we integrate the reduced LCP with a second-order backward differentiation formula (BDF2) time-stepping scheme (\cite{dabaghi2020posteriori, manvelyan2026inference}), evaluating $\eta_{\text{adm}}(t_n)$ at every discrete time step $t_n$ to detect transient impact admissibility violations.

#### 5.3 Quantitative Evaluation Metrics and Acceptance Criteria
To ensure objective assessment across all benchmarks, we predefine four quantitative validation metrics:
1. **Admissibility Detection Recall (Zero False Negatives):** Every test query with $\|R_{\text{red}}\|_\infty \le 10^{-12}$ whose true relative energy-norm error exceeds $5\%$ due to contact inadmissibility must satisfy $\eta_{\text{adm}} > \tau_{\text{adm}}$.
2. **Effectivity Sharpness:** For the certified bound $\eta_{\text{bound}}$ (using $\|\lambda_+\|_2$ as the online multiplier norm proxy), the effectivity index $\theta_{\text{eff}}$ lies within $[1, 10]$ on at least $90\%$ of admissible and borderline test cases, with explicit isolation of bulk-dominated regimes where $\eta_{\text{eq}}$ governs the error.
3. **Monotonic Error Correlation:** The Spearman rank correlation between $\eta_{\text{bound}}$ and the true energy-norm error $\|u - \tilde{u}\|_K$ exceeds $0.90$ across all dual bases and reduced dimensions.
4. **Cost Crossover Accuracy:** The predicted break-even threshold $\mathcal{C}_{\text{ROM}} = \mathcal{C}_{\text{FOM}}$ matches the measured wall-clock break-even point within a factor of $2$ across benchmarks B1–B3.
