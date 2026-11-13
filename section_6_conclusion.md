# Section 6: Conclusion — All Versions

---

## v2 (current)

### 6. Conclusion and Practical Guidelines

#### 6.1 Summary of Contributions
Projection-based reduced-order models for unilateral contact problems face a fundamental reliability hazard: a reduced solver can converge to its projected Karush--Kuhn--Tucker (KKT) conditions to machine precision ($10^{-16}$) while returning a lifted displacement--multiplier pair $(\tilde{u}, \tilde{\lambda}) = (V u_r, U \lambda_r)$ that violates physical non-penetration and traction non-negativity by orders of magnitude \cite{balajewicz2016projection, benaceur2020reduced}. In this work, we addressed this gap through three interconnected contributions:

1. **Computable Lifted Admissibility Estimator:** Adapting constant-sign multiplier residual decompositions \cite{hild2007residual, dabaghi2020posteriori, banz2024generic, fontana2024guaranteed, banz2026posteriori} to lifted reduced-order pairs under general rectangular constraint matrices $B(\mu) \in \mathbb{R}^{n_c \times N}$ (removing the bijective operator restriction of \cite{zhang2014duality}), we proved an energy-norm *a posteriori* upper bound composed of a dual-norm equilibrium residual ($\eta_{\text{eq}}$) and three contact-interface terms measuring multiplier sign violation ($\eta_{\text{sign}}$), primal interpenetration ($\eta_{\text{gap}}$), and complementarity defect ($\eta_{\text{comp}}$). The resulting online admissibility indicator $\eta_{\text{adm}} = \eta_{\text{sign}} + \eta_{\text{gap}} + \sqrt{\eta_{\text{comp}}}$ evaluates in $\mathcal{O}(n_c (r_u + r_\lambda))$ operations without touching the bulk mesh, acting as an effective online circuit breaker against silent reduced-KKT convergence.
2. **Controlled Cross-Basis Benchmark:** Evaluating signed SVD with primal rank repair \cite{leberre2022condition}, non-negative matrix factorization (NMF) \cite{balajewicz2016projection}, cone-projected greedy (CPG) \cite{benaceur2020reduced}, greedy index selection (gIS) \cite{sambataro2026model, kollepara2024sparse}, and reduced integration domain (RID) restriction \cite{fauque2018hybrid, leberre2022condition} on identical snapshot sets clarifies the trade-off between offline construction cost, dual non-negativity, primal rank-repair overhead, and active-set extrapolation robustness.
3. **FLOP- and Iteration-Aware Break-Even Cost Model:** By accounting for dense reduced Schur-complement factorization, unit-vector rank-repair enrichment $r_{\text{rep}}$ \cite{leberre2022condition}, and active-set iteration growth $k_{\text{ROM}}$ relative to sparse full-order solvers \cite{franceschini2022block, pawar2025efficient}, our cost model predicts the exact problem dimensions and active contact-zone fractions where Lagrange-multiplier contact ROMs yield genuine wall-clock speedups.

#### 6.2 Practical Guidelines for Contact ROM Deployment
Based on the theoretical bounds in Section 3 and the cost scaling in Section 4, we summarize four actionable guidelines for practitioners deploying Lagrange-multiplier reduced contact models:

1. **Never rely solely on reduced KKT residuals for online acceptance:**
   Always evaluate the $\mathcal{O}(n_c (r_u + r_\lambda))$ online admissibility indicator $\eta_{\text{adm}}$ after lifting $(u_r, \lambda_r)$. Whenever $\eta_{\text{adm}}$ exceeds the target physical tolerance, flag the reduced prediction and trigger either primal rank-repair enrichment \cite{leberre2022condition} or a fallback solve with the full-order model \cite{dabaghi2020posteriori}.
2. **Enforce dual non-negativity by construction ($\eta_{\text{sign}} = 0$) whenever possible:**
   Avoid truncated signed SVD dual bases unless paired with strict online $\eta_{\text{sign}}$ monitoring and aggressive interface enrichment. When offline NMF cost is acceptable, NMF \cite{balajewicz2016projection} provides compact non-negative dual bases; when fast offline assembly or exact pointwise complementarity at key contact nodes is prioritized, greedy index selection (gIS) \cite{sambataro2026model, kollepara2024sparse} and cone-projected greedy (CPG) \cite{benaceur2020reduced} eliminate sign violations with low offline overhead.
3. **Audit primal rank-repair enrichment against active contact size before deployment:**
   In constant-stiffness problems where $B K^{-1} B^T$ can be precomputed offline \cite{pawar2025efficient}, dense reduced algebra scales cubically with the enriched primal dimension $r_u = r_{u,0} + r_{\text{rep}}$ \cite{leberre2022condition}. Verify via the break-even inequality in Section 4.2 that $r_u$ remains substantially smaller than the active contact count $m_{\text{act}}$; otherwise, a reduced-space active-set solver on the precomputed full-order compliance matrix \cite{franceschini2022block} will outperform the projection ROM.
4. **Use $\eta_{\text{eq}}$ offline to enrich bulk primal modes during greedy training:**
   Because $\eta_{\text{adm}}$ monitors interface admissibility alone, use the full four-component bound $\eta_{\text{bound}}$ (including $\eta_{\text{eq}}$) during offline greedy sampling \cite{benaceur2020reduced, zhang2014duality} so that interior bulk deformation modes are enriched alongside contact interface modes.

#### 6.3 Limitations and Future Work
Two extensions remain open for future investigation. First, while our experiments focus primarily on frictionless unilateral contact, extending the constant-sign multiplier decomposition to Coulomb frictional contact \cite{niakh2024reduced} requires coupling normal traction admissibility with tangential friction-cone projection residuals. Second, integrating the online circuit breaker $\eta_{\text{adm}}$ into an automated, on-the-fly basis enrichment loop for large-scale 3D nonlinear elastoplastic and substructured dynamic contact \cite{koller2021hyper, aresdeparga2023petrov, manvelyan2026inference} represents a natural next step toward certified industrial contact digital twins.

---

## v1 (original)

## 6. Conclusion and Practical Guidelines

### 6.1 Summary of Contributions
Projection-based reduced-order models for unilateral contact problems face a fundamental reliability hazard: a reduced solver can converge to its projected Karush–Kuhn–Tucker (KKT) conditions to machine precision ($10^{-16}$) while returning a lifted displacement–multiplier pair $(\tilde{u}, \tilde{\lambda}) = (V u_r, U \lambda_r)$ that violates physical non-penetration and traction non-negativity by orders of magnitude. In this work, we addressed this gap through three interconnected contributions:

1. **Computable Lifted Admissibility Estimator:** Adapting constant-sign multiplier residual decompositions to lifted reduced-order pairs under general rectangular constraint matrices $B(\mu) \in \mathbb{R}^{n_c \times N}$, we proved an energy-norm *a posteriori* upper bound composed of a dual-norm equilibrium residual ($\eta_{\text{eq}}$) and three contact-interface terms measuring multiplier sign violation ($\eta_{\text{sign}}$), primal interpenetration ($\eta_{\text{gap}}$), and complementarity defect ($\eta_{\text{comp}}$). The resulting online admissibility indicator $\eta_{\text{adm}} = \eta_{\text{sign}} + \eta_{\text{gap}} + \sqrt{\eta_{\text{comp}}}$ evaluates in $\mathcal{O}(n_c (r_u + r_\lambda))$ operations without touching the bulk mesh, acting as an effective online circuit breaker against silent reduced-KKT convergence.
2. **Controlled Cross-Basis Benchmark:** Evaluating signed SVD with primal rank repair, non-negative matrix factorization (NMF), cone-projected greedy (CPG), greedy index selection (gIS), and reduced integration domain (RID) restriction on identical snapshot sets clarifies the trade-off between offline construction cost, dual non-negativity, primal rank-repair overhead, and active-set extrapolation robustness.
3. **FLOP- and Iteration-Aware Break-Even Cost Model:** By accounting for dense reduced Schur-complement factorization, unit-vector rank-repair enrichment $r_{\text{rep}}$, and active-set iteration growth $k_{\text{ROM}}$ relative to sparse full-order solvers, our cost model predicts the exact problem dimensions and active contact-zone fractions where Lagrange-multiplier contact ROMs yield genuine wall-clock speedups.

### 6.2 Practical Guidelines for Contact ROM Deployment
Based on the theoretical bounds in Section 3 and the cost scaling in Section 4, we summarize four actionable guidelines for practitioners deploying Lagrange-multiplier reduced contact models:

1. **Never rely solely on reduced KKT residuals for online acceptance:**
   Always evaluate the $\mathcal{O}(n_c (r_u + r_\lambda))$ online admissibility indicator $\eta_{\text{adm}}$ after lifting $(u_r, \lambda_r)$. Whenever $\eta_{\text{adm}}$ exceeds the target physical tolerance, flag the reduced prediction and trigger either primal rank-repair enrichment or a fallback solve with the full-order model.
2. **Enforce dual non-negativity by construction ($\eta_{\text{sign}} = 0$) whenever possible:**
   Avoid truncated signed SVD dual bases unless paired with strict online $\eta_{\text{sign}}$ monitoring and aggressive interface enrichment. When offline NMF cost is acceptable, NMF provides compact non-negative dual bases; when fast offline assembly or exact pointwise complementarity at key contact nodes is prioritized, greedy index selection (gIS) and cone-projected greedy (CPG) eliminate sign violations with low offline overhead.
3. **Audit primal rank-repair enrichment against active contact size before deployment:**
   In constant-stiffness problems where $B K^{-1} B^T$ can be precomputed offline, dense reduced algebra scales cubically with the enriched primal dimension $r_u = r_{u,0} + r_{\text{rep}}$. Verify via the break-even inequality in Section 4.2 that $r_u$ remains substantially smaller than the active contact count $m_{\text{act}}$; otherwise, a reduced-space active-set solver on the precomputed full-order compliance matrix will outperform the projection ROM.
4. **Use $\eta_{\text{eq}}$ offline to enrich bulk primal modes during greedy training:**
   Because $\eta_{\text{adm}}$ monitors interface admissibility alone, use the full four-component bound $\eta_{\text{bound}}$ (including $\eta_{\text{eq}}$) during offline greedy sampling so that interior bulk deformation modes are enriched alongside contact interface modes.

### 6.3 Limitations and Future Work
Two extensions remain open for future investigation. First, while our experiments focus primarily on frictionless unilateral contact, extending the constant-sign multiplier decomposition to Coulomb frictional contact requires coupling normal traction admissibility with tangential friction-cone projection residuals. Second, integrating the online circuit breaker $\eta_{\text{adm}}$ into an automated, on-the-fly basis enrichment loop for large-scale 3D nonlinear elastoplastic contact represents a natural next step toward certified industrial contact digital twins.
