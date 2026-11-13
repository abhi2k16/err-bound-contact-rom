# Section 4: Basis and Cost — All Versions

---

## v2 (current)

### 4. Dual-Basis Constructions and Break-Even Cost Model

Sections 2 and 3 established that the reliability of a Lagrange-multiplier contact ROM hinges on the interplay between the dual basis $U \in \mathbb{R}^{n_c \times r_\lambda}$, the primal enrichment required to keep $B_r(\mu) = U^T B(\mu) V$ well-conditioned \cite{benaceur2020reduced, leberre2022condition}, and the resulting online computational cost. In this section, we first unify the five major dual-basis constructions within a common algebraic notation and then formulate an explicit FLOP- and iteration-aware break-even cost model comparing the enriched ROM against a sparse full-order contact solver \cite{franceschini2022block, pawar2025efficient}.

#### 4.1 Unified Summary of Dual-Basis Constructions
Let $\Lambda_{\text{snap}} = [\lambda(\mu_1), \dots, \lambda(\mu_{n_s})] \in \mathbb{R}_{\ge 0}^{n_c \times n_s}$ denote the matrix of non-negative multiplier snapshots collected over $\mathcal{P}_{\text{train}}$, and let $V_0 \in \mathbb{R}^{N \times r_{u,0}}$ be an initial orthonormal POD basis of the displacement snapshots $U_{\text{snap}} = [u(\mu_1), \dots, u(\mu_{n_s})]$ \cite{balajewicz2016projection, benaceur2020reduced}. Each of the five dual-basis strategies constructs $U \in \mathbb{R}^{n_c \times r_\lambda}$---and, where necessary, enriches $V_0$ to $V \in \mathbb{R}^{N \times r_u}$---as follows:

1. **Signed SVD with Primal Rank Repair (SVD+Repair):**
   The dual basis $U \in \mathbb{R}^{n_c \times r_\lambda}$ is taken as the leading $r_\lambda$ left singular vectors of $\Lambda_{\text{snap}}$ \cite{balajewicz2016projection, leberre2022condition}. Because $U$ contains mixed signs, $\eta_{\text{sign}} \ge 0$ is not controlled by construction. Moreover, $U^T B V_0$ frequently exhibits rank deficiency or severe ill-conditioning on active contact rows. To restore full row rank of $B_r$, $V_0$ is enriched with $r_{\text{rep}}$ unit vectors (or supremizer modes $\{K^{-1} B^T u_j\}_{j=1}^{r_\lambda}$) localized at poorly resolved contact degrees of freedom, yielding $V = [V_0, V_{\text{rep}}]$ with $r_u = r_{u,0} + r_{\text{rep}}$ \cite{benaceur2020reduced, leberre2022condition}.

2. **Non-Negative Matrix Factorization (NMF):**
   Following Balajewicz et al. \cite{balajewicz2016projection}, $\Lambda_{\text{snap}} \approx U H$ is factored offline via projected alternating least squares or multiplicative updates subject to $U \in \mathbb{R}_{\ge 0}^{n_c \times r_\lambda}$ and $H \in \mathbb{R}_{\ge 0}^{r_\lambda \times n_s}$. Because every entry of $U$ is non-negative and $\lambda_r \ge 0$ is enforced by the reduced KKT solver, the lifted multiplier satisfies $\tilde{\lambda} = U \lambda_r \ge 0$ identically, guaranteeing $\eta_{\text{sign}} = 0$ for all parameters. However, NMF columns are non-orthogonal and generally less compact than SVD modes, often requiring larger $r_\lambda$ and primal enrichment $r_{\text{rep}}$ to prevent locking \cite{balajewicz2016projection, leberre2022condition}.

3. **Cone-Projected Greedy Basis (CPG):**
   Following Benaceur et al. \cite{benaceur2020reduced}, the dual basis is built iteratively by selecting snapshots from $\Lambda_{\text{snap}}$ that maximize a projection residual onto the positive cone $\text{cone}(U^{(k-1)}) = \{ \sum_{j=1}^{k-1} \alpha_j u_j : \alpha_j \ge 0 \}$. Specifically, at iteration $k$, the snapshot $\lambda(\mu^*)$ maximizing $\min_{\alpha \ge 0} \Vert{}\lambda(\mu) - U^{(k-1)} \alpha\Vert{}_2$ is added directly as a normalized column of $U$. Because $U$ consists of positively scaled physical snapshots, $\eta_{\text{sign}} = 0$ holds identically without iterative NMF optimization, though linear independence can degrade as $r_\lambda$ grows \cite{benaceur2020reduced}.

4. **Greedy Index Selection (gIS):**
   Following Sambataro and Ehrlacher \cite{sambataro2026model} (and related sparse contact selection approaches \cite{kollepara2024sparse}), rather than forming dense linear combinations of snapshots, $U$ is chosen as a sparse selection matrix $U = [e_{i_1}, \dots, e_{i_{r_\lambda}}] \in \{0, 1\}^{n_c \times r_\lambda}$ formed from canonical unit vectors at $r_\lambda$ greedily selected contact constraint indices $\mathcal{I} = \{i_1, \dots, i_{r_\lambda}\} \subset \{1, \dots, n_c\}$. This guarantees $\eta_{\text{sign}} = 0$, preserves exact pointwise complementarity $(B \tilde{u} - c)_{i_k} = 0$ on the selected indices $\mathcal{I}$, and eliminates dense dual projection overhead, at the cost of setting $\tilde{\lambda}_i = 0$ for all $i \notin \mathcal{I}$ \cite{sambataro2026model}.

5. **Reduced Integration Domain Restriction (RID):**
   Following Fauque et al. \cite{fauque2018hybrid} and Le Berre et al. \cite{leberre2022condition} (see also empirical cubature surface reduction \cite{koller2021hyper, aresdeparga2023petrov}), the contact interface and bulk integration are restricted to a reduced domain $\Omega_{\text{RID}}$ retaining $n_{c,\text{RID}} \le n_c$ contact rows around high-activity regions, paired with localized primal interface enrichment $V_{\text{RID}}$ (unit vectors on $\Omega_{\text{RID}}$ nodes). Within $\Omega_{\text{RID}}$, pointwise multipliers or localized positive modes enforce $\eta_{\text{sign}} = 0$ on retained rows, while efficiency depends directly on keeping $n_{c,\text{RID}} \ll n_c$ \cite{fauque2018hybrid, leberre2022condition}.

| Dual-Basis Method | Guarantees $\eta_{\text{sign}} = 0$? | Orthogonal $U$? | Offline Build Cost | Typical Primal Enrichment $r_{\text{rep}}$ |
| :--- | :--- | :--- | :--- | :--- |
| **Signed SVD + Repair** \cite{leberre2022condition} | No | Yes | Low (single SVD) | Moderate--High (unit/supremizer modes) |
| **NMF** \cite{balajewicz2016projection} | Yes | No | High (iterative NMF) | Moderate |
| **Cone-Projected Greedy** \cite{benaceur2020reduced} | Yes | No | Moderate (NNLS greedy) | Moderate |
| **Greedy Index Selection** \cite{sambataro2026model, kollepara2024sparse} | Yes | Yes (columns of $I_{n_c}$) | Low--Moderate (index greedy) | Low--Moderate |
| **RID Restriction** \cite{fauque2018hybrid, leberre2022condition} | Yes (on $\Omega_{\text{RID}}$) | Problem-dependent | Moderate (domain + POD) | High (interface nodal enrichment) |

#### 4.2 FLOP- and Iteration-Aware Break-Even Cost Model
A common misconception in reduced-order modeling for contact is that $r_u + r_\lambda \ll N$ automatically guarantees online speedup. In practice, two structural effects erode the computational advantage of the ROM:
1. **Dense vs. Sparse Linear Algebra:** The full-order stiffness matrix $K \in \mathbb{R}^{N \times N}$ is sparse with bandwidth $b_w$ (or fill-in factor $\alpha \in [1.5, 2.0]$ under sparse Cholesky / algebraic multigrid preconditioning \cite{franceschini2022block}), whereas the reduced saddle-point and Schur-complement operators $K_r \in \mathbb{R}^{r_u \times r_u}$ and $B_r K_r^{-1} B_r^T \in \mathbb{R}^{r_\lambda \times r_\lambda}$ are **dense** \cite{balajewicz2016projection, leberre2022condition}.
2. **Primal Rank-Repair & Active-Set Iteration Growth:** Adding $r_{\text{rep}}$ unit-vector or supremizer columns increases the effective primal dimension to $r_u = r_{u,0} + r_{\text{rep}}$ \cite{leberre2022condition}. Simultaneously, poorly conditioned dual bases (such as unorthogonalized NMF or CPG columns) often increase the number of active-set or Lemke LCP iterations $k_{\text{ROM}}$ required for convergence relative to the full-order iteration count $k_{\text{FOM}}$ \cite{leberre2022condition, pawar2025efficient}.

To quantify the exact computational crossover, we model the online floating-point operation (FLOP) counts of both solvers for a single parameter query $\mu \in \mathcal{P}$:

##### Full-Order Online Cost ($\mathcal{C}_{\text{FOM}}$)
Assuming a sparse direct Cholesky factorization of $K$ (costing $c_{\text{fact}} N b_w^2$ FLOPs, or $c_{\text{fact}} N^\alpha$ in 3D \cite{franceschini2022block}) followed by $n_c$ forward/backward substitutions to form the dense contact compliance matrix $B K^{-1} B^T \in \mathbb{R}^{n_c \times n_c}$ and $k_{\text{FOM}}$ active-set updates of average active-set size $m_{\text{act}} \le n_c$ \cite{balajewicz2016projection, leberre2022condition}:
$$\mathcal{C}_{\text{FOM}} = \underbrace{c_{\text{fact}} N b_w^2}_{\text{sparse } K \text{ factorization}} + \underbrace{2 c_{\text{solve}} n_c N b_w}_{\text{form } B K^{-1} B^T} + \underbrace{k_{\text{FOM}} \left( \frac{2}{3} m_{\text{act}}^3 + 2 n_c m_{\text{act}} \right)}_{\text{full-order active-set / LCP solves}}. \tag{4}$$
*(Note: When $K$ is parameter-independent, the first two terms are amortized offline \cite{pawar2025efficient}, leaving only the active-set LCP cost $\mathcal{C}_{\text{FOM}}^{\text{affine}} = k_{\text{FOM}} (\frac{2}{3} m_{\text{act}}^3 + 2 n_c m_{\text{act}})$; when $K(\mu)$ varies with $\mu$, all three terms execute online.)*

##### Reduced-Order Online Cost ($\mathcal{C}_{\text{ROM}}$)
Let $r_u = r_{u,0} + r_{\text{rep}}$ be the enriched primal dimension and $n_{c,\text{eff}} \le n_c$ be the number of online contact rows evaluated ($n_{c,\text{eff}} = n_c$ for global dual bases \cite{balajewicz2016projection, benaceur2020reduced}, $n_{c,\text{eff}} = r_\lambda$ for gIS without global certification \cite{sambataro2026model}, or $n_{c,\text{eff}} = n_{c,\text{RID}}$ for RID \cite{fauque2018hybrid}). Forming the dense reduced operators $K_r(\mu)$ and $B_r(\mu)$ (via affine decomposition or hyper-reduction \cite{benaceur2020reduced, koller2021hyper}), factoring dense $K_r \in \mathbb{R}^{r_u \times r_u}$, computing the dense Schur complement $B_r K_r^{-1} B_r^T \in \mathbb{R}^{r_\lambda \times r_\lambda}$, running $k_{\text{ROM}}$ active-set iterations of average reduced active size $m_{r,\text{act}} \le r_\lambda$ \cite{leberre2022condition}, and evaluating the online admissibility indicator $\eta_{\text{adm}}$ (Section 3.3) require:
$$\mathcal{C}_{\text{ROM}} = \underbrace{\mathcal{C}_{\text{proj}}(N_{\text{eff}}, r_u, r_\lambda)}_{\text{reduced operator assembly}} + \underbrace{\frac{1}{3} r_u^3 + 2 r_\lambda r_u^2 + 2 r_\lambda^2 r_u}_{\text{dense } K_r^{-1} \text{ \& } B_r K_r^{-1} B_r^T} + \underbrace{k_{\text{ROM}} \left( \frac{2}{3} m_{r,\text{act}}^3 + 2 r_\lambda m_{r,\text{act}} \right)}_{\text{reduced active-set / LCP solves}} + \underbrace{2 n_c (r_u + r_\lambda)}_{\text{online } \eta_{\text{adm}} \text{ check}}, \tag{5}$$
where $\mathcal{C}_{\text{proj}} = \mathcal{O}(Q_K r_u^2 + Q_B r_u r_\lambda)$ for affine problems with $Q_K, Q_B$ affine terms \cite{benaceur2020reduced}, and $\mathcal{C}_{\text{proj}} = \mathcal{O}(\text{nnz}(K) r_u + \text{nnz}(B) (r_u + r_\lambda))$ when full-order matrices are projected online \cite{koller2021hyper, aresdeparga2023petrov}.

##### Break-Even Conditions
Defining the theoretical speedup ratio $S := \mathcal{C}_{\text{FOM}} / \mathcal{C}_{\text{ROM}}$, the reduced contact model achieves computational break-even ($S > 1$) if and only if $\mathcal{C}_{\text{ROM}} < \mathcal{C}_{\text{FOM}}$. Equations (4)--(5) reveal two distinct break-even regimes tested in Section 5:

1. **Parameter-Varying Stiffness Regime ($K = K(\mu)$):**
   Because $\mathcal{C}_{\text{FOM}}$ is dominated by the $\mathcal{O}(N b_w^2)$ sparse factorization and $\mathcal{O}(n_c N b_w)$ compliance assembly \cite{franceschini2022block}, break-even holds whenever
   $$r_u^3 + 6 r_\lambda r_u^2 + 2 k_{\text{ROM}} m_{r,\text{act}}^3 \ll 3 c_{\text{fact}} N b_w^2 + 6 c_{\text{solve}} n_c N b_w.$$
   In this regime, substantial speedups ($S \sim 10^1\text{--}10^3$) persist even with moderate rank-repair enrichment $r_{\text{rep}}$, provided $r_u \ll (N b_w^2)^{1/3}$ \cite{koller2021hyper, pawar2025efficient}.
2. **Parameter-Independent Stiffness Regime ($K$ constant, varying loads/gaps $f(\mu), c(\mu)$):**
   When $B K^{-1} B^T$ is precomputed offline for the full-order solver \cite{pawar2025efficient}, $\mathcal{C}_{\text{FOM}}^{\text{affine}}$ depends solely on the active contact count $m_{\text{act}}$. Break-even then reduces to the sharp interface-compression condition
   $$r_u^3 + 6 r_\lambda r_u^2 + 2 k_{\text{ROM}} m_{r,\text{act}}^3 + 6 n_c (r_u + r_\lambda) < 2 k_{\text{FOM}} m_{\text{act}}^3.$$
   Consequently, if the physical contact patch $m_{\text{act}}$ is small (e.g., localized point or line contact with $m_{\text{act}} \lesssim 30$) and rank repair requires adding $r_{\text{rep}} \sim 15\text{--}30$ unit vectors (such as $r_{\text{rep}} = 14$ on the 60-spring chain \cite{leberre2022condition}), the cubic dense cost $r_u^3$ immediately exceeds $2 k_{\text{FOM}} m_{\text{act}}^3$, causing $S < 1$. Thus, **Lagrange-multiplier contact ROMs are computationally advantageous in the constant-$K$ regime only when the active contact zone $m_{\text{act}}$ is large relative to the enriched reduced dimension $r_u = r_{u,0} + r_{\text{rep}}$.**

---

## v1 (original)

## 4. Dual-Basis Constructions and Break-Even Cost Model

Sections 2 and 3 established that the reliability of a Lagrange-multiplier contact ROM hinges on the interplay between the dual basis $U \in \mathbb{R}^{n_c \times r_\lambda}$, the primal enrichment required to keep $B_r(\mu) = U^T B(\mu) V$ well-conditioned, and the resulting online computational cost. In this section, we first unify the five major dual-basis constructions within a common algebraic notation and then formulate an explicit FLOP- and iteration-aware break-even cost model comparing the enriched ROM against a sparse full-order contact solver.

### 4.1 Unified Summary of Dual-Basis Constructions
Let $\Lambda_{\text{snap}} = [\lambda(\mu_1), \dots, \lambda(\mu_{n_s})] \in \mathbb{R}_{\ge 0}^{n_c \times n_s}$ denote the matrix of non-negative multiplier snapshots collected over $\mathcal{P}_{\text{train}}$, and let $V_0 \in \mathbb{R}^{N \times r_{u,0}}$ be an initial orthonormal POD basis of the displacement snapshots $U_{\text{snap}} = [u(\mu_1), \dots, u(\mu_{n_s})]$. Each of the five dual-basis strategies constructs $U \in \mathbb{R}^{n_c \times r_\lambda}$—and, where necessary, enriches $V_0$ to $V \in \mathbb{R}^{N \times r_u}$—as follows:

1. **Signed SVD with Primal Rank Repair (SVD+Repair):**
   The dual basis $U \in \mathbb{R}^{n_c \times r_\lambda}$ is taken as the leading $r_\lambda$ left singular vectors of $\Lambda_{\text{snap}}$. Because $U$ contains mixed signs, $\eta_{\text{sign}} \ge 0$ is not controlled by construction. Moreover, $U^T B V_0$ frequently exhibits rank deficiency or severe ill-conditioning on active contact rows. To restore full row rank of $B_r$, $V_0$ is enriched with $r_{\text{rep}}$ unit vectors (or supremizer modes $\{K^{-1} B^T u_j\}_{j=1}^{r_\lambda}$) localized at poorly resolved contact degrees of freedom, yielding $V = [V_0, V_{\text{rep}}]$ with $r_u = r_{u,0} + r_{\text{rep}}$ \cite{leberre2022condition}.

2. **Non-Negative Matrix Factorization (NMF):**
   Following Balajewicz et al. \cite{balajewicz2016projection}, $\Lambda_{\text{snap}} \approx U H$ is factored offline via projected alternating least squares or multiplicative updates subject to $U \in \mathbb{R}_{\ge 0}^{n_c \times r_\lambda}$ and $H \in \mathbb{R}_{\ge 0}^{r_\lambda \times n_s}$. Because every entry of $U$ is non-negative and $\lambda_r \ge 0$ is enforced by the reduced KKT solver, the lifted multiplier satisfies $\tilde{\lambda} = U \lambda_r \ge 0$ identically, guaranteeing $\eta_{\text{sign}} = 0$ for all parameters. However, NMF columns are non-orthogonal and generally less compact than SVD modes, often requiring larger $r_\lambda$ and primal enrichment $r_{\text{rep}}$ to prevent locking.

3. **Cone-Projected Greedy Basis (CPG):**
   Following Benaceur et al. \cite{benaceur2020reduced}, the dual basis is built iteratively by selecting snapshots from $\Lambda_{\text{snap}}$ that maximize a projection residual onto the positive cone $\text{cone}(U^{(k-1)}) = \{ \sum_{j=1}^{k-1} \alpha_j u_j : \alpha_j \ge 0 \}$. Specifically, at iteration $k$, the snapshot $\lambda(\mu^*)$ maximizing $\min_{\alpha \ge 0} \|\lambda(\mu) - U^{(k-1)} \alpha\|_2$ is added directly as a normalized column of $U$. Because $U$ consists of positively scaled physical snapshots, $\eta_{\text{sign}} = 0$ holds identically without iterative NMF optimization, though linear independence can degrade as $r_\lambda$ grows.

4. **Greedy Index Selection (gIS):**
   Following Sambataro and Ehrlacher \cite{sambataro2026model}, rather than forming dense linear combinations of snapshots, $U$ is chosen as a sparse selection matrix $U = [e_{i_1}, \dots, e_{i_{r_\lambda}}] \in \{0, 1\}^{n_c \times r_\lambda}$ formed from canonical unit vectors at $r_\lambda$ greedily selected contact constraint indices $\mathcal{I} = \{i_1, \dots, i_{r_\lambda}\} \subset \{1, \dots, n_c\}$. This guarantees $\eta_{\text{sign}} = 0$, preserves exact pointwise complementarity $(B \tilde{u} - c)_{i_k} = 0$ on the selected indices $\mathcal{I}$, and eliminates dense dual projection overhead, at the cost of setting $\tilde{\lambda}_i = 0$ for all $i \notin \mathcal{I}$.

5. **Reduced Integration Domain Restriction (RID):**
   Following Fauque et al. \cite{fauque2018hybrid} and Le Berre et al. \cite{leberre2022condition}, the contact interface and bulk integration are restricted to a reduced domain $\Omega_{\text{RID}}$ retaining $n_{c,\text{RID}} \le n_c$ contact rows around high-activity regions, paired with localized primal interface enrichment $V_{\text{RID}}$ (unit vectors on $\Omega_{\text{RID}}$ nodes). Within $\Omega_{\text{RID}}$, pointwise multipliers or localized positive modes enforce $\eta_{\text{sign}} = 0$ on retained rows, while efficiency depends directly on keeping $n_{c,\text{RID}} \ll n_c$.

| Dual-Basis Method | Guarantees $\eta_{\text{sign}} = 0$? | Orthogonal $U$? | Offline Build Cost | Typical Primal Enrichment $r_{\text{rep}}$ |
| :--- | :--- | :--- | :--- | :--- |
| **Signed SVD + Repair** | No | Yes | Low (single SVD) | Moderate–High (unit/supremizer modes) |
| **NMF** \cite{balajewicz2016projection} | Yes | No | High (iterative NMF) | Moderate |
| **Cone-Projected Greedy** \cite{benaceur2020reduced} | Yes | No | Moderate (NNLS greedy) | Moderate |
| **Greedy Index Selection** \cite{sambataro2026model} | Yes | Yes (columns of $I_{n_c}$) | Low–Moderate (index greedy) | Low–Moderate |
| **RID Restriction** \cite{fauque2018hybrid, leberre2022condition} | Yes (on $\Omega_{\text{RID}}$) | Problem-dependent | Moderate (domain + POD) | High (interface nodal enrichment) |

### 4.2 FLOP- and Iteration-Aware Break-Even Cost Model
A common misconception in reduced-order modeling for contact is that $r_u + r_\lambda \ll N$ automatically guarantees online speedup. In practice, two structural effects erode the computational advantage of the ROM:
1. **Dense vs. Sparse Linear Algebra:** The full-order stiffness matrix $K \in \mathbb{R}^{N \times N}$ is sparse with bandwidth $b_w$ (or fill-in factor $\alpha \in [1.5, 2.0]$ under sparse Cholesky / algebraic multigrid preconditioning \cite{franceschini2022block}), whereas the reduced saddle-point and Schur-complement operators $K_r \in \mathbb{R}^{r_u \times r_u}$ and $B_r K_r^{-1} B_r^T \in \mathbb{R}^{r_\lambda \times r_\lambda}$ are **dense**.
2. **Primal Rank-Repair & Active-Set Iteration Growth:** Adding $r_{\text{rep}}$ unit-vector or supremizer columns increases the effective primal dimension to $r_u = r_{u,0} + r_{\text{rep}}$. Simultaneously, poorly conditioned dual bases (such as unorthogonalized NMF or CPG columns) often increase the number of active-set or Lemke LCP iterations $k_{\text{ROM}}$ required for convergence relative to the full-order iteration count $k_{\text{FOM}}$ \cite{leberre2022condition}.

To quantify the exact computational crossover, we model the online floating-point operation (FLOP) counts of both solvers for a single parameter query $\mu \in \mathcal{P}$:

#### Full-Order Online Cost ($\mathcal{C}_{\text{FOM}}$)
Assuming a sparse direct Cholesky factorization of $K$ (costing $c_{\text{fact}} N b_w^2$ FLOPs, or $c_{\text{fact}} N^\alpha$ in 3D) followed by $n_c$ forward/backward substitutions to form the dense contact compliance matrix $B K^{-1} B^T \in \mathbb{R}^{n_c \times n_c}$ and $k_{\text{FOM}}$ active-set updates of average active-set size $m_{\text{act}} \le n_c$:
$$\mathcal{C}_{\text{FOM}} = \underbrace{c_{\text{fact}} N b_w^2}_{\text{sparse } K \text{ factorization}} + \underbrace{2 c_{\text{solve}} n_c N b_w}_{\text{form } B K^{-1} B^T} + \underbrace{k_{\text{FOM}} \left( \frac{2}{3} m_{\text{act}}^3 + 2 n_c m_{\text{act}} \right)}_{\text{full-order active-set / LCP solves}}. \tag{4}$$
*(Note: When $K$ is parameter-independent, the first two terms are amortized offline, leaving only the active-set LCP cost $\mathcal{C}_{\text{FOM}}^{\text{affine}} = k_{\text{FOM}} (\frac{2}{3} m_{\text{act}}^3 + 2 n_c m_{\text{act}})$; when $K(\mu)$ varies with $\mu$, all three terms execute online.)*

#### Reduced-Order Online Cost ($\mathcal{C}_{\text{ROM}}$)
Let $r_u = r_{u,0} + r_{\text{rep}}$ be the enriched primal dimension and $n_{c,\text{eff}} \le n_c$ be the number of online contact rows evaluated ($n_{c,\text{eff}} = n_c$ for global dual bases, $n_{c,\text{eff}} = r_\lambda$ for gIS without global certification, or $n_{c,\text{eff}} = n_{c,\text{RID}}$ for RID). Forming the dense reduced operators $K_r(\mu)$ and $B_r(\mu)$ (in the general non-affine case, or assembling from affine components), factoring dense $K_r \in \mathbb{R}^{r_u \times r_u}$, computing the dense Schur complement $B_r K_r^{-1} B_r^T \in \mathbb{R}^{r_\lambda \times r_\lambda}$, running $k_{\text{ROM}}$ active-set iterations of average reduced active size $m_{r,\text{act}} \le r_\lambda$, and evaluating the online admissibility indicator $\eta_{\text{adm}}$ require:
$$\mathcal{C}_{\text{ROM}} = \underbrace{\mathcal{C}_{\text{proj}}(N_{\text{eff}}, r_u, r_\lambda)}_{\text{reduced operator assembly}} + \underbrace{\frac{1}{3} r_u^3 + 2 r_\lambda r_u^2 + 2 r_\lambda^2 r_u}_{\text{dense } K_r^{-1} \text{ \& } B_r K_r^{-1} B_r^T} + \underbrace{k_{\text{ROM}} \left( \frac{2}{3} m_{r,\text{act}}^3 + 2 r_\lambda m_{r,\text{act}} \right)}_{\text{reduced active-set / LCP solves}} + \underbrace{2 n_c (r_u + r_\lambda)}_{\text{online } \eta_{\text{adm}} \text{ check}}, \tag{5}$$
where $\mathcal{C}_{\text{proj}} = \mathcal{O}(Q_K r_u^2 + Q_B r_u r_\lambda)$ for affine problems with $Q_K, Q_B$ affine terms, and $\mathcal{C}_{\text{proj}} = \mathcal{O}(\text{nnz}(K) r_u + \text{nnz}(B) (r_u + r_\lambda))$ when full-order matrices are projected online.

#### Break-Even Conditions
Defining the theoretical speedup ratio $S := \mathcal{C}_{\text{FOM}} / \mathcal{C}_{\text{ROM}}$, the reduced contact model achieves computational break-even ($S > 1$) if and only if $\mathcal{C}_{\text{ROM}} < \mathcal{C}_{\text{FOM}}$. Equations (4)–(5) reveal two distinct break-even regimes tested in Section 5:

1. **Parameter-Varying Stiffness Regime ($K = K(\mu)$):**
   Because $\mathcal{C}_{\text{FOM}}$ is dominated by the $\mathcal{O}(N b_w^2)$ sparse factorization and $\mathcal{O}(n_c N b_w)$ compliance assembly, break-even holds whenever
   $$r_u^3 + 6 r_\lambda r_u^2 + 2 k_{\text{ROM}} m_{r,\text{act}}^3 \ll 3 c_{\text{fact}} N b_w^2 + 6 c_{\text{solve}} n_c N b_w.$$
   In this regime, substantial speedups ($S \sim 10^1\text{–}10^3$) persist even with moderate rank-repair enrichment $r_{\text{rep}}$, provided $r_u \ll (N b_w^2)^{1/3}$.
2. **Parameter-Independent Stiffness Regime ($K$ constant, varying loads/gaps $f(\mu), c(\mu)$):**
   When $B K^{-1} B^T$ is precomputed offline for the full-order solver, $\mathcal{C}_{\text{FOM}}^{\text{affine}}$ depends solely on the active contact count $m_{\text{act}}$. Break-even then reduces to the sharp interface-compression condition
   $$r_u^3 + 6 r_\lambda r_u^2 + 2 k_{\text{ROM}} m_{r,\text{act}}^3 + 6 n_c (r_u + r_\lambda) < 2 k_{\text{FOM}} m_{\text{act}}^3.$$
   Consequently, if the physical contact patch $m_{\text{act}}$ is small (e.g., localized point or line contact with $m_{\text{act}} \lesssim 30$) and rank repair requires adding $r_{\text{rep}} \sim 15\text{–}30$ unit vectors (such as $r_{\text{rep}} = 14$ on the 60-spring chain), the cubic dense cost $r_u^3$ immediately exceeds $2 k_{\text{FOM}} m_{\text{act}}^3$, causing $S < 1$. Thus, **Lagrange-multiplier contact ROMs are computationally advantageous in the constant-$K$ regime only when the active contact zone $m_{\text{act}}$ is large relative to the enriched reduced dimension $r_u = r_{u,0} + r_{\text{rep}}$.**
