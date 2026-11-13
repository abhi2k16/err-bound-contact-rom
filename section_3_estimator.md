# Section 3: Estimator — All Versions

---

## v2 (current)

# 3. Computable *A Posteriori* Lifted Admissibility Estimator

To detect silent reduced-KKT convergence without solving the full-order contact problem, we derive a residual-based *a posteriori* error bound for arbitrary lifted primal--dual pairs $(\tilde{u}, \tilde{\lambda}) = (V u_r, U \lambda_r)$. Unlike existing reduced-basis estimators for variational inequalities that require a bijective constraint operator and strictly feasible dual spaces \cite{zhang2014duality}, our analysis applies to general rectangular constraint matrices $B(\mu) \in \mathbb{R}^{n_c \times N}$ and makes no assumptions on how $(u_r, \lambda_r)$ is computed or whether $U$ preserves non-negativity \cite{banz2024generic, banz2026posteriori}.

## 3.1 Norms, Multiplier Splitting, and Estimator Components
For a fixed parameter $\mu \in \mathcal{P}$, suppress the explicit $\mu$-dependence for brevity and write $K = K(\mu)$, $B = B(\mu)$, $c = c(\mu)$, and $f = f(\mu)$. Since $K \in \mathbb{R}^{N \times N}$ is SPD, it induces the standard energy norm and dual energy norm \cite{hild2007residual, zhang2014duality}
$$
\|v\|_K := \left( v^T K v \right)^{1/2}, \qquad \|r\|_{K^{-1}} := \left( r^T K^{-1} r \right)^{1/2} = \sup_{v \in \mathbb{R}^N \setminus \{0\}} \frac{r^T v}{\|v\|_K}.
$$
Similarly, define the discrete trace operator norm of $B$ with respect to the energy norm \cite{coorevits2000posteriori, hild2007residual}:
$$
\|B\|_{K \to 2} := \sup_{v \in \mathbb{R}^N \setminus \{0\}} \frac{\|B v\|_2}{\|v\|_K} = \left\| B K^{-1} B^T \right\|_2^{1/2}.
$$
For any vector $z \in \mathbb{R}^{n_c}$, denote its componentwise positive and negative parts by $z_+ := \max(z, 0) \ge 0$ and $z_- := \min(z, 0) \le 0$, so that $z = z_+ + z_-$ and $z_+^T z_- = 0$.

Following the constant-sign multiplier decomposition and non-Galerkin residual splitting introduced for full-order mixed contact and variational inequalities \cite{hild2007residual, dabaghi2020posteriori, banz2024generic, banz2026posteriori}, we split the lifted Lagrange multiplier $\tilde{\lambda} = U \lambda_r$ into its physically admissible (compressive) and inadmissible (adhesive) parts:
$$
\tilde{\lambda} = \lambda_+ + \lambda_-, \qquad \lambda_+ := \max(\tilde{\lambda}, 0), \qquad \lambda_- := \min(\tilde{\lambda}, 0).
$$
Using $\lambda_+$ as a feasible dual proxy \cite{hild2007residual, fontana2024guaranteed}, we define four computable residual indicators for the lifted pair $(\tilde{u}, \tilde{\lambda})$:
$$
\begin{aligned}
\eta_{\mathrm{eq}} &:= \left\| K \tilde{u} + B^T \lambda_+ - f \right\|_{K^{-1}}, & \eta_{\mathrm{sign}} &:= \left\| \lambda_- \right\|_2, \\
\eta_{\mathrm{gap}} &:= \left\| (B \tilde{u} - c)_+ \right\|_2, & \eta_{\mathrm{comp}} &:= \left| \lambda_+^T (B \tilde{u} - c) \right|.
\end{aligned}
$$
Each component isolates a distinct physical or algebraic defect of the lifted solution \cite{dabaghi2020posteriori, banz2024generic, fontana2024guaranteed}:
* $\eta_{\mathrm{eq}}$ measures the bulk force-balance residual when only the non-negative part $\lambda_+$ of the lifted traction is applied \cite{hild2007residual};
* $\eta_{\mathrm{sign}}$ measures the Euclidean norm of non-physical adhesive contact tractions ($\tilde{\lambda}_i < 0$) \cite{balajewicz2016projection, banz2026posteriori};
* $\eta_{\mathrm{gap}}$ measures the Euclidean norm of positive interpenetration across all $n_c$ candidate contact constraints \cite{coorevits2000posteriori, hild2007residual};
* $\eta_{\mathrm{comp}}$ measures the violation of complementarity slackness between the compressive tractions $\lambda_+$ and the physical gap vector $B \tilde{u} - c$ \cite{dabaghi2020posteriori, banz2024generic}.

## 3.2 Certified Energy-Norm Upper Bound
Adapting the primal--dual contact error decomposition of Hild and Nicaise \cite{hild2007residual} and Banz et al. \cite{banz2024generic, banz2026posteriori} to lifted reduced-order pairs without Galerkin orthogonality, the following theorem establishes that the energy-norm displacement error $\|u - \tilde{u}\|_K$ is rigorously bounded from above by these four quantities.

**Theorem 1 (Computable *A Posteriori* Error Bound for Lifted Contact Pairs).**
*Let $(u, \lambda) \in \mathbb{R}^N \times \mathbb{R}^{n_c}$ be the exact full-order KKT solution of the discrete contact problem, and let $(\tilde{u}, \tilde{\lambda}) = (V u_r, U \lambda_r)$ be any lifted reduced-order approximation. Then the energy-norm displacement error satisfies*
$$
\|u - \tilde{u}\|_K^2 \le \eta_{\mathrm{eq}}^2 + 2\,\eta_{\mathrm{comp}} + 2\,\|\lambda\|_2\,\eta_{\mathrm{gap}}. \tag{1}
$$
*Furthermore, letting $\tilde{r}_{\mathrm{eq}} := \|K \tilde{u} + B^T \tilde{\lambda} - f\|_{K^{-1}}$ denote the equilibrium residual evaluated with the unprojected lifted multiplier $\tilde{\lambda}$, the error bound decomposes explicitly into all four indicators:*
$$
\|u - \tilde{u}\|_K^2 \le 2\,\tilde{r}_{\mathrm{eq}}^2 + 2\,\|B\|_{K \to 2}^2\,\eta_{\mathrm{sign}}^2 + 2\,\eta_{\mathrm{comp}} + 2\,\|\lambda\|_2\,\eta_{\mathrm{gap}}. \tag{2}
$$

*Proof.* Let $e := u - \tilde{u} \in \mathbb{R}^N$. Using the exact equilibrium relation $K u = f - B^T \lambda$ and adding and subtracting $B^T \lambda_+$, we expand the squared energy-norm error as \cite{hild2007residual, banz2024generic}
$$
\begin{aligned}
\|e\|_K^2 &= e^T K (u - \tilde{u}) \\
&= e^T \left( f - B^T \lambda - K \tilde{u} \right) \\
&= e^T \left( f - K \tilde{u} - B^T \lambda_+ \right) + (B e)^T (\lambda_+ - \lambda).
\end{aligned}
$$
For the first term, applying the Cauchy--Schwarz inequality in the $K$-inner product followed by Young's inequality ($a b \le \frac{1}{2} a^2 + \frac{1}{2} b^2$) yields
$$
e^T \left( f - K \tilde{u} - B^T \lambda_+ \right) \le \|e\|_K \, \left\| K \tilde{u} + B^T \lambda_+ - f \right\|_{K^{-1}} = \|e\|_K \, \eta_{\mathrm{eq}} \le \frac{1}{2} \|e\|_K^2 + \frac{1}{2} \eta_{\mathrm{eq}}^2.
$$
For the second term, rewrite $B e = (B u - c) - (B \tilde{u} - c)$ and expand the inner product into four contributions \cite{hild2007residual, banz2026posteriori}:
$$
(B e)^T (\lambda_+ - \lambda) = (B u - c)^T \lambda_+ - (B u - c)^T \lambda - (B \tilde{u} - c)^T \lambda_+ + (B \tilde{u} - c)^T \lambda.
$$
By full-order primal feasibility ($B u - c \le 0$) and non-negativity of $\lambda_+$ ($\lambda_+ \ge 0$), the first inner product is non-positive: $(B u - c)^T \lambda_+ \le 0$. By exact full-order complementarity, the second inner product vanishes identically: $(B u - c)^T \lambda = 0$. Dropping these two non-positive terms gives
$$
(B e)^T (\lambda_+ - \lambda) \le -(B \tilde{u} - c)^T \lambda_+ + (B \tilde{u} - c)^T \lambda.
$$
The third term is bounded immediately by $\left| \lambda_+^T (B \tilde{u} - c) \right| = \eta_{\mathrm{comp}}$. For the fourth term, decompose $B \tilde{u} - c = (B \tilde{u} - c)_+ + (B \tilde{u} - c)_-$. Since $\lambda \ge 0$ and $(B \tilde{u} - c)_- \le 0$, we have $(B \tilde{u} - c)_-^T \lambda \le 0$, and hence by the Euclidean Cauchy--Schwarz inequality,
$$
(B \tilde{u} - c)^T \lambda \le (B \tilde{u} - c)_+^T \lambda \le \|\lambda\|_2 \, \left\| (B \tilde{u} - c)_+ \right\|_2 = \|\lambda\|_2 \, \eta_{\mathrm{gap}}.
$$
Combining these bounds into the expansion of $\|e\|_K^2$ gives
$$
\|e\|_K^2 \le \frac{1}{2} \|e\|_K^2 + \frac{1}{2} \eta_{\mathrm{eq}}^2 + \eta_{\mathrm{comp}} + \|\lambda\|_2 \, \eta_{\mathrm{gap}}.
$$
Subtracting $\frac{1}{2} \|e\|_K^2$ from both sides and multiplying by $2$ establishes (1).

To obtain (2), substitute $\lambda_+ = \tilde{\lambda} - \lambda_-$ into $\eta_{\mathrm{eq}}$ and apply the triangle inequality together with $(a + b)^2 \le 2 a^2 + 2 b^2$:
$$
\eta_{\mathrm{eq}}^2 = \left\| \left( K \tilde{u} + B^T \tilde{\lambda} - f \right) - B^T \lambda_- \right\|_{K^{-1}}^2 \le 2\,\tilde{r}_{\mathrm{eq}}^2 + 2\,\left\| B^T \lambda_- \right\|_{K^{-1}}^2.
$$
By definition of the dual energy norm and the trace norm $\|B\|_{K \to 2}$ \cite{coorevits2000posteriori, franceschini2022block},
$$
\left\| B^T \lambda_- \right\|_{K^{-1}} = \sup_{v \in \mathbb{R}^N \setminus \{0\}} \frac{\lambda_-^T B v}{\|v\|_K} \le \sup_{v \in \mathbb{R}^N \setminus \{0\}} \frac{\|\lambda_-\|_2 \, \|B v\|_2}{\|v\|_K} = \|B\|_{K \to 2} \, \eta_{\mathrm{sign}}.
$$
Inserting this estimate into (1) yields (2). $\blacksquare$

**Remark 1 (Practical Evaluation of Constants and $\|\lambda\|_2$).**
In (1)--(2), the exact multiplier norm $\|\lambda\|_2$ is unknown online. When $B$ has full row rank with full-order discrete inf-sup constant $\beta > 0$ (i.e., $\inf_{\mu_c \neq 0} \sup_{v \neq 0} \mu_c^T B v / (\|\mu_c\|_2 \|v\|_K) \ge \beta$), standard saddle-point stability \cite{hild2007residual, franceschini2022block} bounds $\|\lambda\|_2 \le \beta^{-1} (\|f\|_{K^{-1}} + \|K u\|_K) \le 2 \beta^{-1} \|f\|_{K^{-1}}$. In practical online diagnostics, replacing $\|\lambda\|_2$ with the computable proxy $\|\lambda_+\|_2$ (or a precomputed offline upper bound $C_\lambda = \max_{\mu \in \mathcal{P}_{\mathrm{train}}} \|\lambda(\mu)\|_2$) preserves sharp scaling across parameter regimes \cite{fontana2024guaranteed}.

**Remark 2 (Why Signed SVD Fails Silently).**
Bound (2) directly explains the empirical breakdown documented in Section 2.3. Even if a reduced solver drives the projected residual $\tilde{r}_{\mathrm{eq}}$ to a small value, a signed SVD dual basis permits large negative lifted multipliers $\lambda_- \neq 0$. The resulting sign violation $\eta_{\mathrm{sign}} = \|\lambda_-\|_2$ enters the upper bound multiplied by $2 \|B\|_{K \to 2}^2$, simultaneously inducing positive penetration $\eta_{\mathrm{gap}}$ and complementarity error $\eta_{\mathrm{comp}}$. Conversely, when a strictly non-negative dual basis is used (such as NMF \cite{balajewicz2016projection}, cone-projected greedy \cite{benaceur2020reduced}, or greedy index selection \cite{sambataro2026model}), $\eta_{\mathrm{sign}} = 0$ identically, collapsing (2) to $\|u - \tilde{u}\|_K^2 \le \eta_{\mathrm{eq}}^2 + 2\,\eta_{\mathrm{comp}} + 2\,\|\lambda\|_2\,\eta_{\mathrm{gap}}$.

## 3.3 Computational Complexity and Online Admissibility Indicator
Evaluating the four estimator components splits naturally into a bulk equilibrium term and three contact-interface admissibility terms:

| Estimator Term | Definition | Online Operation Count | Bulk $N$-Dependence |
| :--- | :--- | :--- | :--- |
| $\eta_{\mathrm{sign}}$ | $\|\min(U \lambda_r, 0)\|_2$ | $\mathcal{O}(n_c r_\lambda)$ | None |
| $\eta_{\mathrm{gap}}$ | $\|\max((B V) u_r - c, 0)\|_2$ | $\mathcal{O}(n_c r_u)$ | None (if $B V$ precomputed) |
| $\eta_{\mathrm{comp}}$ | $|\lambda_+^T ((B V) u_r - c)|$ | $\mathcal{O}(n_c)$ | None |
| $\eta_{\mathrm{eq}}$ | $\|K \tilde{u} + B^T \lambda_+ - f\|_{K^{-1}}$ | $\mathcal{O}(N^\alpha)$ (general) / $\mathcal{O}((r_u + r_\lambda)^2)$ (affine) | Full-order solve unless affine |

When the parameter dependence of $B(\mu)$ and $c(\mu)$ is affine (or parameter-independent, as in fixed-interface contact with varying external loads), the matrix products $B V \in \mathbb{R}^{n_c \times r_u}$ are precomputed once offline \cite{benaceur2020reduced, zhang2014duality}. Online evaluation of $\eta_{\mathrm{sign}}$, $\eta_{\mathrm{gap}}$, and $\eta_{\mathrm{comp}}$ then requires only $\mathcal{O}(n_c (r_u + r_\lambda))$ floating-point operations---completely independent of the bulk mesh dimension $N$.

For $\eta_{\mathrm{eq}}$, computing the dual norm $\|r\|_{K^{-1}}$ in general non-affine problems requires a full-order linear solve with $K$, making $\eta_{\mathrm{eq}}$ an offline validation metric or greedy snapshot-selection criterion rather than an online diagnostic. In strictly affine settings, $\eta_{\mathrm{eq}}^2$ can be evaluated online in $\mathcal{O}((r_u + r_\lambda)^2)$ operations via standard offline--online Riesz representation of the residual components \cite{zhang2014duality, benaceur2020reduced}, though we reserve $\eta_{\mathrm{eq}}$ for validation in our general framework.

To provide a fast, dimensionally consistent online "circuit breaker" that operates solely on the contact interface \cite{dabaghi2020posteriori, fauque2018hybrid}, we combine the three $\mathcal{O}(n_c (r_u + r_\lambda))$ terms into the **Online Admissibility Indicator**:
$$
\eta_{\mathrm{adm}}(\tilde{u}, \tilde{\lambda}) := \eta_{\mathrm{sign}} + \eta_{\mathrm{gap}} + \sqrt{\eta_{\mathrm{comp}}}. \tag{3}
$$
Taking the square root of $\eta_{\mathrm{comp}}$ ensures that all three terms scale linearly with the first power of the solution error $\mathcal{O}(\|e\|_K)$, matching the scaling of $\|u - \tilde{u}\|_K$ in Theorem 1. Whenever $\eta_{\mathrm{adm}}(\tilde{u}, \tilde{\lambda})$ exceeds a user-specified admissibility tolerance $\tau_{\mathrm{adm}}$, the reduced prediction is flagged as physically inadmissible---even if the reduced KKT residual is at machine precision---triggering either adaptive basis enrichment or fallback to the full-order solver.

---

## v1 (original)

### 3. Computable *A Posteriori* Lifted Admissibility Estimator

To detect silent reduced-KKT convergence without solving the full-order contact problem, we derive a residual-based *a posteriori* error bound for arbitrary lifted primal–dual pairs $(\tilde{u}, \tilde{\lambda}) = (V u_r, U \lambda_r)$. Unlike existing reduced-basis estimators for variational inequalities that require a bijective constraint operator and strictly feasible dual spaces \cite{zhang2014duality}, our analysis applies to general rectangular constraint matrices $B(\mu) \in \mathbb{R}^{n_c \times N}$ and makes no assumptions on how $(u_r, \lambda_r)$ is computed or whether $U$ preserves non-negativity.

#### 3.1 Norms, Multiplier Splitting, and Estimator Components
For a fixed parameter $\mu \in \mathcal{P}$, suppress the explicit $\mu$-dependence for brevity and write $K = K(\mu)$, $B = B(\mu)$, $c = c(\mu)$, and $f = f(\mu)$. Since $K \in \mathbb{R}^{N \times N}$ is SPD, it induces the energy norm and dual energy norm
$$\|v\|_K := \left( v^T K v \right)^{1/2}, \qquad \|r\|_{K^{-1}} := \left( r^T K^{-1} r \right)^{1/2} = \sup_{v \in \mathbb{R}^N \setminus \{0\}} \frac{r^T v}{\|v\|_K}.$$
Similarly, define the discrete trace operator norm of $B$ with respect to the energy norm:
$$\|B\|_{K \to 2} := \sup_{v \in \mathbb{R}^N \setminus \{0\}} \frac{\|B v\|_2}{\|v\|_K} = \left\| B K^{-1} B^T \right\|_2^{1/2}.$$
For any vector $z \in \mathbb{R}^{n_c}$, denote its componentwise positive and negative parts by $z_+ := \max(z, 0) \ge 0$ and $z_- := \min(z, 0) \le 0$, so that $z = z_+ + z_-$ and $z_+^T z_- = 0$.

Following the constant-sign multiplier decomposition introduced for full-order mixed contact problems \cite{hild2007residual, banz2024generic, banz2026posteriori}, we split the lifted Lagrange multiplier $\tilde{\lambda} = U \lambda_r$ into its physically admissible (compressive) and inadmissible (adhesive) parts:
$$\tilde{\lambda} = \lambda_+ + \lambda_-, \qquad \lambda_+ := \max(\tilde{\lambda}, 0), \qquad \lambda_- := \min(\tilde{\lambda}, 0).$$
Using $\lambda_+$ as a feasible dual proxy, we define four computable residual indicators for the lifted pair $(\tilde{u}, \tilde{\lambda})$:
$$\begin{aligned} \eta_{\text{eq}} &:= \left\| K \tilde{u} + B^T \lambda_+ - f \right\|_{K^{-1}}, & \eta_{\text{sign}} &:= \left\| \lambda_- \right\|_2, \\ \eta_{\text{gap}} &:= \left\| (B \tilde{u} - c)_+ \right\|_2, & \eta_{\text{comp}} &:= \left| \lambda_+^T (B \tilde{u} - c) \right|. \end{aligned}$$
Each component isolates a distinct physical or algebraic defect of the lifted solution:
* $\eta_{\text{eq}}$ measures the bulk force-balance residual when only the non-negative part $\lambda_+$ of the lifted traction is applied;
* $\eta_{\text{sign}}$ measures the Euclidean norm of non-physical adhesive contact tractions ($\tilde{\lambda}_i < 0$);
* $\eta_{\text{gap}}$ measures the Euclidean norm of positive interpenetration across all $n_c$ candidate contact constraints;
* $\eta_{\text{comp}}$ measures the violation of complementarity slackness between the compressive tractions $\lambda_+$ and the physical gap vector $B \tilde{u} - c$.

#### 3.2 Certified Energy-Norm Upper Bound
The following theorem establishes that the energy-norm displacement error $\|u - \tilde{u}\|_K$ is rigorously bounded from above by these four quantities.

**Theorem 1 (Computable *A Posteriori* Error Bound for Lifted Contact Pairs).**
*Let $(u, \lambda) \in \mathbb{R}^N \times \mathbb{R}^{n_c}$ be the exact full-order KKT solution of the discrete contact problem, and let $(\tilde{u}, \tilde{\lambda}) = (V u_r, U \lambda_r)$ be any lifted reduced-order approximation. Then the energy-norm displacement error satisfies*
$$\|u - \tilde{u}\|_K^2 \le \eta_{\text{eq}}^2 + 2\,\eta_{\text{comp}} + 2\,\|\lambda\|_2\,\eta_{\text{gap}}. \tag{1}$$
*Furthermore, letting $\tilde{r}_{\text{eq}} := \|K \tilde{u} + B^T \tilde{\lambda} - f\|_{K^{-1}}$ denote the equilibrium residual evaluated with the unprojected lifted multiplier $\tilde{\lambda}$, the error bound decomposes explicitly into all four indicators:*
$$\|u - \tilde{u}\|_K^2 \le 2\,\tilde{r}_{\text{eq}}^2 + 2\,\|B\|_{K \to 2}^2\,\eta_{\text{sign}}^2 + 2\,\eta_{\text{comp}} + 2\,\|\lambda\|_2\,\eta_{\text{gap}}. \tag{2}$$

*Proof.* Let $e := u - \tilde{u} \in \mathbb{R}^N$. Using the exact equilibrium relation $K u = f - B^T \lambda$ and adding and subtracting $B^T \lambda_+$, we expand the squared energy-norm error as
$$\begin{aligned} \|e\|_K^2 &= e^T K (u - \tilde{u}) \\ &= e^T \left( f - B^T \lambda - K \tilde{u} \right) \\ &= e^T \left( f - K \tilde{u} - B^T \lambda_+ \right) + (B e)^T (\lambda_+ - \lambda). \end{aligned}$$
For the first term, applying the Cauchy–Schwarz inequality in the $K$-inner product followed by Young's inequality ($a b \le \frac{1}{2} a^2 + \frac{1}{2} b^2$) yields
$$e^T \left( f - K \tilde{u} - B^T \lambda_+ \right) \le \|e\|_K \, \left\| K \tilde{u} + B^T \lambda_+ - f \right\|_{K^{-1}} = \|e\|_K \, \eta_{\text{eq}} \le \frac{1}{2} \|e\|_K^2 + \frac{1}{2} \eta_{\text{eq}}^2.$$
For the second term, rewrite $B e = (B u - c) - (B \tilde{u} - c)$ and expand the inner product into four contributions:
$$(B e)^T (\lambda_+ - \lambda) = (B u - c)^T \lambda_+ - (B u - c)^T \lambda - (B \tilde{u} - c)^T \lambda_+ + (B \tilde{u} - c)^T \lambda.$$
By full-order primal feasibility ($B u - c \le 0$) and non-negativity of $\lambda_+$ ($\lambda_+ \ge 0$), the first inner product is non-positive: $(B u - c)^T \lambda_+ \le 0$. By exact full-order complementarity, the second inner product vanishes identically: $(B u - c)^T \lambda = 0$. Dropping these two non-positive terms gives
$$(B e)^T (\lambda_+ - \lambda) \le -(B \tilde{u} - c)^T \lambda_+ + (B \tilde{u} - c)^T \lambda.$$
The third term is bounded immediately by $\left| \lambda_+^T (B \tilde{u} - c) \right| = \eta_{\text{comp}}$. For the fourth term, decompose $B \tilde{u} - c = (B \tilde{u} - c)_+ + (B \tilde{u} - c)_-$. Since $\lambda \ge 0$ and $(B \tilde{u} - c)_- \le 0$, we have $(B \tilde{u} - c)_-^T \lambda \le 0$, and hence by the Euclidean Cauchy–Schwarz inequality,
$$(B \tilde{u} - c)^T \lambda \le (B \tilde{u} - c)_+^T \lambda \le \|\lambda\|_2 \, \left\| (B \tilde{u} - c)_+ \right\|_2 = \|\lambda\|_2 \, \eta_{\text{gap}}.$$
Combining these bounds into the expansion of $\|e\|_K^2$ gives
$$\|e\|_K^2 \le \frac{1}{2} \|e\|_K^2 + \frac{1}{2} \eta_{\text{eq}}^2 + \eta_{\text{comp}} + \|\lambda\|_2 \, \eta_{\text{gap}}.$$
Subtracting $\frac{1}{2} \|e\|_K^2$ from both sides and multiplying by $2$ establishes (1).

To obtain (2), substitute $\lambda_+ = \tilde{\lambda} - \lambda_-$ into $\eta_{\text{eq}}$ and apply the triangle inequality together with $(a + b)^2 \le 2 a^2 + 2 b^2$:
$$\eta_{\text{eq}}^2 = \left\| \left( K \tilde{u} + B^T \tilde{\lambda} - f \right) - B^T \lambda_- \right\|_{K^{-1}}^2 \le 2\,\tilde{r}_{\text{eq}}^2 + 2\,\left\| B^T \lambda_- \right\|_{K^{-1}}^2.$$
By definition of the dual energy norm and the trace norm $\|B\|_{K \to 2}$,
$$\left\| B^T \lambda_- \right\|_{K^{-1}} = \sup_{v \in \mathbb{R}^N \setminus \{0\}} \frac{\lambda_-^T B v}{\|v\|_K} \le \sup_{v \in \mathbb{R}^N \setminus \{0\}} \frac{\|\lambda_-\|_2 \, \|B v\|_2}{\|v\|_K} = \|B\|_{K \to 2} \, \eta_{\text{sign}}.$$
Inserting this estimate into (1) yields (2). $\blacksquare$

**Remark 1 (Practical Evaluation of Constants and $\|\lambda\|_2$).**
In (1)–(2), the exact multiplier norm $\|\lambda\|_2$ is unknown online. When $B$ has full row rank with full-order discrete inf-sup constant $\beta > 0$ (i.e., $\inf_{\mu_c \neq 0} \sup_{v \neq 0} \mu_c^T B v / (\|\mu_c\|_2 \|v\|_K) \ge \beta$), standard saddle-point stability bounds $\|\lambda\|_2 \le \beta^{-1} (\|f\|_{K^{-1}} + \|K u\|_K) \le 2 \beta^{-1} \|f\|_{K^{-1}}$. In practical online diagnostics, replacing $\|\lambda\|_2$ with the computable proxy $\|\lambda_+\|_2$ (or a precomputed offline upper bound $C_\lambda = \max_{\mu \in \mathcal{P}_{\text{train}}} \|\lambda(\mu)\|_2$) preserves sharp scaling across parameter regimes.

**Remark 2 (Why Signed SVD Fails Silently).**
Bound (2) directly explains the empirical breakdown documented in Section 2.3. Even if a reduced solver drives the projected residual $\tilde{r}_{\text{eq}}$ to a small value, a signed SVD dual basis permits large negative lifted multipliers $\lambda_- \neq 0$. The resulting sign violation $\eta_{\text{sign}} = \|\lambda_-\|_2$ enters the upper bound multiplied by $2 \|B\|_{K \to 2}^2$, simultaneously inducing positive penetration $\eta_{\text{gap}}$ and complementarity error $\eta_{\text{comp}}$. Conversely, when a strictly non-negative dual basis is used (such as NMF, cone-projected greedy, or greedy index selection), $\eta_{\text{sign}} = 0$ identically, collapsing (2) to $\|u - \tilde{u}\|_K^2 \le \eta_{\text{eq}}^2 + 2\,\eta_{\text{comp}} + 2\,\|\lambda\|_2\,\eta_{\text{gap}}$.

#### 3.3 Computational Complexity and Online Admissibility Indicator
Evaluating the four estimator components splits naturally into a bulk equilibrium term and three contact-interface admissibility terms:

| Estimator Term | Definition | Online Operation Count | Bulk $N$-Dependence |
| :--- | :--- | :--- | :--- |
| $\eta_{\text{sign}}$ | $\|\min(U \lambda_r, 0)\|_2$ | $\mathcal{O}(n_c r_\lambda)$ | None |
| $\eta_{\text{gap}}$ | $\|\max((B V) u_r - c, 0)\|_2$ | $\mathcal{O}(n_c r_u)$ | None (if $B V$ precomputed) |
| $\eta_{\text{comp}}$ | $|\lambda_+^T ((B V) u_r - c)|$ | $\mathcal{O}(n_c)$ | None |
| $\eta_{\text{eq}}$ | $\|K \tilde{u} + B^T \lambda_+ - f\|_{K^{-1}}$ | $\mathcal{O}(N^\alpha)$ (general) / $\mathcal{O}((r_u + r_\lambda)^2)$ (affine) | Full-order solve unless affine |

When the parameter dependence of $B(\mu)$ and $c(\mu)$ is affine (or parameter-independent, as in fixed-interface contact with varying external loads), the matrix products $B V \in \mathbb{R}^{n_c \times r_u}$ are precomputed once offline. Online evaluation of $\eta_{\text{sign}}$, $\eta_{\text{gap}}$, and $\eta_{\text{comp}}$ then requires only $\mathcal{O}(n_c (r_u + r_\lambda))$ floating-point operations—completely independent of the bulk mesh dimension $N$.

For $\eta_{\text{eq}}$, computing the dual norm $\|r\|_{K^{-1}}$ in general non-affine problems requires a full-order linear solve with $K$, making $\eta_{\text{eq}}$ an offline validation metric or greedy snapshot-selection criterion rather than an online diagnostic. (In strictly affine settings, $\eta_{\text{eq}}^2$ can be evaluated online in $\mathcal{O}((r_u + r_\lambda)^2)$ operations via standard offline–online Riesz representation of the residual components, though we reserve $\eta_{\text{eq}}$ for validation in our general framework.)

To provide a fast, dimensionally consistent online "circuit breaker" that operates solely on the contact interface, we combine the three $\mathcal{O}(n_c (r_u + r_\lambda))$ terms into the **Online Admissibility Indicator**:
$$\eta_{\text{adm}}(\tilde{u}, \tilde{\lambda}) := \eta_{\text{sign}} + \eta_{\text{gap}} + \sqrt{\eta_{\text{comp}}}. \tag{3}$$
Taking the square root of $\eta_{\text{comp}}$ ensures that all three terms scale linearly with the first power of the solution error $\mathcal{O}(\|e\|_K)$, matching the scaling of $\|u - \tilde{u}\|_K$ in Theorem 1. Whenever $\eta_{\text{adm}}(\tilde{u}, \tilde{\lambda})$ exceeds a user-specified admissibility tolerance $\tau_{\text{adm}}$, the reduced prediction is flagged as physically inadmissible—even if the reduced KKT residual is at machine precision—triggering either adaptive basis enrichment or fallback to the full-order solver.
