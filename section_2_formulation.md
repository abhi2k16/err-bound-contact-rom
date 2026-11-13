# Section 2: Formulation — All Versions

---

## v2 (current)

### 2. Problem Formulation and Silent Reduced-KKT Failure

#### 2.1 Full-Order Discrete Unilateral Contact Problem
Consider a parametrized linear-elastic or quasi-static discrete mechanical system with $N$ displacement degrees of freedom and $n_c$ unilateral contact constraints ($n_c \le N$) \cite{coorevits2000posteriori, hild2007residual}. For a parameter vector $\mu \in \mathcal{P} \subset \mathbb{R}^{n_p}$, the full-order equilibrium and non-penetration conditions are governed by the constrained quadratic minimization problem \cite{balajewicz2016projection, benaceur2020reduced}:
$$\min_{u \in \mathcal{K}(\mu)} \frac{1}{2} u^T K(\mu) u - f(\mu)^T u, \qquad \mathcal{K}(\mu) := \left\{ v \in \mathbb{R}^N : B(\mu) v - c(\mu) \le 0 \right\},$$
where $K(\mu) \in \mathbb{R}^{N \times N}$ is the symmetric positive definite (SPD) stiffness matrix, $f(\mu) \in \mathbb{R}^N$ is the external load vector, $B(\mu) \in \mathbb{R}^{n_c \times N}$ is a rectangular constraint Jacobian mapping nodal displacements to normal contact gaps (arising from node-to-segment or mortar contact discretizations \cite{coorevits2000posteriori, franceschini2022block}), and $c(\mu) \in \mathbb{R}^{n_c}$ is the initial gap vector. Throughout this work, inequalities on vectors are understood componentwise.

Introducing the vector of contact Lagrange multipliers $\lambda \in \mathbb{R}^{n_c}$ representing normal contact tractions, the first-order optimality conditions take the mixed Karush--Kuhn--Tucker (KKT) saddle-point form \cite{hild2007residual, banz2024generic}: find $(u(\mu), \lambda(\mu)) \in \mathbb{R}^N \times \mathbb{R}^{n_c}$ such that
$$\begin{aligned}
K(\mu) u(\mu) + B(\mu)^T \lambda(\mu) &= f(\mu), \\
\lambda(\mu) &\ge 0, \\
B(\mu) u(\mu) - c(\mu) &\le 0, \\
\lambda(\mu)^T \bigl(B(\mu) u(\mu) - c(\mu)\bigr) &= 0.
\end{aligned}$$
Equivalently, eliminating $u(\mu) = K(\mu)^{-1}(f(\mu) - B(\mu)^T \lambda(\mu))$ via Schur complement condensation yields the dense full-order linear complementarity problem (LCP) for the multipliers \cite{balajewicz2016projection, leberre2022condition}:
$$w(\mu) := B(\mu) K(\mu)^{-1} B(\mu)^T \lambda(\mu) + \bigl(c(\mu) - B(\mu) K(\mu)^{-1} f(\mu)\bigr) \ge 0, \qquad \lambda(\mu) \ge 0, \qquad \lambda(\mu)^T w(\mu) = 0.$$
Because $K(\mu)$ is SPD, a unique primal solution $u(\mu)$ exists for every $\mu \in \mathcal{P}$. Uniqueness of $\lambda(\mu)$ holds whenever the active rows of $B(\mu)$ have full row rank, which is satisfied on standard non-redundant contact interfaces \cite{franceschini2022block}.

#### 2.2 Galerkin--Petrov Mixed Reduced-Order Model
To construct a projection-based reduced-order model, we collect $n_s$ full-order snapshot pairs $\{(u(\mu_i), \lambda(\mu_i))\}_{i=1}^{n_s}$ over a training set $\mathcal{P}_{\text{train}} \subset \mathcal{P}$ \cite{balajewicz2016projection, benaceur2020reduced}. From these snapshots, we build a primal displacement basis $V \in \mathbb{R}^{N \times r_u}$ with orthonormal columns ($V^T V = I_{r_u}$, $r_u \ll N$) and a dual multiplier basis $U \in \mathbb{R}^{n_c \times r_\lambda}$ ($r_\lambda \le n_c$) \cite{balajewicz2016projection, benaceur2020reduced, fauque2018hybrid, sambataro2026model}.

Approximating the primal and dual variables in the reduced subspaces via
$$u(\mu) \approx \tilde{u}(\mu) := V u_r(\mu), \qquad \lambda(\mu) \approx \tilde{\lambda}(\mu) := U \lambda_r(\mu),$$
and testing the equilibrium and gap equations with $V^T$ and $U^T$, respectively, yields the reduced mixed KKT system \cite{balajewicz2016projection, leberre2022condition}: find $(u_r(\mu), \lambda_r(\mu)) \in \mathbb{R}^{r_u} \times \mathbb{R}^{r_\lambda}$ such that
$$\begin{aligned}
K_r(\mu) u_r(\mu) + B_r(\mu)^T \lambda_r(\mu) &= f_r(\mu), \\
\lambda_r(\mu) &\ge 0, \\
B_r(\mu) u_r(\mu) - c_r(\mu) &\le 0, \\
\lambda_r(\mu)^T \bigl(B_r(\mu) u_r(\mu) - c_r(\mu)\bigr) &= 0,
\end{aligned}$$
where the reduced operators are defined by
$$K_r(\mu) := V^T K(\mu) V \in \mathbb{R}^{r_u \times r_u}, \quad B_r(\mu) := U^T B(\mu) V \in \mathbb{R}^{r_\lambda \times r_u}, \quad f_r(\mu) := V^T f(\mu) \in \mathbb{R}^{r_u}, \quad c_r(\mu) := U^T c(\mu) \in \mathbb{R}^{r_\lambda}.$$
Because $V$ has full column rank and $K(\mu)$ is SPD, $K_r(\mu)$ is SPD, and the reduced system admits the Schur-complement LCP formulation \cite{balajewicz2016projection, leberre2022condition}:
$$w_r(\mu) := B_r(\mu) K_r(\mu)^{-1} B_r(\mu)^T \lambda_r(\mu) + \bigl(c_r(\mu) - B_r(\mu) K_r(\mu)^{-1} f_r(\mu)\bigr) \ge 0, \qquad \lambda_r(\mu) \ge 0, \qquad \lambda_r(\mu)^T w_r(\mu) = 0.$$

#### 2.3 Mechanisms of Silent Reduced-KKT Failure
A fundamental pitfall of the reduced formulation above is that enforcing the reduced KKT conditions to machine precision does **not** imply physical admissibility of the lifted pair $(\tilde{u}(\mu), \tilde{\lambda}(\mu))$. Specifically, define the reduced KKT residual vector:
$$R_{\text{red}}(u_r, \lambda_r) := \begin{bmatrix}
\Vert K_r u_r + B_r^T \lambda_r - f_r \Vert_2 \\
\Vert \min(\lambda_r, 0) \Vert_2 \\
\Vert \max(B_r u_r - c_r, 0) \Vert_2 \\
\vert \lambda_r^T (B_r u_r - c_r) \vert
\end{bmatrix}.$$
Even when $\Vert R_{\text{red}}(u_r, \lambda_r) \Vert_\infty \le 10^{-16}$, the lifted solution can fail catastrophically through three coupled algebraic mechanisms:

1. **Loss of Dual Cone Preservation (Signed Dual Bases):**
   When $U$ is obtained via standard truncated SVD or POD of multiplier snapshots, its columns contain mixed positive and negative entries \cite{balajewicz2016projection, benaceur2020reduced}. Consequently, $\lambda_r \ge 0$ does not guarantee $\tilde{\lambda} = U \lambda_r \ge 0$. Negative components in $\tilde{\lambda}$ act as artificial adhesive tractions that pull the elastic body across the contact boundary to satisfy reduced equilibrium $V^T(K \tilde{u} + B^T \tilde{\lambda} - f) = 0$.
2. **Projected Penetration Masking:**
   The reduced non-penetration condition $U^T (B \tilde{u} - c) \le 0$ constrains only $r_\lambda$ weighted averages of the $n_c$ physical gaps \cite{fauque2018hybrid, sambataro2026model}. When $r_\lambda < n_c$, positive penetrations $(B \tilde{u} - c)_i > 0$ at active contact nodes can cancel against large negative gaps $(B \tilde{u} - c)_j < 0$ at inactive nodes in the inner product $U^T(B\tilde{u} - c)$.
3. **Reduced Inf-Sup Rank Deficiency:**
   Even when $B(\mu)$ has full row rank on the active contact set, the projected matrix $B_r(\mu) = U^T B(\mu) V$ can become rank-deficient or severely ill-conditioned if $V$ lacks sufficient normal-displacement resolution on the contact interface \cite{benaceur2020reduced, leberre2022condition}. In that case, $B_r K_r^{-1} B_r^T$ is singular or nearly singular, allowing large spurious vectors in $\ker(B_r^T)$ to enter $\lambda_r$ with negligible change to the reduced residual \cite{leberre2022condition}.

##### Illustrative Benchmark Evidence
To make this failure concrete, consider the parametrized 2D obstacle membrane ($N = 324$) \cite{benaceur2020reduced} and 1D discrete spring-chain contact ($N = 60$) \cite{pawar2025efficient} benchmarks investigated in our preliminary measurements. When reduced with identical primal and dual dimensions $(r_u, r_\lambda)$ using a standard signed SVD dual basis, active-set and Lemke LCP solvers routinely converge to reduced KKT residuals below $10^{-16}$. Upon lifting, however, the signed SVD model exhibits lifted multipliers reaching $\min_i \tilde{\lambda}_i \approx -1.4 \times 10^2$ and relative displacement errors of $3 \times 10^2$ to $4 \times 10^2$---four orders of magnitude larger than the $3 \times 10^{-2}$ to $9 \times 10^{-2}$ error achieved by a non-negative matrix factorization (NMF) dual basis \cite{balajewicz2016projection} of identical rank. Crucially, the reduced solver alone provides no warning that this breakdown has occurred, motivating the computable lifted admissibility estimator derived in Section 3.

---

## v1 (original)

## 2. Problem Formulation and Silent Reduced-KKT Failure

### 2.1 Full-Order Discrete Unilateral Contact Problem

Consider a parametrized linear-elastic or quasi-static discrete mechanical system with $N$ displacement degrees of freedom and $n_c$ unilateral contact constraints ($n_c \le N$). For a parameter vector $\mu \in \mathcal{P} \subset \mathbb{R}^{n_p}$, the full-order equilibrium and non-penetration conditions are governed by the constrained quadratic minimization problem

$$
\min_{u \in \mathcal{K}(\mu)} \frac{1}{2} u^T K(\mu) u - f(\mu)^T u, \qquad \mathcal{K}(\mu) := \left\{ v \in \mathbb{R}^N : B(\mu) v - c(\mu) \le 0 \right\},
$$

where $K(\mu) \in \mathbb{R}^{N \times N}$ is the symmetric positive definite (SPD) stiffness matrix, $f(\mu) \in \mathbb{R}^N$ is the external load vector, $B(\mu) \in \mathbb{R}^{n_c \times N}$ is a rectangular constraint Jacobian mapping nodal displacements to normal contact gaps, and $c(\mu) \in \mathbb{R}^{n_c}$ is the initial gap vector. Throughout this work, inequalities on vectors are understood componentwise.

Introducing the vector of contact Lagrange multipliers $\lambda \in \mathbb{R}^{n_c}$ representing normal contact tractions, the first-order optimality conditions take the mixed Karush--Kuhn--Tucker (KKT) saddle-point form: find $(u(\mu), \lambda(\mu)) \in \mathbb{R}^N \times \mathbb{R}^{n_c}$ such that

$$
\begin{aligned}
K(\mu) u(\mu) + B(\mu)^T \lambda(\mu) &= f(\mu), \\
\lambda(\mu) &\ge 0, \\
B(\mu) u(\mu) - c(\mu) &\le 0, \\
\lambda(\mu)^T \bigl(B(\mu) u(\mu) - c(\mu)\bigr) &= 0.
\end{aligned}
$$

Equivalently, eliminating $u(\mu) = K(\mu)^{-1}(f(\mu) - B(\mu)^T \lambda(\mu))$ yields the dense full-order linear complementarity problem (LCP) for the multipliers:

$$
w(\mu) := B(\mu) K(\mu)^{-1} B(\mu)^T \lambda(\mu) + \bigl(c(\mu) - B(\mu) K(\mu)^{-1} f(\mu)\bigr) \ge 0, \qquad \lambda(\mu) \ge 0, \qquad \lambda(\mu)^T w(\mu) = 0.
$$

Because $K(\mu)$ is SPD, a unique primal solution $u(\mu)$ exists for every $\mu \in \mathcal{P}$. Uniqueness of $\lambda(\mu)$ holds whenever the active rows of $B(\mu)$ have full row rank, which is satisfied on standard non-redundant contact interfaces.

### 2.2 Galerkin--Petrov Mixed Reduced-Order Model

To construct a projection-based reduced-order model, we collect $n_s$ full-order snapshot pairs $\{(u(\mu_i), \lambda(\mu_i))\}_{i=1}^{n_s}$ over a training set $\mathcal{P}_{\text{train}} \subset \mathcal{P}$. From these snapshots, we build a primal displacement basis $V \in \mathbb{R}^{N \times r_u}$ with orthonormal columns ($V^T V = I_{r_u}$, $r_u \ll N$) and a dual multiplier basis $U \in \mathbb{R}^{n_c \times r_\lambda}$ ($r_\lambda \le n_c$).

Approximating the primal and dual variables in the reduced subspaces via

$$
u(\mu) \approx \tilde{u}(\mu) := V u_r(\mu), \qquad \lambda(\mu) \approx \tilde{\lambda}(\mu) := U \lambda_r(\mu),
$$

and testing the equilibrium and gap equations with $V^T$ and $U^T$, respectively, yields the reduced mixed KKT system: find $(u_r(\mu), \lambda_r(\mu)) \in \mathbb{R}^{r_u} \times \mathbb{R}^{r_\lambda}$ such that

$$
\begin{aligned}
K_r(\mu) u_r(\mu) + B_r(\mu)^T \lambda_r(\mu) &= f_r(\mu), \\
\lambda_r(\mu) &\ge 0, \\
B_r(\mu) u_r(\mu) - c_r(\mu) &\le 0, \\
\lambda_r(\mu)^T \bigl(B_r(\mu) u_r(\mu) - c_r(\mu)\bigr) &= 0,
\end{aligned}
$$

where the reduced operators are defined by

$$
K_r(\mu) := V^T K(\mu) V \in \mathbb{R}^{r_u \times r_u}, \quad B_r(\mu) := U^T B(\mu) V \in \mathbb{R}^{r_\lambda \times r_u}, \quad f_r(\mu) := V^T f(\mu) \in \mathbb{R}^{r_u}, \quad c_r(\mu) := U^T c(\mu) \in \mathbb{R}^{r_\lambda}.
$$

Because $V$ has full column rank and $K(\mu)$ is SPD, $K_r(\mu)$ is SPD, and the reduced system admits the Schur-complement LCP formulation

$$
w_r(\mu) := B_r(\mu) K_r(\mu)^{-1} B_r(\mu)^T \lambda_r(\mu) + \bigl(c_r(\mu) - B_r(\mu) K_r(\mu)^{-1} f_r(\mu)\bigr) \ge 0, \qquad \lambda_r(\mu) \ge 0, \qquad \lambda_r(\mu)^T w_r(\mu) = 0.
$$

### 2.3 Mechanisms of Silent Reduced-KKT Failure

A fundamental pitfall of the reduced formulation above is that enforcing the reduced KKT conditions to machine precision does **not** imply physical admissibility of the lifted pair $(\tilde{u}(\mu), \tilde{\lambda}(\mu))$. Specifically, define the reduced KKT residual vector

$$
R_{\text{red}}(u_r, \lambda_r) := \begin{bmatrix} \|K_r u_r + B_r^T \lambda_r - f_r\|_2 \\ \|\min(\lambda_r, 0)\|_2 \\ \|\max(B_r u_r - c_r, 0)\|_2 \\ |\lambda_r^T (B_r u_r - c_r)| \end{bmatrix}.
$$

Even when $\|R_{\text{red}}(u_r, \lambda_r)\|_\infty \le 10^{-16}$, the lifted solution can fail catastrophically through three coupled algebraic mechanisms:

1. **Loss of Dual Cone Preservation (Signed Dual Bases):** When $U$ is obtained via standard truncated SVD or POD of multiplier snapshots, its columns contain mixed positive and negative entries. Consequently, $\lambda_r \ge 0$ does not guarantee $\tilde{\lambda} = U \lambda_r \ge 0$. Negative components in $\tilde{\lambda}$ act as artificial adhesive tractions that pull the elastic body across the contact boundary to satisfy reduced equilibrium $V^T(K \tilde{u} + B^T \tilde{\lambda} - f) = 0$.
2. **Projected Penetration Masking:** The reduced non-penetration condition $U^T (B \tilde{u} - c) \le 0$ constrains only $r_\lambda$ weighted averages of the $n_c$ physical gaps. When $r_\lambda < n_c$, positive penetrations $(B \tilde{u} - c)_i > 0$ at active contact nodes can cancel against large negative gaps $(B \tilde{u} - c)_j < 0$ at inactive nodes in the inner product $U^T(B\tilde{u} - c)$.
3. **Reduced Inf-Sup Rank Deficiency:** Even when $B(\mu)$ has full row rank on the active contact set, the projected matrix $B_r(\mu) = U^T B(\mu) V$ can become rank-deficient or severely ill-conditioned if $V$ lacks sufficient normal-displacement resolution on the contact interface. In that case, $B_r K_r^{-1} B_r^T$ is singular or nearly singular, allowing large spurious vectors in $\ker(B_r^T)$ to enter $\lambda_r$ with negligible change to the reduced residual.

#### Illustrative Benchmark Evidence

To make this failure concrete, consider the parametrized 2D obstacle membrane ($N = 324$) and 1D discrete spring-chain contact ($N = 60$) benchmarks investigated in our preliminary measurements. When reduced with identical primal and dual dimensions $(r_u, r_\lambda)$ using a standard signed SVD dual basis, active-set and Lemke LCP solvers routinely converge to reduced KKT residuals below $10^{-16}$. Upon lifting, however, the signed SVD model exhibits lifted multipliers reaching $\min_i \tilde{\lambda}_i \approx -1.4 \times 10^2$ and relative displacement errors of $3 \times 10^2$ to $4 \times 10^2$—four orders of magnitude larger than the $3 \times 10^{-2}$ to $9 \times 10^{-2}$ error achieved by a non-negative matrix factorization (NMF) dual basis of identical rank. Crucially, the reduced solver alone provides no warning that this breakdown has occurred, motivating the computable lifted admissibility estimator derived in Section 3.
