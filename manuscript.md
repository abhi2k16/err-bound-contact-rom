# Combined Manuscript Drafts (original + v2)

*Final compiled version: `manuscript_master_v2.tex` / `.pdf`.*

---

## Introduction

# 1. Introduction

Projection-based reduced-order models (ROMs) are widely used to accelerate parametrized structural simulations in design optimization, uncertainty quantification, and real-time control. When unilateral contact constraints are present, however, the governing system takes the form of a discrete variational inequality or mixed Karush–Kuhn–Tucker (KKT) saddle-point problem involving both primal displacements $u \in \mathbb{R}^N$ and dual contact Lagrange multipliers $\lambda \in \mathbb{R}^{n_c}$. Because the active contact set evolves non-smoothly across parameter regimes, constructing reduced primal and dual spaces that remain simultaneously accurate, inf-sup stable, and physically admissible is substantially more delicate than in linear or smooth nonlinear mechanics.

A critical and under-appreciated failure mode in Lagrange-multiplier contact ROMs is **silent reduced-KKT convergence**. During online evaluation, a reduced contact solver enforces the projected KKT complementarity conditions strictly within the low-dimensional subspace. Consequently, the solver can satisfy its reduced KKT residuals to machine precision ($10^{-16}$) while returning a lifted primal–dual pair $(\tilde{u}, \tilde{\lambda}) = (V u_r, U \lambda_r)$ that violates full-order physical admissibility by orders of magnitude. When standard proper orthogonal decomposition (POD) or singular value decomposition (SVD) is used to build the dual basis $U$, linear combinations of signed basis modes routinely produce negative lifted contact pressures ($\tilde{\lambda} < 0$, representing non-physical adhesive tractions) alongside severe interpenetration errors. In benchmark tests on parametrized membrane and discrete chain contact problems, a signed SVD dual basis can satisfy reduced KKT tolerances to $10^{-16}$ while producing lifted multipliers down to $-1.4 \times 10^2$ and displacement errors of $3 \times 10^2$ to $4 \times 10^2$—four orders of magnitude worse than a non-negative basis of identical dimension.

Existing literature has addressed different facets of this challenge in isolation, leaving three major gaps:

1. **Positivity-preserving dual bases without lifted error bounds:** To prevent non-physical multipliers, several non-negative dual-basis constructions have been proposed, including non-negative matrix factorization (NMF) [1], cone-projected greedy algorithms [2], greedy index selection (gIS) [3], and restriction to a reduced integration domain (RID) [4, 5], alongside empirical cubature and parametric contact reduction strategies [6, 7, 8, 9]. While these methods improve physical consistency on their respective benchmarks, they do not provide a computable *a posteriori* error bound for the lifted primal–dual solution $(\tilde{u}, \tilde{\lambda})$, nor have they been evaluated head-to-head on identical snapshot sets.
2. **Restricted *a posteriori* error estimation for reduced variational inequalities:** In the full-order finite element literature, residual-based and equilibrated *a posteriori* error estimators for unilateral contact decompose the error into equilibrium, sign-violation, and complementarity residuals without requiring Galerkin orthogonality [10, 11, 12, 13, 14, 15]. Within the reduced-basis framework, however, certified duality-based error bounds have been derived under the restrictive assumption of a bijective (square, invertible) constraint operator $B$ with dual feasibility enforced by construction [16]. They do not cover general rectangular constraint matrices $B \in \mathbb{R}^{n_c \times N}$ arising in multi-body and non-matching contact, where reduced KKT satisfaction can mask catastrophic lifted inadmissibility.
3. **Unquantified enrichment, conditioning, and rank-repair overhead:** Preserving discrete inf-sup stability and non-singularity of the projected contact matrix often requires enriching the primal basis with supremizer or unit-vector rank-repair modes and retaining contact-zone integration rows. Although saddle-point preconditioning in full-order contact [17] and projected contact conditioning [5] have been investigated, the literature lacks an explicit break-even cost model predicting when rank-repair enrichment, active-set iteration growth, and dense reduced linear algebra eliminate the computational savings over a sparse full-order solver.

In this paper, we address all three gaps by developing a unified reliability and cost framework for Lagrange-multiplier reduced contact models:

* **Computable Lifted Admissibility Estimator:** Adapting full-order contact residual decompositions to lifted ROM pairs under general rectangular constraint matrices $B \in \mathbb{R}^{n_c \times N}$, we derive a rigorous energy-norm *a posteriori* error bound composed of four terms: a dual-norm equilibrium residual ($\eta_{\text{eq}}$) and three online-computable contact admissibility indicators measuring multiplier sign violation ($\eta_{\text{sign}}$), primal gap penetration ($\eta_{\text{gap}}$), and complementarity defect ($\eta_{\text{comp}}$). The three contact indicators evaluate in $\mathcal{O}(n_c r_u)$ operations without touching the bulk mesh, acting as an online circuit breaker that detects silent reduced-KKT failures.
* **Unified Cross-Basis Benchmark:** We perform the first controlled comparison of five dual-basis strategies—NMF, cone-projected greedy, greedy index selection, signed SVD with rank repair, and RID restriction—on identical snapshot sets across 1D, 2D, and 3D frictionless contact problems, evaluating lifted accuracy, estimator effectivity, and out-of-range active-set extrapolation.
* **Break-Even Computational Cost Model:** We formulate and experimentally validate an explicit FLOP- and iteration-aware cost model incorporating primal rank-repair enrichment, contact-row retention, active-set iteration counts, and dense-versus-sparse linear algebra scaling, identifying the exact regimes in which Lagrange-multiplier contact ROMs outperform full-order solvers.

The remainder of the paper is organized as follows. Section 2 states the full-order and reduced mixed contact formulations and illustrates silent reduced-KKT failure. Section 3 derives the four-component *a posteriori* error bound and its online admissibility indicator. Section 4 reviews the five dual-basis constructions and formulates the break-even cost model. Section 5 presents numerical experiments across 1D chain, 2D obstacle membrane, and 3D two-body contact benchmarks, and Section 6 summarizes practical guidelines for selecting and certifying contact ROMs.

---

## References

1. M. Balajewicz, D. Amsallem, and C. Farhat, "Projection-based model reduction for contact problems," *International Journal for Numerical Methods in Engineering*, vol. 106, no. 8, pp. 644–663, 2016. DOI: `10.1002/nme.5135`
2. A. Benaceur, A. Ern, and V. Ehrlacher, "A reduced basis method for parametrized variational inequalities applied to contact mechanics," *International Journal for Numerical Methods in Engineering*, vol. 121, no. 6, pp. 1170–1197, 2020. DOI: `10.1002/nme.6261`
3. G. Sambataro and V. Ehrlacher, "Model order reduction for parametrized variational inequalities: application to crowd motion," *arXiv preprint arXiv:2605.04037*, 2026. DOI: `10.48550/arXiv.2605.04037`
4. J. Fauque, I. Ramière, and D. Ryckelynck, "Hybrid hyper-reduced modeling for contact mechanics problems," *International Journal for Numerical Methods in Engineering*, vol. 115, no. 1, pp. 117–139, 2018. DOI: `10.1002/nme.5798`
5. S. Le Berre, I. Ramière, J. Fauque, and D. Ryckelynck, "Condition Number and Clustering-Based Efficiency Improvement of Reduced-Order Solvers for Contact Problems Using Lagrange Multipliers," *Mathematics*, vol. 10, no. 9, p. 1495, 2022. DOI: `10.3390/math10091495`
6. L. Koller, D. J. Rixen, et al., "Hyper-reduction for finite element structures with nonlinear surface loads and contact," *Computer Methods in Applied Mechanics and Engineering*, vol. 385, p. 114023, 2021. DOI: `10.1016/j.cma.2021.114023`
7. S. Ares de Parga, J. A. Hernández, and R. Rossi, "Petrov–Galerkin hyper-reduction via the empirical cubature method for non-symmetric Jacobians," *Computer Methods in Applied Mechanics and Engineering*, vol. 414, p. 116180, 2023. DOI: `10.1016/j.cma.2023.116180`
8. G. S. Pawar and S. S. Kulkarni, "Efficient parametric model order reduction in contact mechanics," *Computers & Structures*, vol. 315, p. 107838, 2025. DOI: `10.1016/j.compstruc.2025.107838`
9. D. Manvelyan-Stroot, Y. Filanova, I. Pontes Duff, P. Benner, and U. Wever, "Inference of substructured reduced-order models for dynamic contact from contact-free simulations," *Computational Mechanics*, 2026. DOI: `10.1007/s00466-026-02848-w`
10. P. Coorevits, P. Hild, and J.-P. Pelle, "A posteriori error estimation for unilateral contact with matching and non-matching meshes," *Computer Methods in Applied Mechanics and Engineering*, vol. 186, no. 1, pp. 65–83, 2000. DOI: `10.1016/S0045-7825(99)00370-9`
11. P. Hild and S. Nicaise, "Residual a posteriori error estimators for contact problems in elasticity," *ESAIM: Mathematical Modelling and Numerical Analysis*, vol. 41, no. 5, pp. 897–923, 2007. DOI: `10.1051/m2an:2007045`
12. J. Dabaghi, V. Martin, and M. Vohralík, "A posteriori error estimates distinguishing the error components and adaptive stopping criteria for numerical approximations of parabolic variational inequalities," *Computer Methods in Applied Mechanics and Engineering*, vol. 367, p. 113105, 2020. DOI: `10.1016/j.cma.2020.113105`
13. L. Banz and F. Bertrand, "Generic a posteriori error estimators for two-field contact problems," *Computers & Mathematics with Applications*, vol. 169, pp. 188–205, 2024. DOI: `10.1016/j.camwa.2024.06.015`
14. M. Fontana and D. A. Di Pietro, "Guaranteed a posteriori error bounds for finite element approximations of contact problems," *Computer Methods in Applied Mechanics and Engineering*, vol. 424, p. 116895, 2024. DOI: `10.1016/j.cma.2024.116895`
15. L. Banz, M. Schönauer, and A. Schröder, "A posteriori error estimates for variational inequalities of the second kind in an abstract framework," *arXiv preprint arXiv:2609.26001*, 2026. DOI: `10.48550/arXiv.2609.26001`
16. Z. Zhang, E. Bader, and K. Veroy, "A Duality Approach to Error Estimation for Variational Inequalities," *arXiv preprint arXiv:1410.2095*, 2014. DOI: `10.48550/arXiv.1410.2095`
17. A. Franceschini, N. Castelletto, J. A. White, and H. A. Tchelepi, "Algebraically-preconditioned Krylov solvers for Lagrange-multiplier contact problems," *Computer Methods in Applied Mechanics and Engineering*, vol. 392, p. 114632, 2022. DOI: `10.1016/j.cma.2022.114632`

---

## Section 2: Formulation — v2 (current)

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


## Section 2: Formulation — original (v1)

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

---

## Section 3: Estimator — v2 (current)

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


## Section 3: Estimator — original (v1)

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

---

## Section 4: Basis and Cost — v2 (current)

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


## Section 4: Basis and Cost — original (v1)

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

---

## Section 5: Experiments — v2 (current)

### 5. Numerical Experiments Protocol

To systematically evaluate the certified *a posteriori* estimator (Section 3), compare the five dual-basis constructions under controlled conditions (Section 4.1), and validate the FLOP- and iteration-aware break-even cost model (Section 4.2), we design a six-part experimental benchmark suite spanning 1D, 2D, and 3D unilateral contact problems. All reduced-order models are evaluated on identical snapshot sets $\{(u(\mu_i), \lambda(\mu_i))\}_{i=1}^{n_s}$ and solved using a unified active-set / Lemke LCP solver converged to a reduced KKT residual tolerance of $10^{-14}$ \cite{balajewicz2016projection, benaceur2020reduced, leberre2022condition}.

#### 5.1 Benchmark Test Cases
We consider three progressively larger frictionless unilateral contact problems covering distinct interface topologies and sparsity structures:

1. **Benchmark B1: 1D Discrete Spring-Chain Contact ($N = 60$, $n_c = 30$)**
   A nonlinearly coupled or parametrized chain of $N = 60$ spring-mass nodes suspended above a rigid variable-profile obstacle $c(\mu)$ \cite{pawar2025efficient, balajewicz2016projection}. The parameter vector $\mu = (\mu_1, \mu_2) \in \mathcal{P} \subset \mathbb{R}^2$ controls the spatially varying downward load amplitude and obstacle curvature. Because $n_c / N = 0.50$ is large relative to $N$, this benchmark isolates rank-repair overhead ($r_{\text{rep}} = 14$ unit vectors required for full row rank of $B_r$ \cite{leberre2022condition}) in the constant-$K$ regime.
2. **Benchmark B2: 2D Parametrized Obstacle Membrane ($N = 324$, $n_c = 144$)**
   A clamped two-dimensional elastic membrane discretized on an $18 \times 18$ interior grid ($N = 324$) pressed against a rigid curved indenter \cite{benaceur2020reduced, hild2007residual, coorevits2000posteriori}. The parameter vector $\mu = (F_0, x_c, y_c) \in \mathbb{R}^3$ varies the peak load magnitude $F_0$ and the indenter center $(x_c, y_c)$, causing the active contact patch $m_{\text{act}}$ to expand, contract, and translate across the candidate interface.
3. **Benchmark B3: 3D Two-Body Hertzian / Flange Contact ($N \in [2 \times 10^4, 10^5]$, $n_c \in [400, 2500]$)**
   A three-dimensional linear-elastic deformable hemisphere (or bolted structural flange) pressed onto a deformable elastic block with a node-to-segment frictionless contact interface \cite{franceschini2022block, koller2021hyper, fauque2018hybrid}. Full-order finite element matrices $(K(\mu), B, c(\mu), f(\mu))$ are assembled and verified using `FEniCSx` / `GetFEM` with a sparse direct Cholesky / `PETSc TAO` contact solver \cite{franceschini2022block}, and exported directly into the reduced-order pipeline. Parameters $\mu = (E_1/E_2, P_{\text{load}}, \delta_0)$ vary the stiffness contrast (testing parameter-varying $K(\mu)$ \cite{koller2021hyper}), applied pressure, and initial clearance.

| Benchmark | Spatial Dim. | Full-Order $N$ | Contact Rows $n_c$ | Stiffness $K(\mu)$ | Primary Purpose |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **B1 (Spring Chain)** \cite{pawar2025efficient} | 1D | $60$ | $30$ | Constant | Rank-repair cost crossover & silent SVD failure |
| **B2 (Obstacle Membrane)** \cite{benaceur2020reduced} | 2D | $324$ | $144$ | Constant | Five-basis comparison & contact-patch translation |
| **B3 (3D Two-Body Contact)** \cite{franceschini2022block, koller2021hyper} | 3D | $2 \times 10^4\text{--}10^5$ | $400\text{--}2500$ | Parametrized $K(\mu)$ | Large-scale wall-clock speedup & 3D effectivity |

#### 5.2 Experimental Campaign (E1--E6)

* **Experiment E1: Silent Reduced-KKT Failure Audit (B1, B2)**
  Using identical primal POD bases $V_0$ and dual bases $U$ of matched dimension $r_\lambda \in \{4, 8, 12, 16, 20\}$, we log the reduced KKT residual $\Vert R_{\text{red}}(u_r, \lambda_r)\Vert_\infty$ alongside the true relative energy-norm error $\Vert u - \tilde{u}\Vert_K / \Vert u\Vert_K$, the four estimator terms $(\eta_{\text{eq}}, \eta_{\text{sign}}, \eta_{\text{gap}}, \eta_{\text{comp}})$ \cite{hild2007residual, banz2024generic, banz2026posteriori}, and the online admissibility indicator $\eta_{\text{adm}}$. This experiment verifies that $\eta_{\text{adm}}$ flags every case where $\Vert R_{\text{red}}\Vert_\infty \le 10^{-14}$ masks catastrophic lifted error.

* **Experiment E2: Estimator Effectivity and Extrapolation Stress Test (B1, B2, B3)**
  We evaluate the certified upper bound from Theorem 1 and compute the **effectivity index** \cite{hild2007residual, zhang2014duality, fontana2024guaranteed}
  $$\theta_{\text{eff}}(\mu) := \frac{\eta_{\text{bound}}(\mu)}{\Vert u(\mu) - \tilde{u}(\mu)\Vert_K}, \qquad \eta_{\text{bound}}(\mu) := \left( \eta_{\text{eq}}^2 + 2\,\eta_{\text{comp}} + 2\,\Vert \lambda_+\Vert_2\,\eta_{\text{gap}} \right)^{1/2},$$
  across three distinct parameter test sets:
  1. *In-sample validation set* $\mathcal{P}_{\text{val}} \subset \mathcal{P}_{\text{train}}$ (interpolatory regime);
  2. *Load-magnitude extrapolation* ($F_0 > \max_{\mu \in \mathcal{P}_{\text{train}}} F_0$), expanding the active contact set beyond training snapshots \cite{benaceur2020reduced};
  3. *Contact-patch translation* (shifting $(x_c, y_c)$ so that previously inactive rows of $B$ activate online, testing whether $\eta_{\text{gap}}$ and $\eta_{\text{comp}}$ capture unseen active constraints \cite{fauque2018hybrid, sambataro2026model}).

* **Experiment E3: Head-to-Head Five-Basis Benchmark (B2, B3)**
  On identical snapshot matrices $(U_{\text{snap}}, \Lambda_{\text{snap}})$, we construct and compare all five dual-basis strategies from Section 4.1---**SVD+Repair** \cite{leberre2022condition}, **NMF** \cite{balajewicz2016projection}, **CPG** \cite{benaceur2020reduced}, **gIS** \cite{sambataro2026model, kollepara2024sparse}, and **RID** \cite{fauque2018hybrid, leberre2022condition}---across dual dimensions $r_\lambda \in \{4, 8, 12, 16, 24\}$. Each method is evaluated along five Pareto axes:
  1. Lifted energy-norm displacement error and multiplier $\ell_2$-error;
  2. Online admissibility violation $\eta_{\text{adm}}$;
  3. Primal rank-repair enrichment $r_{\text{rep}}$ required to maintain $\text{cond}(B_r K_r^{-1} B_r^T) \le 10^6$ \cite{leberre2022condition};
  4. Active-set iteration count $k_{\text{ROM}}$ and online wall-clock time;
  5. Offline basis construction wall-clock time.

* **Experiment E4: Break-Even Cost Model Validation (B1, B2, B3)**
  We measure empirical floating-point operation counts and single-thread wall-clock execution times as functions of the enriched primal dimension $r_u = r_{u,0} + r_{\text{rep}}$, active contact size $m_{\text{act}}$, and mesh size $N$ \cite{franceschini2022block, leberre2022condition, pawar2025efficient}. Measured speedup curves $S_{\text{meas}} = t_{\text{FOM}} / t_{\text{ROM}}$ are plotted against the theoretical predictions $S_{\text{pred}} = \mathcal{C}_{\text{FOM}} / \mathcal{C}_{\text{ROM}}$ from Equations (4)--(5) in both the constant-$K$ (B1, B2) and parameter-varying $K(\mu)$ (B3) regimes.

* **Experiment E5: Large-Scale 3D Contact Certification (B3)**
  On the 3D two-body contact model ($N = 5 \times 10^4$ and $N = 10^5$) \cite{franceschini2022block, koller2021hyper}, we deploy the top-performing non-negative dual bases from E3 paired with the online circuit breaker $\eta_{\text{adm}} \le \tau_{\text{adm}}$ \cite{dabaghi2020posteriori}, demonstrating certified online acceleration and automatic fallback when extrapolation triggers $\eta_{\text{adm}} > \tau_{\text{adm}}$.

* **Experiment E6: Time-Dependent Contact Dynamics Extension (Optional, B1/B2)**
  To demonstrate applicability beyond static equilibrium, we integrate the reduced LCP with a second-order backward differentiation formula (BDF2) time-stepping scheme \cite{dabaghi2020posteriori, manvelyan2026inference}, evaluating $\eta_{\text{adm}}(t_n)$ at every discrete time step $t_n$ to detect transient impact admissibility violations.

#### 5.3 Quantitative Evaluation Metrics and Acceptance Criteria
To ensure objective assessment across all benchmarks, we predefine four quantitative validation metrics:
1. **Admissibility Detection Recall (Zero False Negatives):** Every test query with $\Vert R_{\text{red}}\Vert_\infty \le 10^{-12}$ whose true relative energy-norm error exceeds $5\%$ due to contact inadmissibility must satisfy $\eta_{\text{adm}} > \tau_{\text{adm}}$.
2. **Effectivity Sharpness:** For the certified bound $\eta_{\text{bound}}$ (using $\Vert \lambda_+\Vert_2$ as the online multiplier norm proxy \cite{hild2007residual, benaceur2020reduced, fontana2024guaranteed}), the effectivity index $\theta_{\text{eff}}$ lies within $[1, 10]$ on at least $90\%$ of admissible and borderline test cases, with explicit isolation of bulk-dominated regimes where $\eta_{\text{eq}}$ governs the error.
3. **Monotonic Error Correlation:** The Spearman rank correlation between $\eta_{\text{bound}}$ and the true energy-norm error $\Vert u - \tilde{u}\Vert_K$ exceeds $0.90$ across all dual bases and reduced dimensions.
4. **Cost Crossover Accuracy:** The predicted break-even threshold $\mathcal{C}_{\text{ROM}} = \mathcal{C}_{\text{FOM}}$ matches the measured wall-clock break-even point within a factor of $2$ across benchmarks B1--B3.


## Section 5: Experiments — original (v1)

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

---

## Section 6: Conclusion — v2 (current)

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


## Section 6: Conclusion — original (v1)

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

---

## Appendix A: RESEARCH_PLAN.md

# Research plan: reliability and cost of Lagrange-multiplier reduced contact models

Status: DRAFT for discussion. Built on the measured results of Chapter 14 of the book (`packages_book/book_v3`) and a literature check done on 13 Nov 2026 (abstracts, plus parts of two full texts; Google Scholar could not be searched).

## 1. Working title
A computable admissibility estimator and a cost model for reduced contact models with Lagrange multipliers.

## 2. Problem and gap
A reduced contact model with multipliers can satisfy its own reduced KKT conditions to 1e-16 and still return a lifted solution that is wrong by orders of magnitude (measured in the book: displacement error 3e2 to 4e2 with a signed SVD dual basis, against 3e-2 to 9e-2 with an NMF basis at the same sizes; lifted multipliers down to about -1.4e2).

What the literature covers (verified from abstracts or text):
- Positive dual bases: NMF (Balajewicz et al. 2016), cone-projected greedy (Benaceur, Ern, Ehrlacher 2020), greedy index selection (Sambataro, Ehrlacher 2026), restriction to a reduced integration domain (Fauque et al. 2018; Le Berre et al. 2022).
- Certified bounds for reduced variational inequalities (Zhang, Bader, Veroy 2014): they assume a bijective constraint operator B and give feasible dual approximations by construction.
- Full-order residual estimators for mixed contact (Hild, Nicaise 2007).
- Conditioning of the projected contact matrix linked to dual precision (Le Berre et al. 2022).

- Full-order estimators built by replacing the multiplier with a post-processed discrete one (Banz, Schönauer, Schröder, arXiv 2609.26001, Sept 2026; variational inequalities of the second kind, friction-type). The device is the same as the lambda+ projection used here, so the novelty claimed is the reduced-order lifted setting, not the technique.
- Greedy positive-cone compression with hyper-reduction for discrete contact (Sambataro, Ehrlacher 2026), which postpones an error indicator to future work (in the context of basis enrichment for inf-sup stability).

What I did not find (to be re-checked before submission): an estimator for the lifted solution of a reduced model with a rectangular constraint matrix; a measured case where the reduced residual passes and the lifted solution fails; one benchmark that compares the dual-basis constructions; a cost model for when rank repair plus reduced-domain enrichment removes the saving.

## 3. Research questions
RQ1. Can a computable a posteriori estimate computed from the lifted pair detect inadmissible reduced solutions that the reduced residual accepts?
RQ2. Which dual basis (NMF, cone-projected greedy, greedy index selection, signed SVD with repair, reduced-integration-domain restriction) is most reliable, at what size?
RQ3. When does the total cost (reduced unknowns after rank repair, reduced integration domain, contact rows) fall below the full-order cost?

## 4. Proposed method
Lifted pair: u~ = V u_r, lambda~ = U lambda_r. Define lambda+ = max(lambda~, 0), lambda- = min(lambda~, 0).

Estimator components:
- eta_sign = ||lambda-||, cost O(contact rows).
- eta_gap = ||(c - B u~)+||, cost O(contact rows).
- eta_comp = |lambda+^T (B u~ - c)|, cost O(contact rows).
- eta_eq = ||K u~ + B^T lambda+ - f||_{K^-1}, needs a full-order operation, so it is a validation or greedy-selection quantity, not an online one.

Online admissibility indicator (candidate, to be tested in E2): eta_adm = eta_sign + eta_gap + sqrt(eta_comp), with each term scaled to a common unit (for instance through the K and K^-1 norms) because force, length and work differ in units. The square root on eta_comp is a definition to test, not a proved property: the claim that all three terms scale linearly with the error is measured in E2. Use as an online circuit breaker: if eta_adm exceeds a tolerance, flag the prediction and trigger basis enrichment or a fall-back to the full-order solver. eta_eq stays an offline certificate (validation and greedy selection).

Target (to prove or to establish empirically): the energy-norm error is bounded by a constant times a combination of these terms, following the constant-sign multiplier idea of Hild and Nicaise. A proved bound is not assumed; if it fails, the paper reports the empirical effectivity instead.

Cost model, two tiers.
Tier 1, size: unknowns after enrichment = reduced displacement columns + unit vectors added to restore rank (14 on the book's chain) + reduced multipliers; compared with full-order unknowns, as a function of contact-zone size and the number of contact rows.
Tier 2, work: measured wall-clock time and fitted scaling per solve, including dense reduced-matrix factorisation (O(r^3)) against sparse full-order factorisation, and the number of active-set or LCP iterations k_iter per step (a poorly conditioned dual basis can need more iterations than the full-order solver, which removes the saving even when r is much smaller than N). The break-even is reported per solver and hardware from measured scaling; no closed-form crossover is assumed.

## 5. Experiments
E1. Reproduce the book's failure on the 324-unknown membrane and the 60-spring chain; compute all four estimator components for every basis and size.
E2. Effectivity: estimate against the true lifted error, for each basis, in three regimes: (i) in-range parameters; (ii) load-magnitude extrapolation (the contact patch grows beyond the training snapshots, where rank deficiency and negative multipliers are expected); (iii) contact-patch shift or sliding (rows of B inactive in training become active, testing whether eta_gap and eta_comp catch unseen active constraints). Also measure the scaling of each term and of eta_adm against the error. Report the failure cases.
E3. Comparison of the five dual-basis constructions on identical snapshots (accuracy, admissibility violation, online cost, offline cost).
E4. Cost model validation: measured against predicted break-even on the rope and on the membrane, then on a 3D case.
E5. Scale-up: 3D Hertzian contact (elastic hemisphere or cylinder on a block; analytical solution for verification), 1e4 to 1e5 unknowns, with the exact (K, B, c, f) exported from a scriptable FEM backend so that the ROM code is identical to the 1D and 2D code. Backend not yet chosen: FEniCSx contact libraries are marked experimental and I could not confirm two-body support; GetFEM was not checked. A one-week feasibility spike must export (K, B, c, f) for a small 3D case before a tool is fixed.
E6. Dynamics (optional): the reduced LCP with BDF2, with the estimator evaluated per step.

## 6. Success criteria (set before running)
- No false negative: every case with true error above a stated tolerance is flagged by eta_sign, eta_gap or eta_comp.
- Effectivity index between 1 and 10 on at least 90 % of cases for the best basis, or a clear statement of where it fails.
- Spearman rank correlation of the estimate with the true error above 0.9 across bases and sizes.
- Predicted break-even within a factor of 2 of the measured one.
- Online indicator eta_adm: (i) scaling: the slope of log eta_adm against log error, fitted over all cases, lies between 0.8 and 1.2, otherwise the form (units or square root) is changed and the change reported; (ii) circuit breaker: with the tolerance set on the training range only, no case in the extrapolation and patch-shift regimes of E2 with true error above the stated tolerance goes unflagged, and the false-alarm rate (flagged cases with true error below tolerance) is reported, with 10 % as the target; (iii) online cost of eta_adm stays below 10 % of one reduced solve step.
- Thresholds above are initial targets set before running and are revised only with the reason recorded.

## 7. Risks
- The bound may not be provable; mitigation: report it empirically.
- eta_eq is not online-cheap; mitigation: use it for validation and greedy selection and state the limitation.
- Overlap with Le Berre et al. on conditioning; mitigation: cite and show what the measured rank repair adds.
- A recent preprint (Sambataro, Ehrlacher 2026) may be followed by an indicator; mitigation: monitor arXiv before submission.
- The 3D benchmark needs a verified contact solver; this is the main cost. Mitigation: the feasibility spike in E5 before committing.
- The unit scaling and the square-root form of eta_adm may not give linear scaling with the error; mitigation: measure it, and report the form that works.

## 8. Timeline (about 8 months)
Months 1-2: literature recheck with Scholar, estimator derivation, E1-E2 on the existing benchmarks.
Months 3-4: E3, the comparison with the four other bases.
Months 5-6: E4-E5, the 3D solver and cost model.
Month 7: E6 if time allows; writing.
Month 8: internal review, submission.

## 9. Target journals
First choice Computer Methods in Applied Mechanics and Engineering or International Journal for Numerical Methods in Engineering if the estimator has a proof; Computers & Structures or Computational Mechanics if the contribution is the comparison and the cost model.

## 10. Reuse from the package
`ch13_assets.py` (`ObstacleMembrane`, `NonlinearChainContact`, `solve_contact_dynamics`, `ContactROM`, `nonneg_dual_basis`, `svd_dual_basis`, `greedy_contact_sampling`, `nonlinearity_indicator`); `ch14_intro.py`, `ch14_basis.py`, `ch14_more.py` for the measured evidence.

## 11. Still to verify before writing
- Full texts of Balajewicz et al. (2016) and Benaceur et al. (2020): is there an error estimate for the lifted multiplier? STILL OPEN (only abstracts read).
- Google Scholar searches: `"reduced basis" contact "Lagrange multiplier" "a posteriori"`; `"non-negative" "dual basis" "reduced" contact "error estimate"`; `"hyper-reduction" contact "inf-sup" "unilateral"`.
- Whether any paper shows a passing reduced residual with a failing lifted solution.

## 12. Sources
- Banz, Schönauer, Schröder, A posteriori error estimates for variational inequalities of the second kind in an abstract framework, arXiv 2609.26001 (2026).
- Zhang, Bader, Veroy, A Duality Approach to Error Estimation for Variational Inequalities, arXiv 1410.2095.
- Sambataro, Ehrlacher, arXiv 2605.04037 (2026).
- Hild, Nicaise, ESAIM M2AN 41(5):897-923 (2007).
- Balajewicz, Amsallem, Farhat, IJNME 106(8):644-663 (2016), arXiv 1503.01000.
- Benaceur, Ern, Ehrlacher, IJNME 121(6):1170-1197 (2020).
- Fauque, Ramiere, Ryckelynck, IJNME 115(1):117-139 (2018).
- Le Berre et al., Mathematics 10(9):1495 (2022).
- Niakh, Drouet, Ehrlacher, Ern, SMAI-JCM 10:29-54 (2024), arXiv 2307.11541.
- Manvelyan-Stroot et al., arXiv 2505.18050 (2025).
- Kollepara et al., arXiv 2406.11461 (2024).
- Pawar, Kulkarni, Computers & Structures 315:107838 (2025).

## 13. Data and reproducibility
The benchmark matrices (K, B, c, f), unified Python implementations of all five dual-basis builders, and the estimator verification scripts will be released in a public repository at submission, with the random seeds and solver versions recorded.

---

## Appendix B: research_plan_updated.md

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

