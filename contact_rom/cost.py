"""FLOP- and iteration-aware break-even cost model (manuscript Section 4.2, Eqs. 16-17)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp


@dataclass(frozen=True)
class ProblemSizes:
    """Full-order quantities entering Eq. (16) and the projection term of Eq. (17)."""

    n_dof: int          # N
    n_contact: int      # n_c
    bandwidth: int      # b_w of K
    nnz_K: int
    nnz_B: int


@dataclass(frozen=True)
class ReducedSizes:
    """Reduced quantities entering Eq. (17)."""

    r_u: int            # enriched primal dimension r_{u,0} + r_rep
    r_lambda: int       # dual dimension actually used
    k_rom: float        # active-set iterations (reduced)
    m_r_act: float      # reduced active-set size


def matrix_bandwidth(K_mat: sp.spmatrix) -> int:
    """Half bandwidth ``max |i - j|`` over the nonzeros of ``K_mat``."""
    coo = sp.coo_matrix(K_mat)
    return int(np.max(np.abs(coo.row - coo.col))) if coo.nnz else 0


def problem_sizes(K_mat: sp.spmatrix, B_mat: sp.spmatrix) -> ProblemSizes:
    """Collect :class:`ProblemSizes` from the sparse operators."""
    return ProblemSizes(K_mat.shape[0], B_mat.shape[0], matrix_bandwidth(K_mat),
                        int(sp.csr_matrix(K_mat).nnz), int(sp.csr_matrix(B_mat).nnz))


def flops_fom(sizes: ProblemSizes, k_fom: float, m_act: float, c_fact: float, c_solve: float,
              k_varies: bool) -> float:
    """Eq. (16). With constant K the factorisation and compliance assembly are offline (affine form)."""
    lcp = k_fom * ((2.0 / 3.0) * m_act**3 + 2.0 * sizes.n_contact * m_act)
    if not k_varies:
        return float(lcp)
    factor = c_fact * sizes.n_dof * sizes.bandwidth**2
    compliance = 2.0 * c_solve * sizes.n_contact * sizes.n_dof * sizes.bandwidth
    return float(factor + compliance + lcp)


def flops_projection(sizes: ProblemSizes, red: ReducedSizes, c_proj: float, k_varies: bool) -> float:
    """``C_proj``: zero when reduced operators are precomputed (constant K); otherwise
    ``c_proj (nnz(K) r_u + nnz(B)(r_u + r_lambda))`` for online projection of full-order matrices."""
    if not k_varies:
        return 0.0
    return float(c_proj * (sizes.nnz_K * red.r_u + sizes.nnz_B * (red.r_u + red.r_lambda)))


def flops_rom(sizes: ProblemSizes, red: ReducedSizes, c_proj: float, k_varies: bool) -> float:
    """Eq. (17): projection + dense solves + reduced active-set + online eta_adm check.

    Constant K: the Cholesky factor of ``A_r``, ``A_r^-1 B_r^T`` and the reduced Schur complement are offline, so
    the online dense work is two triangular solves for ``A_r^-1 b`` and one product with ``B_r``. Varying K:
    factorisation, ``A_r^-1 B_r^T`` and the Schur complement are recomputed for every query.
    """
    r_u, r_l = red.r_u, red.r_lambda
    if k_varies:
        dense = r_u**3 / 3.0 + 2.0 * r_l * r_u**2 + 2.0 * r_l**2 * r_u
    else:
        dense = 2.0 * r_u**2 + 2.0 * r_l * r_u
    lcp = red.k_rom * ((2.0 / 3.0) * red.m_r_act**3 + 2.0 * r_l * red.m_r_act)
    check = 2.0 * sizes.n_contact * (r_u + r_l)
    return float(flops_projection(sizes, red, c_proj, k_varies) + dense + lcp + check)


def speedup_pred(sizes: ProblemSizes, red: ReducedSizes, k_fom: float, m_act: float, c_fact: float,
                 c_solve: float, c_proj: float, k_varies: bool) -> float:
    """Predicted speedup ``S = C_FOM / C_ROM``."""
    return flops_fom(sizes, k_fom, m_act, c_fact, c_solve, k_varies) / flops_rom(sizes, red, c_proj, k_varies)


def find_crossover(x_vals: np.ndarray, speedups: np.ndarray) -> float:
    """Smallest ``x`` where ``S`` crosses 1 (log-log interpolation); NaN if it never crosses in range."""
    log_s = np.log(speedups)
    for i in range(len(x_vals) - 1):
        if log_s[i] <= 0.0 < log_s[i + 1] or log_s[i] >= 0.0 > log_s[i + 1]:
            frac = -log_s[i] / (log_s[i + 1] - log_s[i])
            return float(np.exp(np.log(x_vals[i]) + frac * (np.log(x_vals[i + 1]) - np.log(x_vals[i]))))
    return float("nan")


def find_break_even(x_vals: np.ndarray, speedups: np.ndarray) -> float:
    """Break-even size: the size from which ``S`` stays above 1 (log-log interpolation after the last point with ``S <= 1``).

    NaN if ``S > 1`` over the whole range or ``S <= 1`` at the largest size. Unlike ``find_crossover`` this ignores
    early noisy sign changes: for monotone curves the two agree.
    """
    log_s = np.log(speedups)
    below = np.where(log_s <= 0.0)[0]
    if below.size == 0 or below[-1] == len(x_vals) - 1:
        return float("nan")
    i = int(below[-1])
    frac = -log_s[i] / (log_s[i + 1] - log_s[i]) if log_s[i + 1] != log_s[i] else 0.0
    return float(np.exp(np.log(x_vals[i]) + frac * (np.log(x_vals[i + 1]) - np.log(x_vals[i]))))
