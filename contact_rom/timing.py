"""Wall-clock measurement of full-order and reduced online solves (single thread, Python)."""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from rom_engine.complementarity import lcp_active_set

from .benchmarks import ContactProblem


@dataclass
class FomTimings:
    """Per-query full-order measurements for both stiffness regimes."""

    t_affine: np.ndarray    # K constant: K^-1 and B K^-1 B^T precomputed offline
    t_full: np.ndarray      # K parameter-varying: factorisation + compliance assembly online
    k_fom: np.ndarray       # active-set iterations
    m_act: np.ndarray       # number of active constraints


def best_of(func, repeats: int) -> float:
    """Minimum wall time of ``repeats`` calls to ``func`` (seconds)."""
    best = np.inf
    for _ in range(repeats):
        start = time.perf_counter()
        func()
        best = min(best, time.perf_counter() - start)
    return float(best)


def fom_solve_with_factor(lu: spla.SuperLU, Kinv_Bt: np.ndarray, A_lcp: np.ndarray, B_s: sp.spmatrix,
                          f_vec: np.ndarray, c_s: np.ndarray) -> tuple:
    """Online full-order solve given offline factors; returns ``(u, lam, iterations)``."""
    Kinv_f = lu.solve(f_vec)
    lam, _, iters = lcp_active_set(A_lcp, B_s @ Kinv_f - c_s)
    return Kinv_f + Kinv_Bt @ lam, lam, iters


def assemble_factors(K_mat: sp.spmatrix, B_s: sp.spmatrix) -> tuple:
    """``(lu, K^-1 B^T, A = B K^-1 B^T)`` for the suite sign convention."""
    lu = spla.splu(sp.csc_matrix(K_mat))
    Kinv_Bt = lu.solve(B_s.T.toarray())
    A_lcp = B_s @ Kinv_Bt
    return lu, Kinv_Bt, 0.5 * (A_lcp + A_lcp.T)


def time_fom_queries(problem: ContactProblem, mu_list: np.ndarray, repeats: int) -> FomTimings:
    """Time the FOM per query in the affine and the K-varying regime (best of ``repeats``)."""
    K_mat, B_s = problem.K_mat(mu_list[0]), -problem.B_mat(mu_list[0])
    offline = assemble_factors(K_mat, B_s)
    t_aff, t_full, k_fom, m_act = [], [], [], []
    for mu in mu_list:
        f_vec, c_s = problem.f_vec(mu), -problem.c_vec(mu)
        t_aff.append(best_of(lambda: fom_solve_with_factor(*offline, B_s, f_vec, c_s), repeats))
        t_full.append(best_of(lambda: fom_solve_with_factor(*assemble_factors(K_mat, B_s), B_s, f_vec, c_s), repeats))
        _, lam, iters = fom_solve_with_factor(*offline, B_s, f_vec, c_s)
        k_fom.append(iters)
        m_act.append(int(np.sum(lam > 1e-12)))
    return FomTimings(np.array(t_aff), np.array(t_full), np.array(k_fom, float), np.array(m_act, float))


def time_eta_check(B_mat: sp.spmatrix, V_mat: np.ndarray, U_mat: np.ndarray, repeats: int) -> float:
    """Online ``eta_adm`` check with ``B V`` precomputed: lift, gap, positive/negative parts, four scalars."""
    BV = B_mat @ V_mat
    u_r, lam_r = np.ones(V_mat.shape[1]), np.ones(U_mat.shape[1])
    c_vec = np.ones(B_mat.shape[0])

    def check() -> float:
        lam = U_mat @ lam_r
        gap = BV @ u_r - c_vec
        lam_pos = np.maximum(lam, 0.0)
        return (np.linalg.norm(np.minimum(lam, 0.0)) + np.linalg.norm(np.maximum(gap, 0.0))
                + np.sqrt(abs(lam_pos @ gap)))

    return best_of(check, repeats)
