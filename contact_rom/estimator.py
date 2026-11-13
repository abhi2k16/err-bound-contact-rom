"""A posteriori admissibility estimator and Theorem 1 bounds (manuscript Section 3)."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def trace_norm_B(K_lu: spla.SuperLU, B_mat: sp.spmatrix) -> float:
    """Trace constant ``||B||_{K->2} = ||B K^-1 B^T||_2^{1/2}`` (dense; small N only)."""
    compliance = B_mat @ K_lu.solve(B_mat.T.toarray())
    return float(np.sqrt(np.linalg.eigvalsh(0.5 * (compliance + compliance.T)).max()))


def dual_energy_norm(K_lu: spla.SuperLU, r_vec: np.ndarray) -> float:
    """``||r||_{K^-1} = (r^T K^-1 r)^{1/2}``."""
    return float(np.sqrt(max(r_vec @ K_lu.solve(r_vec), 0.0)))


def estimator_terms(K_mat, B_mat, f_vec, c_vec, u_lift, lam_lift, K_lu, B_norm) -> dict:
    """Computable indicators of the lifted pair ``(u~, lam~)`` (paper convention).

    Returns eta_eq, eta_sign, eta_gap, eta_comp, eta_adm (Eq. 3) and the
    residual ``r_eq_tilde`` with the unprojected multiplier.
    """
    lam_pos, lam_neg = np.maximum(lam_lift, 0.0), np.minimum(lam_lift, 0.0)
    gap_raw = B_mat @ u_lift - c_vec                      # > 0 means penetration
    res_pos = K_mat @ u_lift + B_mat.T @ lam_pos - f_vec
    res_full = K_mat @ u_lift + B_mat.T @ lam_lift - f_vec
    terms = {
        "eta_eq": dual_energy_norm(K_lu, res_pos),
        "eta_sign": float(np.linalg.norm(lam_neg)),
        "eta_gap": float(np.linalg.norm(np.maximum(gap_raw, 0.0))),
        "eta_comp": float(abs(lam_pos @ gap_raw)),
        "r_eq_tilde": dual_energy_norm(K_lu, res_full),
        "lam_pos_norm": float(np.linalg.norm(lam_pos)),
        "B_norm": B_norm,
    }
    terms["eta_adm"] = terms["eta_sign"] + terms["eta_gap"] + float(np.sqrt(terms["eta_comp"]))
    return terms


def error_bounds(terms: dict, lam_exact_norm: float, C_lambda: float) -> dict:
    """Theorem 1 bounds (1), (2) with exact ``||lambda||_2``.

    ``bound_Clam`` is the online-usable version of (1): ``||lambda||_2`` replaced by the
    offline constant ``C_lambda = max_train ||lambda(mu)||_2`` (Remark 1). It is
    rigorous whenever ``||lambda(mu)||_2 <= C_lambda`` (interpolation regime).
    """
    def bound1(lam_norm: float) -> float:
        return float(np.sqrt(terms["eta_eq"] ** 2 + 2 * terms["eta_comp"] + 2 * lam_norm * terms["eta_gap"]))

    bound2_sq = (2 * terms["r_eq_tilde"] ** 2 + 2 * terms["B_norm"] ** 2 * terms["eta_sign"] ** 2
                 + 2 * terms["eta_comp"] + 2 * lam_exact_norm * terms["eta_gap"])
    return {"bound1": bound1(lam_exact_norm), "bound2": float(np.sqrt(bound2_sq)),
            "bound_Clam": bound1(C_lambda)}


def eta_adm_v2(terms: dict, C_lambda: float) -> float:
    """Square-root-consistent admissibility indicator (energy-norm units).

    ``(2 ||B||^2 eta_sign^2 + 2 eta_comp + 2 C_lambda eta_gap)^{1/2}``: the
    interface part of bound (2) without the equilibrium residual, with the exact
    ``||lambda||_2`` replaced by the offline constant ``C_lambda`` and the trace
    constant ``||B||_{K->2}`` precomputed offline. Online cost stays
    ``O(n_c (r_u + r_lambda))``. Unlike Eq. (3) it is on the same scale as the
    energy-norm error, so it can share ``tau_adm``. Not an upper bound (no eta_eq).
    """
    value = (2.0 * terms["B_norm"] ** 2 * terms["eta_sign"] ** 2
             + 2.0 * terms["eta_comp"] + 2.0 * C_lambda * terms["eta_gap"])
    return float(np.sqrt(value))


def tau_adm_from_snapshots(U_snap: np.ndarray, K_mat: sp.spmatrix, err_silent: float) -> float:
    """Admissibility tolerance: ``err_silent`` times the mean training energy norm ``||u||_K``."""
    norms = np.sqrt(np.einsum("ij,ij->j", U_snap, K_mat @ U_snap))
    return float(err_silent * norms.mean())
