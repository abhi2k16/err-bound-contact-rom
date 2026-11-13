"""Primal rank repair by supremizer enrichment (manuscript Section 4.2, ``r_rep``)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass
class RepairResult:
    """Outcome of :func:`rank_repair`."""

    V_mat: np.ndarray       # (N, r_u + r_rep) enriched orthonormal primal basis
    r_rep: int              # number of enrichment columns added
    cond_before: float
    cond_after: float


def schur_condition(V_mat: np.ndarray, U_mat: np.ndarray, K_mat: sp.spmatrix, B_mat: sp.spmatrix) -> tuple:
    """Return ``(cond, eigvals, eigvecs)`` of ``S = B_r K_r^-1 B_r^T`` with ``B_r = U^T B V``."""
    B_red = U_mat.T @ (B_mat @ V_mat)
    K_red = V_mat.T @ (K_mat @ V_mat)
    S_mat = B_red @ np.linalg.solve(K_red, B_red.T)
    eigvals, eigvecs = np.linalg.eigh(0.5 * (S_mat + S_mat.T))
    lam_max = eigvals[-1]
    cond = np.inf if eigvals[0] <= 1e-14 * max(lam_max, 1e-300) else float(lam_max / eigvals[0])
    return cond, eigvals, eigvecs


def enrich_with_supremizer(V_mat: np.ndarray, vec_dual: np.ndarray, K_lu: spla.SuperLU,
                           B_mat: sp.spmatrix, U_mat: np.ndarray) -> np.ndarray | None:
    """Append ``K^-1 B^T U w`` (orthonormalised against ``V_mat``); None if it adds nothing."""
    sup = K_lu.solve(B_mat.T @ (U_mat @ vec_dual))
    sup = sup - V_mat @ (V_mat.T @ sup)
    sup = sup - V_mat @ (V_mat.T @ sup)             # re-orthogonalise
    norm = np.linalg.norm(sup)
    if norm < 1e-12:
        return None
    return np.column_stack([V_mat, sup / norm])


def rank_repair(V_mat: np.ndarray, U_mat: np.ndarray, K_mat: sp.spmatrix, B_mat: sp.spmatrix,
                cond_max: float) -> RepairResult:
    """Add supremizer columns until ``cond(B_r K_r^-1 B_r^T) <= cond_max`` (at most ``r_lambda`` columns).

    Each step enriches along the weakest eigen-direction of the reduced Schur
    complement, so ``r_rep <= r_lambda``.
    """
    K_lu = spla.splu(sp.csc_matrix(K_mat))
    cond, _, eigvecs = schur_condition(V_mat, U_mat, K_mat, B_mat)
    cond_before, V_cur = cond, V_mat
    for _ in range(U_mat.shape[1]):
        if cond <= cond_max:
            break
        enriched = enrich_with_supremizer(V_cur, eigvecs[:, 0], K_lu, B_mat, U_mat)
        if enriched is None:
            break
        V_cur = enriched
        cond, _, eigvecs = schur_condition(V_cur, U_mat, K_mat, B_mat)
    return RepairResult(V_cur, V_cur.shape[1] - V_mat.shape[1], cond_before, cond)
