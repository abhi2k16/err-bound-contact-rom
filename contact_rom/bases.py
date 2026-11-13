"""Snapshots, POD primal basis and dual-basis constructions."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .benchmarks import ContactProblem
from .query_context import solve_fom_factored


@dataclass
class Snapshots:
    """Training snapshots: columns are parameter samples."""

    mu_set: np.ndarray      # (n_s, n_p)
    U_snap: np.ndarray      # (N, n_s) displacements
    Lam_snap: np.ndarray    # (n_c, n_s) multipliers, >= 0


def build_snapshots(problem: ContactProblem, mu_set: np.ndarray, method: str) -> Snapshots:
    """Solve the FOM at every row of ``mu_set`` and stack the results."""
    sols = [solve_fom_factored(problem, mu)[0] for mu in mu_set]
    assert all(sol.converged for sol in sols), "FOM LCP failed during snapshot generation"
    U_snap = np.column_stack([sol.u_vec for sol in sols])
    Lam_snap = np.column_stack([sol.lam_vec for sol in sols])
    assert Lam_snap.min() >= -1e-12, "FOM multipliers must be non-negative"
    return Snapshots(mu_set, U_snap, np.maximum(Lam_snap, 0.0))


def pod_primal_basis(U_snap: np.ndarray, r_primal: int) -> np.ndarray:
    """Orthonormal POD basis V_mat, shape (N, r_u), from displacement snapshots."""
    left_vecs, _, _ = np.linalg.svd(U_snap, full_matrices=False)
    V_mat = left_vecs[:, :r_primal]
    assert np.allclose(V_mat.T @ V_mat, np.eye(V_mat.shape[1]), atol=1e-10)
    return V_mat
