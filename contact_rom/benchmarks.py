"""Benchmark problems (paper convention: K u + B^T lam = f, B u <= c)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import scipy.sparse as sp

from .config import B1Config, B2Config


@dataclass
class ContactProblem:
    """Parametrized full-order contact problem in the paper's convention.

    Attributes
    ----------
    K_mat : callable ``mu -> (N, N)`` sparse SPD stiffness.
    f_vec : callable ``mu -> (N,)`` load vector.
    B_mat : callable ``mu -> (n_c, N)`` sparse constraint Jacobian.
    c_vec : callable ``mu -> (n_c,)`` gap offset (clearance).
    sample_mu : callable ``(rng, n) -> (n, n_p)`` parameter sampler.
    """

    name: str
    n_dof: int
    n_contact: int
    K_mat: Callable[[np.ndarray], sp.spmatrix]
    f_vec: Callable[[np.ndarray], np.ndarray]
    B_mat: Callable[[np.ndarray], sp.spmatrix]
    c_vec: Callable[[np.ndarray], np.ndarray]
    sample_mu: Callable[[np.random.Generator, int], np.ndarray]
    elem_dofs: np.ndarray | None = None   # (n_elem, dofs_per_elem) connectivity (used by RID)
    k_varies: bool = False                # True if K depends on mu (reduced operators rebuilt per query)
    sample_mu_extrap: Callable[[np.random.Generator, int], np.ndarray] | None = None  # out-of-range sampler


def chain_stiffness(n_dof: int, k_spring: float, k_ground: float) -> sp.csr_matrix:
    """Tridiagonal spring-chain stiffness with weak ground springs, (N, N)."""
    main = np.full(n_dof, 2.0 * k_spring + k_ground)
    off = np.full(n_dof - 1, -k_spring)
    return sp.diags([off, main, off], [-1, 0, 1], format="csr")


def contact_selector(n_dof: int, n_contact: int) -> sp.csr_matrix:
    """Node-selection Jacobian B, shape (n_c, N): picks evenly spaced nodes."""
    stride = n_dof // n_contact
    cols = np.arange(n_contact) * stride + stride // 2
    rows = np.arange(n_contact)
    return sp.csr_matrix((np.ones(n_contact), (rows, cols)), shape=(n_contact, n_dof))


def build_b1(cfg_b1: B1Config) -> ContactProblem:
    """Benchmark B1: N=60 chain, n_c=30 candidate contact nodes."""
    n_dof, n_c = cfg_b1.n_dof, cfg_b1.n_contact
    K_const = chain_stiffness(n_dof, cfg_b1.k_spring, cfg_b1.k_ground)
    B_const = contact_selector(n_dof, n_c)
    x_contact = np.linspace(-1.0, 1.0, n_c)

    def f_vec(mu: np.ndarray) -> np.ndarray:
        return np.full(n_dof, float(mu[0]))

    def c_vec(mu: np.ndarray) -> np.ndarray:
        return cfg_b1.clearance0 + float(mu[1]) * x_contact**2

    def sample_mu(rng: np.random.Generator, n_samples: int) -> np.ndarray:
        lo = np.array([cfg_b1.load_range[0], cfg_b1.curvature_range[0]])
        hi = np.array([cfg_b1.load_range[1], cfg_b1.curvature_range[1]])
        return lo + rng.random((n_samples, 2)) * (hi - lo)

    return ContactProblem(
        "B1", n_dof, n_c, lambda mu: K_const, f_vec, lambda mu: B_const, c_vec, sample_mu,
        elem_dofs=np.column_stack([np.arange(n_dof - 1), np.arange(1, n_dof)]),
    )


def assert_problem_wellposed(problem: ContactProblem, mu: np.ndarray, tol_sym: float) -> None:
    """Sanity checks: shapes, K symmetric and positive definite."""
    K_mat, B_mat = problem.K_mat(mu), problem.B_mat(mu)
    f_vec, c_vec = problem.f_vec(mu), problem.c_vec(mu)
    assert K_mat.shape == (problem.n_dof, problem.n_dof)
    assert B_mat.shape == (problem.n_contact, problem.n_dof)
    assert f_vec.shape == (problem.n_dof,) and c_vec.shape == (problem.n_contact,)
    asym = abs(K_mat - K_mat.T).max()
    assert asym <= tol_sym, f"K not symmetric: {asym:.2e}"
    np.linalg.cholesky(K_mat.toarray())  # raises LinAlgError unless SPD


def membrane_stiffness(n_grid: int) -> sp.csr_matrix:
    """5-point Dirichlet Laplacian on an ``n_grid x n_grid`` interior grid, (N, N), SPD."""
    tri = sp.diags([-1.0, 2.0, -1.0], [-1, 0, 1], shape=(n_grid, n_grid))
    eye = sp.identity(n_grid)
    return (sp.kron(tri, eye) + sp.kron(eye, tri)).tocsr()


def membrane_selector(n_grid: int, n_sub: int) -> tuple[sp.csr_matrix, np.ndarray]:
    """Node-selection Jacobian for the central ``n_sub x n_sub`` block; returns (B, node coords (n_c, 2))."""
    start = (n_grid - n_sub) // 2
    ij = [(i, j) for i in range(start, start + n_sub) for j in range(start, start + n_sub)]
    cols = np.array([i * n_grid + j for i, j in ij])
    B_mat = sp.csr_matrix((np.ones(len(cols)), (np.arange(len(cols)), cols)), shape=(len(cols), n_grid**2))
    coords = (np.array(ij, dtype=float) + 1.0) / (n_grid + 1.0)
    return B_mat, coords


def membrane_cells(n_grid: int) -> np.ndarray:
    """Quad-cell connectivity, shape ((n_grid-1)^2, 4)."""
    idx = np.arange(n_grid**2).reshape(n_grid, n_grid)
    return np.column_stack([idx[:-1, :-1].ravel(), idx[:-1, 1:].ravel(), idx[1:, 1:].ravel(), idx[1:, :-1].ravel()])


def build_b2(cfg_b2: B2Config) -> ContactProblem:
    """Benchmark B2: N=324 membrane, n_c=144; mu = (F0, xc, yc)."""
    n_grid = cfg_b2.n_grid
    K_const = membrane_stiffness(n_grid)
    B_const, xy_contact = membrane_selector(n_grid, cfg_b2.n_sub)
    xy_all = (np.array([(i, j) for i in range(n_grid) for j in range(n_grid)], dtype=float) + 1.0) / (n_grid + 1.0)
    bump = np.exp(-((xy_all[:, 0] - 0.5) ** 2 + (xy_all[:, 1] - 0.5) ** 2) / 0.08)

    def f_vec(mu: np.ndarray) -> np.ndarray:
        return float(mu[0]) * bump

    def c_vec(mu: np.ndarray) -> np.ndarray:
        dist2 = (xy_contact[:, 0] - mu[1]) ** 2 + (xy_contact[:, 1] - mu[2]) ** 2
        return cfg_b2.clearance0 + cfg_b2.curvature * dist2

    def sample_mu(rng: np.random.Generator, n_samples: int) -> np.ndarray:
        lo = np.array([cfg_b2.load_range[0], cfg_b2.center_range[0], cfg_b2.center_range[0]])
        hi = np.array([cfg_b2.load_range[1], cfg_b2.center_range[1], cfg_b2.center_range[1]])
        return lo + rng.random((n_samples, 3)) * (hi - lo)

    return ContactProblem("B2", n_grid**2, B_const.shape[0], lambda mu: K_const, f_vec, lambda mu: B_const,
                          c_vec, sample_mu, elem_dofs=membrane_cells(n_grid))
