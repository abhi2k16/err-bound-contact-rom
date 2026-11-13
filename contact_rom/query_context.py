"""Per-query full-order data: one factorisation shared by the exact solve and every estimator."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from rom_engine.complementarity import lcp_active_set

from .benchmarks import ContactProblem
from .fom import FomSolution
from .timing import assemble_factors


@dataclass
class QueryContext:
    """Everything about one parameter query that does not depend on the reduced model."""

    mu: np.ndarray
    K_mat: sp.spmatrix
    B_mat: sp.spmatrix
    f_vec: np.ndarray
    c_vec: np.ndarray
    K_lu: spla.SuperLU
    B_norm: float               # trace constant ||B||_{K->2} = sqrt(lambda_max(B K^-1 B^T))
    fom: FomSolution
    fom_iterations: int


def solve_fom_factored(problem: ContactProblem, mu: np.ndarray) -> tuple:
    """Exact contact solve via one sparse LU of ``K(mu)``; returns ``(FomSolution, lu, compliance, iterations)``.

    Same LCP as ``fea_engine.contact_solvers.solve_contact_static`` (paper convention flipped to the
    suite's), but the factorisation and the compliance ``A = B K^-1 B^T`` are kept for reuse.
    """
    K_mat, B_mat = problem.K_mat(mu), problem.B_mat(mu)
    f_vec, c_vec = problem.f_vec(mu), problem.c_vec(mu)
    lu, Kinv_Bt, A_lcp = assemble_factors(K_mat, -B_mat)
    Kinv_f = lu.solve(f_vec)
    lam, ok, iters = lcp_active_set(A_lcp, -(B_mat @ Kinv_f) + c_vec)
    u_vec = Kinv_f + Kinv_Bt @ lam
    sol = FomSolution(u_vec, lam, c_vec - B_mat @ u_vec, int(np.sum(lam > 1e-12)), bool(ok))
    return sol, lu, A_lcp, int(iters)


def build_query_context(problem: ContactProblem, mu: np.ndarray) -> QueryContext:
    """Solve the FOM at ``mu`` and collect the shared factorisation and constants."""
    sol, lu, A_lcp, iters = solve_fom_factored(problem, mu)
    assert sol.converged, "FOM LCP did not converge"
    B_norm = float(np.sqrt(max(np.linalg.eigvalsh(A_lcp)[-1], 0.0)))
    return QueryContext(mu, problem.K_mat(mu), problem.B_mat(mu), problem.f_vec(mu), problem.c_vec(mu),
                        lu, B_norm, sol, iters)


def cross_check_sparse_solver(problem: ContactProblem, mu: np.ndarray, tol: float = 1e-8) -> float:
    """Max ``|u_ours - u_suite|`` against ``fea_engine``'s sparse primal-dual contact solver (asserts ``<= tol``)."""
    from fea_engine.contact_solvers import solve_contact_static_sparse

    ours = solve_fom_factored(problem, mu)[0]
    theirs = solve_contact_static_sparse(problem.K_mat(mu), problem.f_vec(mu), -problem.B_mat(mu), -problem.c_vec(mu))
    assert theirs.converged, f"suite sparse solver did not converge: {theirs.message}"
    diff = float(np.max(np.abs(ours.u_vec - theirs.u)))
    assert diff <= tol * max(1.0, float(np.max(np.abs(ours.u_vec)))), f"FOM mismatch vs suite sparse solver: {diff:.2e}"
    return diff
