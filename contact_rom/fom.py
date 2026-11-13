"""Full-order contact solve (adapter over fea_engine) and KKT verification."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from fea_engine.contact_solvers import solve_contact_static

from .benchmarks import ContactProblem, assert_problem_wellposed
from .config import Config


@dataclass
class FomSolution:
    """Full-order contact solution in the paper's convention."""

    u_vec: np.ndarray       # (N,) displacements
    lam_vec: np.ndarray     # (n_c,) multipliers, >= 0
    gap_vec: np.ndarray     # (n_c,) c - B u, >= 0
    n_active: int
    converged: bool


def solve_fom(problem: ContactProblem, mu: np.ndarray, method: str) -> FomSolution:
    """Solve ``K u + B^T lam = f, B u <= c`` via the suite (signs flipped)."""
    K_mat, B_mat = problem.K_mat(mu), problem.B_mat(mu)
    f_vec, c_vec = problem.f_vec(mu), problem.c_vec(mu)
    result = solve_contact_static(K_mat, f_vec, -B_mat, -c_vec, method=method)
    gap_vec = c_vec - B_mat @ result.u
    return FomSolution(result.u, result.lam, gap_vec, int(result.active.size), bool(result.converged))


def kkt_residuals(problem: ContactProblem, mu: np.ndarray, sol: FomSolution) -> dict:
    """Return the four KKT residuals (each should be ~ 0 or >= 0)."""
    K_mat, B_mat = problem.K_mat(mu), problem.B_mat(mu)
    stationarity = K_mat @ sol.u_vec + B_mat.T @ sol.lam_vec - problem.f_vec(mu)
    return {
        "stationarity": float(np.linalg.norm(stationarity, np.inf)),
        "gap_violation": float(max(0.0, -sol.gap_vec.min())),
        "lam_violation": float(max(0.0, -sol.lam_vec.min())),
        "complementarity": float(np.max(np.abs(sol.lam_vec * sol.gap_vec))),
    }


def assert_fom_valid(residuals: dict, tol: float) -> None:
    """Raise if any KKT residual exceeds ``tol``."""
    for name, value in residuals.items():
        assert value <= tol, f"FOM KKT check failed: {name} = {value:.2e} > {tol:.1e}"


def verify_fom_on_samples(problem: ContactProblem, cfg: Config) -> dict:
    """Solve B1 at random parameters; verify KKT, count active constraints."""
    rng = np.random.default_rng(cfg.seed)
    mu_set = problem.sample_mu(rng, cfg.n_train)
    assert_problem_wellposed(problem, mu_set[0], cfg.tol_sym)
    worst = {"stationarity": 0.0, "gap_violation": 0.0, "lam_violation": 0.0, "complementarity": 0.0}
    n_active = np.zeros(len(mu_set), dtype=int)
    for idx, mu in enumerate(mu_set):
        sol = solve_fom(problem, mu, cfg.lcp_method)
        assert sol.converged, f"LCP did not converge at sample {idx}"
        res = kkt_residuals(problem, mu, sol)
        assert_fom_valid(res, cfg.tol_kkt_fom)
        worst = {k: max(worst[k], res[k]) for k in worst}
        n_active[idx] = sol.n_active
    return {"worst": worst, "n_active": n_active, "mu_set": mu_set}


def cross_check_methods(problem: ContactProblem, mu: np.ndarray) -> float:
    """Max |u_active_set - u_lemke| at one parameter (independent LCP solvers)."""
    sol_a = solve_fom(problem, mu, "active_set")
    sol_b = solve_fom(problem, mu, "lemke")
    return float(np.max(np.abs(sol_a.u_vec - sol_b.u_vec)))
