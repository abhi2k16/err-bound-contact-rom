"""Experiment runner: E1 (silent-failure audit) and E3 (five-basis comparison) on B1/B2."""

from __future__ import annotations

import argparse
import json

import numpy as np

from .audit import run_basis_audit_sets, run_silent_failure_audit
from .bases import Snapshots, build_snapshots
from .benchmarks import ContactProblem, build_b1, build_b2
from .calibration import calibrate_all
from .config import Config
from .dual_bases import DUAL_KINDS
from .environment import environment_info, save_run_metadata
from .estimator import tau_adm_from_snapshots
from .fom import cross_check_methods, verify_fom_on_samples
from .metrics import assert_theorem1, criterion1_verdict, criterion2_verdict
from .report import (print_audit_table, print_calibrated_summary, print_calibrated_table, print_comparison_table,
                     print_estimator_table, print_step1_summary, save_audit)

BENCHMARKS = {"b1": lambda cfg: build_b1(cfg.b1), "b2": lambda cfg: build_b2(cfg.b2)}


def prepare(cfg: Config, problem: ContactProblem) -> tuple:
    """Return ``(snapshots, mu_test, mu_val)`` drawn from one seeded generator (train, test, validation)."""
    rng = np.random.default_rng(cfg.seed + 1)
    snaps = build_snapshots(problem, problem.sample_mu(rng, cfg.n_train), cfg.lcp_method)
    mu_test = problem.sample_mu(rng, cfg.n_test)
    return snaps, mu_test, problem.sample_mu(rng, cfg.n_val)


def print_criterion1(audit: dict, cfg: Config, tau_adm: float) -> None:
    """Criterion 1 (missed / silent failures, pooled over r_lambda) for Eq. (3) and eta_adm_v2."""
    print("criterion 1 missed/silent, Eq.(3) :", criterion1_verdict(audit, cfg, tau_adm, "eta_adm"))
    print("criterion 1 missed/silent, v2     :", criterion1_verdict(audit, cfg, tau_adm, "eta_adm_v2"))


def run_e1(cfg: Config, problem: ContactProblem) -> None:
    """E1 plus step-1 FOM verification and the estimator audit (signed SVD vs NMF, no rank repair)."""
    out = verify_fom_on_samples(problem, cfg)
    print_step1_summary(problem, out, cross_check_methods(problem, out["mu_set"][0]))
    snaps, mu_test, _ = prepare(cfg, problem)
    audit = run_silent_failure_audit(problem, cfg, snaps, mu_test)
    print_audit_table(audit, cfg)
    assert_theorem1(audit)
    tau_adm = tau_adm_from_snapshots(snaps.U_snap, problem.K_mat(snaps.mu_set[0]), cfg.err_silent)
    print_estimator_table(audit, cfg, tau_adm)
    print("criterion 2 (scoped, pooled over r_lambda):", criterion2_verdict(audit, cfg))
    print_criterion1(audit, cfg, tau_adm)
    print("saved:", save_audit(audit, cfg, f"e1_{problem.name.lower()}"))


def run_e3(cfg: Config, problem: ContactProblem) -> None:
    """E3: five dual bases with rank repair on identical snapshots."""
    snaps, mu_test, mu_val = prepare(cfg, problem)
    audits = run_basis_audit_sets(problem, cfg, snaps, {"test": mu_test, "val": mu_val}, DUAL_KINDS,
                                  cfg.r_dual_list_e3, repair=True)
    audit = audits["test"]
    assert_theorem1(audit)
    tau_adm = tau_adm_from_snapshots(snaps.U_snap, problem.K_mat(snaps.mu_set[0]), cfg.err_silent)
    print_comparison_table(audit, cfg, tau_adm)
    print_criterion1(audit, cfg, tau_adm)
    calibrated = calibrate_all(audits["val"], audit, cfg)
    print_calibrated_table(calibrated, cfg)
    print_calibrated_summary(calibrated)
    print("saved:", save_audit(audit, cfg, f"e3_{problem.name.lower()}"))


def main() -> None:
    """CLI: ``python -m contact_rom [e1|e3|all] [--bench b1|b2]``."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment", nargs="?", default="all", choices=["e1", "e3", "all"])
    parser.add_argument("--bench", default="b1", choices=sorted(BENCHMARKS))
    args = parser.parse_args()
    cfg = Config(run_id=f"{args.bench}_{args.experiment}")
    print("metadata:", save_run_metadata(cfg))
    print("Environment:", json.dumps(environment_info()))
    problem = BENCHMARKS[args.bench](cfg)
    if args.experiment in ("e1", "all"):
        run_e1(cfg, problem)
    if args.experiment in ("e3", "all"):
        run_e3(cfg, problem)


if __name__ == "__main__":
    main()
