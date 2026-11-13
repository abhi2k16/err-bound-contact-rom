"""Experiment runner: E1 (silent-failure audit) and E3 (five-basis comparison) on B1/B2."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from dataclasses import replace

import numpy as np

from .acceptance import criterion4, print_summary, score_benchmark
from .audit import run_basis_audit_sets, run_silent_failure_audit
from .bases import Snapshots, build_snapshots
from .benchmark_3d import build_b3
from .benchmarks import ContactProblem, build_b1, build_b2
from .calibration import calibrate_all
from .circuit_breaker import breaker_table
from .query_context import cross_check_sparse_solver
from .timing import time_fom_queries
from .config import Config
from .dual_bases import DUAL_KINDS
from .break_even_experiment import (criterion4_verdict, crossovers, run_e4_entries, run_sweep,
                                    sweep_problem_b3)
from .environment import environment_info, save_run_metadata
from .estimator import tau_adm_from_snapshots
from .fom import cross_check_methods, verify_fom_on_samples
from .metrics import assert_theorem1, criterion1_verdict, criterion2_verdict
from .report_breaker import print_breaker_summary, print_breaker_table, save_breaker
from .report_break_even import print_e4_entries, print_regime_summary, print_sweep, save_e4
from .report import (print_audit_table, print_calibrated_summary, print_calibrated_table, print_comparison_table,
                     print_estimator_table, print_step1_summary, save_audit)

BENCHMARKS = {"b1": lambda cfg: build_b1(cfg.b1), "b2": lambda cfg: build_b2(cfg.b2),
              "b3": lambda cfg: build_b3(cfg.b3)}


def prepare(cfg: Config, problem: ContactProblem) -> tuple:
    """Return ``(snapshots, mu_test, mu_val)`` drawn from one seeded generator (train, test, validation)."""
    rng = np.random.default_rng(cfg.seed + 1)
    snaps = build_snapshots(problem, problem.sample_mu(rng, cfg.n_train), cfg.lcp_method)
    mu_test = problem.sample_mu(rng, cfg.n_test)
    return snaps, mu_test, problem.sample_mu(rng, cfg.n_val)


def print_criterion1(audit: dict, cfg: Config, tau_adm: float) -> None:
    """Criterion 1 (missed / silent failures, pooled over r_lambda) for Eq. (15) and eta_adm_v2."""
    print("criterion 1 missed/silent, Eq.(15) :", criterion1_verdict(audit, cfg, tau_adm, "eta_adm"))
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


def run_e4(cfg: Config, problem: ContactProblem) -> None:
    """E4: predicted vs measured speedup per basis plus the mesh-size sweep (B2); B3 runs only its sweep."""
    if problem.name == "B3":
        sweep = run_sweep(cfg, cfg.e4_grids_b3, sweep_problem_b3)
        cross = crossovers(sweep)
        print_sweep(sweep, cross, criterion4_verdict(cross), "B3")
        print("saved:", save_e4(cfg, "e4_sweep_b3", sweep))
        return
    snaps, mu_test, _ = prepare(cfg, problem)
    audits = run_basis_audit_sets(problem, cfg, snaps, {"test": mu_test}, DUAL_KINDS, cfg.r_dual_list_e3, repair=True)
    rows, fom, _ = run_e4_entries(cfg, problem, snaps, mu_test, audits["test"])
    print_e4_entries(rows, float(np.median(fom.k_fom)), float(np.median(fom.m_act)), problem.name)
    print_regime_summary(rows)
    print("saved:", save_e4(cfg, f"e4_{problem.name.lower()}", rows))
    if problem.name == "B2":
        sweep = run_sweep(cfg)
        cross = crossovers(sweep)
        print_sweep(sweep, cross, criterion4_verdict(cross))
        print("saved:", save_e4(cfg, "e4_sweep", sweep))


def override_b3(cfg: Config, n_xy: int | None, nz: int | None) -> Config:
    """Return ``cfg`` with B3 mesh sizes overridden (development runs on smaller meshes)."""
    b3 = cfg.b3
    b3 = replace(b3, n_xy=n_xy or b3.n_xy, nz_up=nz or b3.nz_up, nz_low=nz or b3.nz_low)
    return replace(cfg, b3=b3)


def e5_sets(cfg: Config, problem: ContactProblem) -> tuple:
    """Training snapshots and the query sets ``val``, ``test`` (in range) and ``extrap`` (out of range)."""
    rng = np.random.default_rng(cfg.seed + 1)
    snaps = build_snapshots(problem, problem.sample_mu(rng, cfg.e5_n_train), cfg.lcp_method)
    mu_sets = {"test": problem.sample_mu(rng, cfg.e5_n_test), "val": problem.sample_mu(rng, cfg.e5_n_val),
               "extrap": problem.sample_mu_extrap(rng, cfg.e5_n_extrap)}
    return snaps, mu_sets


def run_e5(cfg: Config, problem: ContactProblem) -> None:
    """E5: circuit breaker on a parametrized-K problem (B3) with in-range and extrapolation queries."""
    if problem.sample_mu_extrap is None:
        raise SystemExit("E5 needs a benchmark with an extrapolation sampler (use --bench b3)")
    snaps, mu_sets = e5_sets(cfg, problem)
    print(f"FOM cross-check vs fea_engine sparse solver: max |du| = {cross_check_sparse_solver(problem, mu_sets['test'][0]):.2e}")
    audits = run_basis_audit_sets(problem, cfg, snaps, mu_sets, cfg.e5_kinds, cfg.e5_r_dual_list, repair=True)
    calibrated = calibrate_all(audits["val"], audits["test"], cfg)
    fom_times = {name: float(np.median(time_fom_queries(problem, mu_sets[name][:cfg.e5_timing_queries],
                                                         cfg.e5_repeats).t_full)) for name in ("test", "extrap")}
    table = breaker_table(audits, calibrated, fom_times, cfg)
    tau_info = {key: {ind: sc["tau"] for ind, sc in per.items()} for key, per in calibrated.items()}
    print(f"FOM time per query (K varying): test {fom_times['test']:.2f} s, extrapolation {fom_times['extrap']:.2f} s")
    print_breaker_table(table, tau_info, cfg)
    print_breaker_summary(table)
    print("saved:", save_breaker(cfg, table))


def run_summary(cfg: Config) -> None:
    """Step 7: acceptance criteria on B1 and B2 (recomputed) and criterion 4 from the saved E4 sweep."""
    results = {name: score_benchmark(cfg, BENCHMARKS[name](cfg), prepare) for name in ("b1", "b2")}
    out = Path(cfg.out_dir)
    print_summary(results, {"B2 membrane": criterion4(out / "b2_e4" / "e4_sweep.npz"),
                            "B3 two-body": criterion4(out / "b3_e4" / "e4_sweep_b3.npz")})


def main() -> None:
    """CLI: ``python -m contact_rom [e1|e3|e4|e5|summary|all] [--bench b1|b2]``."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment", nargs="?", default="all", choices=["e1", "e3", "e4", "e5", "summary", "all"])
    parser.add_argument("--bench", default="b1", choices=sorted(BENCHMARKS))
    parser.add_argument("--b3-n-xy", type=int, default=None, help="override B3 elements per side")
    parser.add_argument("--b3-nz", type=int, default=None, help="override B3 element layers per body")
    args = parser.parse_args()
    cfg = override_b3(Config(run_id=f"{args.bench}_{args.experiment}"), args.b3_n_xy, args.b3_nz)
    if args.experiment == "summary":
        return run_summary(cfg)
    print("metadata:", save_run_metadata(cfg))
    print("Environment:", json.dumps(environment_info()))
    problem = BENCHMARKS[args.bench](cfg)
    if args.experiment in ("e1", "all"):
        run_e1(cfg, problem)
    if args.experiment in ("e3", "all"):
        run_e3(cfg, problem)
    if args.experiment == "e4" or (args.experiment == "all" and problem.name != "B3"):
        run_e4(cfg, problem)
    if args.experiment == "e5":
        run_e5(cfg, problem)


if __name__ == "__main__":
    main()
