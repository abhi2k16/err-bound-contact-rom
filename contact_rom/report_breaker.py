"""Console table and saved output for E5."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .config import Config


def print_breaker_table(table: dict, tau_info: dict, cfg: Config) -> None:
    """One row per (basis, r_lambda, indicator, query set)."""
    header = (f"{'basis':<5}{'r_lam':>6}{'indicator':>12}{'set':>8}{'tau':>10}{'n':>4}{'no-breaker unsafe':>19}"
              f"{'fallback':>10}{'missed':>8}{'max err acc.':>14}{'speedup':>9}")
    print(f"E5 circuit breaker (tau calibrated on validation, unsafe = rel. energy error > {cfg.err_silent:.0%})")
    print(header)
    print("-" * len(header))
    for (kind, r_dual, indicator, set_name), out in table.items():
        tau = tau_info[(kind, r_dual)][indicator]
        print(f"{kind:<5}{r_dual:>6}{indicator:>12}{set_name:>8}{tau:>10.2e}{out['n']:>4d}"
              f"{out['unsafe_no_breaker']:>19.2f}{out['fallback']:>10.2f}{out['missed']:>8d}"
              f"{out['max_err_accepted']:>14.2e}{out['speedup']:>9.2f}")


def print_breaker_summary(table: dict) -> None:
    """Pooled missed / total and mean fallback per indicator and query set."""
    groups: dict = {}
    for (kind, r_dual, indicator, set_name), out in table.items():
        slot = groups.setdefault((indicator, set_name), {"missed": 0, "n": 0, "fb": [], "sp": []})
        slot["missed"] += out["missed"]
        slot["n"] += out["n"]
        slot["fb"].append(out["fallback"])
        slot["sp"].append(out["speedup"])
    print("pooled over bases (missed unsafe answers / queries, mean fallback, mean speedup):")
    for (indicator, set_name), slot in groups.items():
        print(f"  {indicator:<11}{set_name:<8} missed {slot['missed']}/{slot['n']}, "
              f"fallback {np.mean(slot['fb']):.2f}, speedup {np.mean(slot['sp']):.2f}")


def save_breaker(cfg: Config, table: dict) -> Path:
    """Save the E5 table as ``.npz`` (keys ``kind_r_indicator_set_metric``)."""
    path = Path(cfg.out_dir) / cfg.run_id / "e5_circuit_breaker.npz"
    flat = {"_".join(map(str, key)) + f"_{metric}": np.array(value)
            for key, out in table.items() for metric, value in out.items()}
    np.savez(path, **flat)
    return path
