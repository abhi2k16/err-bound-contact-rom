"""E5: online circuit breaker. Accept the reduced answer if the admissibility indicator is below
the calibrated threshold, otherwise fall back to the full-order solve."""

from __future__ import annotations

import numpy as np

from .calibration import score_on_test
from .config import Config


def flagged_queries(entry: dict, key: str, tau: float) -> np.ndarray:
    """Queries sent to the full-order solver. An uncalibrated threshold (NaN) falls back on everything."""
    if np.isnan(tau):
        return np.ones(len(entry[key]), dtype=bool)
    return entry[key] > tau


def reduced_time_per_query(entry: dict) -> np.ndarray:
    """ROM wall time per query when K varies: online projection + reduced solve + eta check."""
    return entry["t_proj"] + entry["t_online"] + entry["t_eta"]


def evaluate_breaker(entry: dict, key: str, tau: float, t_fom: float, cfg: Config) -> dict:
    """Outcome of the circuit breaker on one query set.

    ``fallback`` is the share of queries sent to the full-order solver (their answer is exact);
    ``missed`` counts accepted answers whose relative energy error exceeds ``cfg.err_silent``;
    ``speedup`` is ``t_fom`` over the mean cost, where a fallback pays the ROM attempt plus the FOM.
    """
    flagged = flagged_queries(entry, key, tau)
    accepted = ~flagged
    err = entry["err_energy"]
    cost = reduced_time_per_query(entry) + flagged * t_fom
    return {
        "n": int(len(err)),
        "fallback": float(flagged.mean()),
        "unsafe_no_breaker": float(np.mean(err > cfg.err_silent)),
        "missed": int((accepted & (err > cfg.err_silent)).sum()),
        "max_err_accepted": float(err[accepted].max()) if accepted.any() else float("nan"),
        "speedup": float(t_fom / cost.mean()),
    }


def breaker_table(audits: dict, calibrated: dict, fom_times: dict, cfg: Config) -> dict:
    """``{(kind, r, indicator, set_name): outcome}`` using tau calibrated on the validation set."""
    out: dict = {}
    for entry_key, per_indicator in calibrated.items():
        for indicator, scored in per_indicator.items():
            for set_name in ("test", "extrap"):
                out[(*entry_key, indicator, set_name)] = evaluate_breaker(
                    audits[set_name][entry_key], indicator, scored["tau"], fom_times[set_name], cfg)
    return out


def extrapolation_scores(audits: dict, calibrated: dict, cfg: Config) -> dict:
    """Detection scoring (missed / false alarm) of the calibrated thresholds on the extrapolation set."""
    return {entry_key: {ind: score_on_test(audits["extrap"][entry_key], cfg, ind, sc["tau"])
                        for ind, sc in per_ind.items()}
            for entry_key, per_ind in calibrated.items()}
