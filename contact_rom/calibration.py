"""Validation-set calibration of the admissibility threshold ``tau_adm`` (per indicator, basis, r_lambda)."""

from __future__ import annotations

import numpy as np

from .config import Config
from .metrics import silent_failure_mask

INDICATORS = ("eta_adm", "eta_adm_v2")


def calibrate_tau(entry_val: dict, cfg: Config, key: str) -> float:
    """Threshold fitted on validation queries.

    ``tau = tau_margin * min(indicator over validation silent failures)``, so every
    validation silent failure is flagged with a safety margin. Returns NaN when the
    validation set contains no silent failure (nothing to calibrate against).
    """
    silent = silent_failure_mask(entry_val, cfg)
    if not silent.any():
        return float("nan")
    return float(cfg.tau_margin * entry_val[key][silent].min())


def score_on_test(entry_test: dict, cfg: Config, key: str, tau: float) -> dict:
    """Test-set recall figures for a calibrated threshold (``tau`` NaN -> not scored, ``missed = -1``)."""
    silent = silent_failure_mask(entry_test, cfg)
    good = entry_test["err_energy"] <= cfg.err_silent
    base = {"tau": tau, "n_silent": int(silent.sum()), "n_good": int(good.sum())}
    if np.isnan(tau):
        return {**base, "missed": -1, "false_alarm": float("nan")}
    flagged = entry_test[key] > tau
    false_alarm = float(flagged[good].mean()) if good.any() else float("nan")
    return {**base, "missed": int((silent & ~flagged).sum()), "false_alarm": false_alarm}


def calibrate_all(audit_val: dict, audit_test: dict, cfg: Config) -> dict:
    """``{(kind, r_lambda): {indicator: score dict}}`` with tau fitted on validation, scored on test."""
    out: dict = {}
    for entry_key, entry_val in audit_val.items():
        out[entry_key] = {}
        for key in INDICATORS:
            tau = calibrate_tau(entry_val, cfg, key)
            out[entry_key][key] = score_on_test(audit_test[entry_key], cfg, key, tau)
    return out


def pooled_summary(calibrated: dict) -> dict:
    """Per basis kind and indicator: total missed / silent (calibrated entries only), false alarms, #uncalibrated."""
    summary: dict = {}
    for (kind, _), scores in calibrated.items():
        for key, sc in scores.items():
            slot = summary.setdefault(kind, {}).setdefault(key, {"missed": 0, "silent": 0, "fa": [], "uncal": 0})
            if sc["missed"] < 0:
                slot["uncal"] += 1
                continue
            slot["missed"] += sc["missed"]
            slot["silent"] += sc["n_silent"]
            if not np.isnan(sc["false_alarm"]):
                slot["fa"].append(sc["false_alarm"])
    return summary
