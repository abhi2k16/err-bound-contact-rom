"""Acceptance-criteria metrics computed from audit results."""

from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr

from .config import Config


def silent_failure_mask(entry: dict, cfg: Config) -> np.ndarray:
    """True where reduced KKT looks converged but the lifted error is large."""
    return (entry["res_red"] <= cfg.tol_kkt_red) & (entry["err_energy"] > cfg.err_silent)


def assert_theorem1(audit: dict, rel_slack: float = 1e-8) -> None:
    """Theorem 1 must hold on every query: bounds (1), (2) >= true energy error."""
    for (kind, r_dual), entry in audit.items():
        for name in ("bound1", "bound2"):
            slack = entry[name] - entry["err_abs"] * (1.0 - rel_slack)
            assert slack.min() >= -1e-12, f"Theorem 1 ({name}) violated for {kind}, r_lam={r_dual}: {slack.min():.2e}"
        if kind == "nmf":
            assert entry["eta_sign"].max() <= 1e-12, "NMF basis must give eta_sign = 0"


def effectivity_stats(entry: dict, cfg: Config, bound_key: str = "bound1") -> dict:
    """Effectivity diagnostics for one (basis, r_lambda) entry.

    Acceptance criterion 2 is scoped to queries whose relative energy error
    exceeds ``cfg.err_silent`` (the errors the estimator must catch):
    ``in_scope`` is the share of those with theta in [1, 10], ``n_scope`` their count.
    ``in_all`` is the same share over all queries (diagnostic only). ``bound_key`` selects the bound: ``"bound1"``
    (exact ``||lambda||_2``, a validation quantity) or ``"bound_Clam"`` (offline constant ``C_lambda``, usable online).
    """
    theta = entry[bound_key] / np.maximum(entry["err_abs"], 1e-300)
    in_range = (theta >= 1) & (theta <= 10)
    scope = entry["err_energy"] > cfg.err_silent
    return {
        "theta_med": float(np.median(theta)),
        "in_all": float(in_range.mean()),
        "in_scope": float(in_range[scope].mean()) if scope.any() else float("nan"),
        "n_scope": int(scope.sum()),
        "spearman": float(spearmanr(entry[bound_key], entry["err_abs"])[0]),
        "valid": float(np.mean(entry[bound_key] >= entry["err_abs"] * (1.0 - 1e-8))),
        "clam_ok": float(np.mean(entry["bound_Clam"] >= entry["err_abs"])),
    }


def detection_stats(entry: dict, cfg: Config, tau_adm: float, key: str = "eta_adm") -> dict:
    """Recall on silent failures and false-alarm rate on accurate queries for indicator ``key``."""
    silent = silent_failure_mask(entry, cfg)
    flagged = entry[key] > tau_adm
    good = entry["err_energy"] <= cfg.err_silent
    return {"n_silent": int(silent.sum()), "missed": int((silent & ~flagged).sum()),
            "false_alarm": float(np.mean(flagged[good])) if good.any() else float("nan")}


def criterion1_verdict(audit: dict, cfg: Config, tau_adm: float, key: str) -> dict:
    """Criterion 1 per basis kind: total missed / total silent failures for indicator ``key``."""
    out = {}
    for kind in sorted({k for k, _ in audit}):
        stats = [detection_stats(e, cfg, tau_adm, key) for (k, _), e in audit.items() if k == kind]
        missed, total = sum(s["missed"] for s in stats), sum(s["n_silent"] for s in stats)
        out[kind] = f"{missed}/{total}"
    return out


def criterion2_verdict(audit: dict, cfg: Config, bound_key: str = "bound1") -> dict:
    """Pooled scoped effectivity per basis kind: share of in-scope queries with theta in [1, 10]."""
    out = {}
    for kind in ("svd", "nmf"):
        theta = np.concatenate([e[bound_key] / np.maximum(e["err_abs"], 1e-300)
                                for (k, _), e in audit.items() if k == kind])
        scope = np.concatenate([e["err_energy"] > cfg.err_silent for (k, _), e in audit.items() if k == kind])
        share = float(((theta >= 1) & (theta <= 10))[scope].mean()) if scope.any() else float("nan")
        out[kind] = {"n_scope": int(scope.sum()), "share": share, "pass": bool(share >= 0.90)}
    return out
