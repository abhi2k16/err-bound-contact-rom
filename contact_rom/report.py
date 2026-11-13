"""Console tables and saved outputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .benchmarks import ContactProblem
from .config import Config
from .calibration import pooled_summary
from .metrics import detection_stats, effectivity_stats, silent_failure_mask


def print_step1_summary(problem: ContactProblem, out: dict, cross_err: float) -> None:
    """Clean console table of step-1 metrics."""
    n_act = out["n_active"]
    rows = [
        ("benchmark", f"{problem.name}  (N={problem.n_dof}, n_c={problem.n_contact})"),
        ("samples solved", f"{len(n_act)}"),
        ("active constraints  min/mean/max", f"{n_act.min()} / {n_act.mean():.1f} / {n_act.max()}"),
        ("max |stationarity|", f"{out['worst']['stationarity']:.2e}"),
        ("max gap violation", f"{out['worst']['gap_violation']:.2e}"),
        ("max lam violation", f"{out['worst']['lam_violation']:.2e}"),
        ("max complementarity |lam*gap|", f"{out['worst']['complementarity']:.2e}"),
        ("active_set vs lemke  max |du|", f"{cross_err:.2e}"),
    ]
    width = max(len(name) for name, _ in rows)
    print("\n" + "=" * (width + 40))
    for name, value in rows:
        print(f"{name:<{width}} : {value}")
    print("=" * (width + 40))


def print_audit_table(audit: dict, cfg: Config) -> None:
    """One row per (basis, r_lambda): worst/median metrics and silent-failure count."""
    header = (f"{'basis':<5}{'r_lam':>6}{'max R_red':>11}{'min lam~':>11}{'max gap viol':>14}"
              f"{'med err_E':>11}{'max err_E':>11}{'silent fail':>13}")
    print("\n" + header + "\n" + "-" * len(header))
    for (kind, r_dual), entry in audit.items():
        n_silent = int(silent_failure_mask(entry, cfg).sum())
        print(f"{kind:<5}{r_dual:>6}{entry['res_red'].max():>11.1e}{entry['lam_min'].min():>11.2e}"
              f"{entry['gap_viol'].max():>14.2e}{np.median(entry['err_energy']):>11.2e}"
              f"{entry['err_energy'].max():>11.2e}{n_silent:>7d}/{len(entry['res_red']):<5d}")


def print_estimator_table(audit: dict, cfg: Config, tau_adm: float) -> None:
    """E1/E2 summary: Theorem 1, effectivity (scoped + all), correlation, detection."""
    header = (f"{'basis':<5}{'r_lam':>6}{'thm1':>6}{'med theta':>11}{'in[1,10] scope':>16}{'(n)':>6}"
              f"{'in[1,10] all':>14}{'spearman':>10}{'Clam ok':>9}{'med eta_adm':>13}{'missed':>8}{'FA':>6}{'med v2':>10}{'missed v2':>11}{'FA v2':>7}")
    print(f"\nEstimator audit (tau_adm = {tau_adm:.3e}; scope = rel. energy error > {cfg.err_silent:.0%})")
    print(header + "\n" + "-" * len(header))
    for (kind, r_dual), entry in audit.items():
        eff, det = effectivity_stats(entry, cfg), detection_stats(entry, cfg, tau_adm)
        det2 = detection_stats(entry, cfg, tau_adm, "eta_adm_v2")
        scope_txt = "n/a" if eff["n_scope"] == 0 else f"{eff['in_scope']:.2f}"
        print(f"{kind:<5}{r_dual:>6}{'ok':>6}{eff['theta_med']:>11.2f}{scope_txt:>16}{eff['n_scope']:>6d}"
              f"{eff['in_all']:>14.2f}{eff['spearman']:>10.3f}{eff['clam_ok']:>9.2f}"
              f"{np.median(entry['eta_adm']):>13.2e}{det['missed']:>4d}/{det['n_silent']:<3d}{det['false_alarm']:>6.2f}"
              f"{np.median(entry['eta_adm_v2']):>10.2e}{det2['missed']:>7d}/{det2['n_silent']:<3d}{det2['false_alarm']:>7.2f}")


def save_audit(audit: dict, cfg: Config, name: str) -> Path:
    """Flatten the audit dict into one ``.npz`` named ``name`` (keys ``kind_rlam_metric``)."""
    flat = {f"{kind}_{r_dual}_{key}": vals for (kind, r_dual), entry in audit.items() for key, vals in entry.items()}
    path = Path(cfg.out_dir) / cfg.run_id / f"{name}.npz"
    np.savez(path, **flat)
    return path


def print_comparison_table(audit: dict, cfg: Config, tau_adm: float) -> None:
    """E3 Pareto table: accuracy, admissibility, rank repair, iterations, online and offline time."""
    header = (f"{'basis':<5}{'r_lam':>6}{'dim':>5}{'r_rep':>6}{'med err_E':>11}{'max err_E':>11}{'med err_lam':>13}"
              f"{'med eta_adm':>13}{'missed':>8}{'missed v2':>11}{'k_rom':>7}{'t_on ms':>9}{'t_off s':>9}{'spearman':>10}")
    print(f"\nE3 basis comparison (tau_adm = {tau_adm:.3e})\n" + header + "\n" + "-" * len(header))
    for (kind, r_dual), entry in audit.items():
        det = detection_stats(entry, cfg, tau_adm)
        det2 = detection_stats(entry, cfg, tau_adm, "eta_adm_v2")
        rho = effectivity_stats(entry, cfg)["spearman"]
        print(f"{kind:<5}{r_dual:>6}{int(entry['r_dual_eff'][0]):>5}{int(entry['r_rep'][0]):>6}"
              f"{np.median(entry['err_energy']):>11.2e}{entry['err_energy'].max():>11.2e}"
              f"{np.median(entry['err_lam']):>13.2e}{np.median(entry['eta_adm']):>13.2e}"
              f"{det['missed']:>4d}/{det['n_silent']:<3d}{det2['missed']:>7d}/{det2['n_silent']:<3d}{np.median(entry['k_rom']):>7.1f}"
              f"{1e3 * np.median(entry['t_online']):>9.2f}{entry['t_offline'][0]:>9.2f}{rho:>10.3f}")


def print_calibrated_table(calibrated: dict, cfg: Config) -> None:
    """Per (basis, r_lambda): validation-calibrated tau and test-set misses / false alarms, Eq. (15) vs v2."""
    header = (f"{'basis':<5}{'r_lam':>6}{'silent':>8}{'good':>6}{'tau Eq3':>10}{'miss':>6}{'FA':>7}"
              f"{'tau v2':>10}{'miss':>6}{'FA':>7}")
    print(f"Calibrated tau_adm (fit on {cfg.n_val} validation queries, scored on test; margin {cfg.tau_margin})")
    print(header)
    print("-" * len(header))
    for (kind, r_dual), scores in calibrated.items():
        a, b = scores["eta_adm"], scores["eta_adm_v2"]
        miss_a = "   n/a" if a["missed"] < 0 else f"{a['missed']:>6d}"
        miss_b = "   n/a" if b["missed"] < 0 else f"{b['missed']:>6d}"
        print(f"{kind:<5}{r_dual:>6}{a['n_silent']:>8d}{a['n_good']:>6d}{a['tau']:>10.2e}{miss_a}"
              f"{a['false_alarm']:>7.2f}{b['tau']:>10.2e}{miss_b}{b['false_alarm']:>7.2f}")


def print_calibrated_summary(calibrated: dict) -> None:
    """Pooled test misses and mean false alarm per basis kind (calibrated entries only)."""
    print("pooled over calibrated r_lambda (missed/silent, mean FA, #uncalibrated):")
    for kind, per_key in pooled_summary(calibrated).items():
        parts = []
        for key, slot in per_key.items():
            fa = float(np.mean(slot["fa"])) if slot["fa"] else float("nan")
            parts.append(f"{key}: {slot['missed']}/{slot['silent']}, FA {fa:.2f}, uncal {slot['uncal']}")
        print(f"  {kind:<4}" + " | ".join(parts))
