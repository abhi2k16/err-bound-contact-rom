"""Console tables and saved outputs for E4."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .config import Config


def print_e4_entries(rows: dict, fom_k: float, fom_m: float, bench: str) -> None:
    """Predicted vs measured speedup per (basis, r_lambda), both stiffness regimes."""
    header = (f"{'basis':<5}{'r_lam':>6}{'r_u':>5}{'k_rom':>7}{'m_r':>6}"
              f"{'S_pred aff':>12}{'S_meas aff':>12}{'S_pred var':>12}{'S_meas var':>12}{'meas/pred var':>15}")
    print(f"E4 {bench}: k_FOM = {fom_k:.1f}, m_act = {fom_m:.1f} (medians over test queries)")
    print(header)
    print("-" * len(header))
    for (kind, r_dual), row in rows.items():
        ratio = row["S_meas_var"] / row["S_pred_var"]
        print(f"{kind:<5}{r_dual:>6}{row['r_u']:>5}{row['k_rom']:>7.1f}{row['m_r']:>6.1f}"
              f"{row['S_pred_aff']:>12.3f}{row['S_meas_aff']:>12.3f}{row['S_pred_var']:>12.2f}"
              f"{row['S_meas_var']:>12.2f}{ratio:>15.2f}")


def print_regime_summary(rows: dict) -> None:
    """Share of models with S > 1 (predicted vs measured) and agreement of the break-even decision."""
    for regime in ("aff", "var"):
        pred = np.array([row[f"S_pred_{regime}"] > 1 for row in rows.values()])
        meas = np.array([row[f"S_meas_{regime}"] > 1 for row in rows.values()])
        print(f"regime {regime}: S>1 predicted {pred.sum()}/{len(pred)}, measured {meas.sum()}/{len(meas)}, "
              f"decision agrees on {np.mean(pred == meas):.0%}")


def print_sweep(rows: list, cross: dict, verdict: dict) -> None:
    """N-sweep table and crossover / criterion 4 summary."""
    header = (f"{'N':>6}{'n_c':>6}{'b_w':>5}{'m_act':>7}{'r_u':>5}"
              f"{'S_pred aff':>12}{'S_meas aff':>12}{'S_pred var':>12}{'S_meas var':>12}")
    print("E4 sweep (membrane, NMF r_lambda fixed, rank-repaired)")
    print(header)
    print("-" * len(header))
    for row in rows:
        print(f"{row['N']:>6}{row['n_c']:>6}{row['b_w']:>5}{row['m_act']:>7.1f}{row['r_u']:>5}"
              f"{row['S_pred_aff']:>12.3f}{row['S_meas_aff']:>12.3f}{row['S_pred_var']:>12.2f}{row['S_meas_var']:>12.2f}")
    for regime, values in cross.items():
        print(f"break-even N*  {regime}: predicted {values['N_pred']:.0f}, measured {values['N_meas']:.0f}, "
              f"ratio {values['ratio']:.2f}, criterion 4: {verdict[regime]}")


def save_e4(cfg: Config, name: str, rows: dict | list) -> Path:
    """Save E4 numbers (dict rows keyed by (kind, r) or a list of sweep rows) as ``.npz``."""
    path = Path(cfg.out_dir) / cfg.run_id / f"{name}.npz"
    if isinstance(rows, dict):
        flat = {f"{k[0]}_{k[1]}_{m}": np.array(v) for k, row in rows.items() for m, v in row.items()}
    else:
        flat = {m: np.array([row[m] for row in rows]) for m in rows[0]}
    np.savez(path, **flat)
    return path
