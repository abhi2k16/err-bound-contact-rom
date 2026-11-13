"""Seed robustness of the acceptance criteria 1-3 (B1, B2): rerun the E3 scoring for several random seeds."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from .acceptance import score_benchmark
from .config import Config


def run_seeds(cfg: Config, builders: dict, prepare, seeds: tuple, n_val_list: tuple, benches: tuple = ("b1", "b2")) -> dict:
    """``{bench: {n_val: {seed: score_benchmark(...)}}}`` for every seed and validation-set size."""
    out: dict = {}
    for bench in benches:
        for n_val in n_val_list:
            for seed in seeds:
                run_cfg = replace(cfg, seed=seed, n_val=n_val)
                print(f"robustness: {bench} n_val={n_val} seed={seed}", flush=True)
                out.setdefault(bench, {}).setdefault(n_val, {})[seed] = score_benchmark(
                    run_cfg, builders[bench](run_cfg), prepare)
    return out


def save_robustness(cfg: Config, results: dict) -> Path:
    """Write the raw per-seed scores to ``results/robustness/robustness.json``."""
    path = Path(cfg.out_dir) / "robustness" / "robustness.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(results, indent=1, default=float), encoding="utf-8")
    return path


def summarize(results: dict) -> None:
    """Per basis kind: criterion 1 misses per seed, criterion 2 shares, criterion 3 minimum Spearman across seeds."""
    for bench, per_nval in results.items():
        for n_val, per_seed in per_nval.items():
            seeds = list(per_seed)
            kinds = list(per_seed[seeds[0]]["c1"])
            print(f"\n[{bench}] validation queries = {n_val}, seeds = {seeds}")
            for indicator in ("eta_adm", "eta_adm_v2"):
                print(f"  criterion 1, {indicator}: missed/silent per seed, pooled over r_lambda")
                for kind in kinds:
                    cells = [per_seed[s]["c1"][kind][indicator] for s in seeds]
                    txt = " ".join(f"{c['missed']}/{c['silent']}" for c in cells)
                    uncal = sum(c["uncal"] for c in cells)
                    zero = sum(c["missed"] == 0 and c["silent"] > 0 for c in cells)
                    print(f"    {kind:<4} {txt:<40} zero-miss seeds {zero}/{len(seeds)}, uncalibrated entries {uncal}")
            for kind in kinds:
                mins = [per_seed[s]["c3_min"][kind] for s in seeds]
                meds = [per_seed[s]["c3_median"][kind] for s in seeds]
                c2 = [per_seed[s]["c2"][kind]["share"] for s in seeds if kind in per_seed[s]["c2"]]
                c2txt = f", criterion 2 share {np.nanmin(c2):.2f}-{np.nanmax(c2):.2f}" if c2 else ""
                print(f"  criterion 3 {kind:<4} Spearman min over entries: {np.min(mins):.2f} to {np.max(mins):.2f} "
                      f"(pass in {sum(m > 0.9 for m in mins)}/{len(seeds)} seeds), median {np.median(meds):.2f}{c2txt}")
                minc = [per_seed[s]["c3_min_Clam"][kind] for s in seeds if "c3_min_Clam" in per_seed[s]]
                if minc:
                    valid = [per_seed[s]["valid_Clam"][kind] for s in seeds]
                    c2c = [per_seed[s]["c2_Clam"][kind]["share"] for s in seeds if kind in per_seed[s]["c2_Clam"]]
                    c2ctxt = f", criterion 2 share {np.nanmin(c2c):.2f}-{np.nanmax(c2c):.2f}" if c2c else ""
                    print(f"      with C_lambda bound: min {np.min(minc):.2f} to {np.max(minc):.2f} "
                          f"(pass in {sum(m > 0.9 for m in minc)}/{len(seeds)} seeds), valid on {np.min(valid):.1%} of queries{c2ctxt}")
