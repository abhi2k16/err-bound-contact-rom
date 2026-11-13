"""Library/platform logging and run-metadata saving (reproducibility)."""

from __future__ import annotations

import json
import os
import platform
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
import scipy

import fea_engine
import rom_engine

from .config import Config


def environment_info() -> dict:
    """Library and platform versions for the run log."""
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "fea_engine": fea_engine.__version__,
        "rom_engine": rom_engine.__version__,
        "threads": {v: os.environ.get(v, "unset") for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
    }


def save_run_metadata(cfg: Config, problem=None) -> Path:
    """Write ``config.json`` (config + environment, plus the benchmark's sizes if ``problem`` is given)."""
    run_dir = Path(cfg.out_dir) / cfg.run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    meta = {"config": asdict(cfg), "environment": environment_info()}
    if problem is not None:
        meta["problem"] = {"name": problem.name, "n_dof": int(problem.n_dof), "n_contact": int(problem.n_contact),
                           "k_varies": bool(problem.k_varies)}
    path = run_dir / "config.json"
    path.write_text(json.dumps(meta, indent=2, default=str), encoding="utf-8")
    return path
