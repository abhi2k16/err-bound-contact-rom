"""Repeat the wall-clock timing study several times and extend the membrane sweep to smaller meshes.

Each repetition runs ``contact_rom.timing_study.measure_size`` for every size (the full-order and the reduced solves of
one query are timed back to back, so both regimes see the same machine state) and stores its own JSON in
``results/timing_repeats/``. ``summarize`` then reports, per regime, the break-even size ``N*`` of every repetition and
the spread of the measured speedups between repetitions.

Usage (from the repository root, nothing else running):
    python -m scripts.analysis.timing_repeats run b2 3      # three repetitions of the membrane sweep
    python -m scripts.analysis.timing_repeats run b3 2
    python -m scripts.analysis.timing_repeats summarize
"""

import json
import sys
from pathlib import Path

import numpy as np

from contact_rom.config import Config
from contact_rom import timing_study as ts
from contact_rom.break_even_experiment import sweep_problem, sweep_problem_b3, sweep_problem_b3_code

OUT = Path("results/timing_repeats")
SMALL_B2 = (4,)               # one additional small membrane grid (N = 16); grids 2 and 3 are too small for the contact sub-grid


def jobs(cfg: Config, bench: str) -> list:
    if bench == "b2":
        return [(g, sweep_problem(cfg, g), g) for g in SMALL_B2 + tuple(cfg.e4_grids)]
    out = [(g, sweep_problem_b3(cfg, g), g) for g in cfg.e4_grids_b3]
    return out + [(c, sweep_problem_b3_code(cfg, c), c) for c in cfg.e4_codes_b3_small]


def run(bench: str, n_rep: int, repeats: int = 10) -> None:
    cfg = Config()
    OUT.mkdir(parents=True, exist_ok=True)
    for rep in range(n_rep):
        path = OUT / f"timing_{bench}_rep{rep}.json"
        if path.exists():
            print("cached", path, flush=True)
            continue
        rows = []
        for key, problem, offset in jobs(cfg, bench):
            print(f"{bench} rep {rep}: size parameter {key}", flush=True)
            try:
                row = ts.measure_size(cfg, problem, offset, repeats)
            except Exception as exc:                      # tiny meshes may have no contact
                print("  skipped:", type(exc).__name__, exc, flush=True)
                row = None
            if row is not None:
                rows.append(row)
        rows.sort(key=lambda r: r["N"])
        path.write_text(json.dumps({"bench": bench, "repeats": repeats, "n_queries": cfg.e4_n_test, "rows": rows}), encoding="utf-8")
        print("saved", path, flush=True)


def speedups(data: dict) -> dict:
    out = {"N": [], "aff": [], "var": []}
    for r in data["rows"]:
        a = {k: np.array(r[k]) for k in ("t_aff", "t_full", "t_rom", "t_proj", "t_eta")}
        a["t_rom_aff"] = np.array(r.get("t_rom_aff", r["t_rom"]))
        s = ts.speedup_from(a["t_aff"], a["t_full"], a["t_rom"], a["t_proj"], a["t_eta"], a["t_rom_aff"])
        out["N"].append(r["N"]); out["aff"].append(float(s["aff"])); out["var"].append(float(s["var"]))
    return out


def crossing(n: list, s: list) -> float:
    """First log-log interpolated crossing of S = 1 from below; NaN if S > 1 at the smallest size or never crosses."""
    n, s = np.array(n, float), np.array(s, float)
    if s[0] >= 1.0:
        return float("nan")
    for i in range(1, len(n)):
        if s[i - 1] < 1.0 <= s[i]:
            f = (0 - np.log(s[i - 1])) / (np.log(s[i]) - np.log(s[i - 1]))
            return float(np.exp(np.log(n[i - 1]) + f * (np.log(n[i]) - np.log(n[i - 1]))))
    return float("nan")


def summarize() -> None:
    res = {}
    for bench in ("b2", "b3"):
        files = sorted(OUT.glob(f"timing_{bench}_rep*.json"))
        if not files:
            continue
        sp = [speedups(json.loads(f.read_text())) for f in files]
        res[bench] = {"reps": len(files), "crossing_aff": [crossing(s["N"], s["aff"]) for s in sp],
                      "crossing_var": [crossing(s["N"], s["var"]) for s in sp],
                      "smallest_S_var": [s["var"][0] for s in sp], "smallest_N": [s["N"][0] for s in sp]}
        common = sorted(set.intersection(*[set(s["N"]) for s in sp]))
        spread = {}
        for reg in ("aff", "var"):
            vals = np.array([[dict(zip(s["N"], s[reg]))[n] for n in common] for s in sp])      # reps x sizes
            spread[reg] = {"max_ratio_between_reps": float(np.max(vals.max(axis=0) / vals.min(axis=0))),
                           "median_ratio_between_reps": float(np.median(vals.max(axis=0) / vals.min(axis=0)))}
        res[bench]["spread"] = spread
        res[bench]["speedups"] = [{k: s[k] for k in s} for s in sp]
    Path("results/analysis").mkdir(parents=True, exist_ok=True)
    Path("results/analysis/timing_repeats.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    for bench, r in res.items():
        print(bench, "repetitions", r["reps"])
        print("  N* constant K per repetition:", [round(x, 1) for x in r["crossing_aff"]])
        print("  N* varying K per repetition :", [round(x, 1) for x in r["crossing_var"]])
        print("  smallest N and varying-K speedup there:", r["smallest_N"], [round(x, 2) for x in r["smallest_S_var"]])
        print("  speedup ratio between repetitions (max over sizes / median):", {k: (round(v["max_ratio_between_reps"], 2), round(v["median_ratio_between_reps"], 2)) for k, v in r["spread"].items()})


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], int(sys.argv[3]))
    else:
        summarize()
