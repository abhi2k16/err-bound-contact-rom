# err-bound-contact-rom

Code for the study **"Lifted Error Bounds, Admissibility Screening and Break-Even Analysis for
Lagrange-Multiplier Reduced-Order Contact Models"**.

Reduced-order models (ROMs) of frictionless unilateral contact can satisfy the reduced KKT conditions to machine
precision and still be wrong: a signed (SVD) dual basis lets the reduced multiplier go negative once it is lifted to the
full contact set. This repository contains the `contact_rom` Python package and the saved outputs behind the numerical
results of the manuscript:

- an a posteriori **error bound** for the lifted reduced solution (Theorem 1 of the manuscript) and cheap online
  **admissibility indicators** (`eta_adm`, `eta_adm,v2`) that flag silent failures without a full-order solve;
- five **dual-basis constructions** (signed SVD, NMF, CPG, gIS, RID) with rank repair, compared at matched dual
  dimension and matched accuracy;
- **threshold calibration** from validation queries (several rules), the resulting recall, false-alarm rate and the net
  speedup when rejected queries fall back to the full-order solver;
- a **break-even cost model** (operation counts) checked against repeated single-thread timings, plus mesh-size,
  extrapolation, ablation and storage studies;
- a **circuit breaker** that routes rejected queries to the full-order solver.

Negative findings are part of the results and are reported in the manuscript. Among them: the calibrated indicator does
not reach zero missed failures for every basis, NMF and CPG bounds are loose and weakly ordered, and on the membrane
benchmark screening can cost more than it saves.

## What is in this repository

The repository contains the `contact_rom` package, the saved experiment outputs in `results/` (NumPy `.npz`, pickle and
JSON files, about 14 MB) and this README. The manuscript, its figure and table generators, and working notes are not
part of the repository. Every number in the results can be regenerated with the commands below.

```
contact_rom/            package
  benchmarks.py         B1 spring chain, B2 membrane problem definitions
  benchmark_3d.py       B3 two-body 3D problem (Hex8 solids)
  fom.py                full-order contact solve and KKT checks
  bases.py, dual_bases.py, repair.py
                        snapshots, POD primal basis, the five dual bases, rank repair
  estimator.py          error bound and admissibility indicators
  audit.py, metrics.py  per-query audits (error, residual, indicators, bounds)
  calibration*.py       threshold rules and the calibration study
  cost.py, timing*.py   break-even cost model and timing studies
  break_even_experiment.py, scale_study.py, ablation.py, basis_comparison.py
                        break-even sweeps, mesh/extrapolation studies, ablation, matched-basis comparison
  circuit_breaker.py    screening with full-order fallback
  run.py                command-line entry point (python -m contact_rom ...)
results/                saved outputs, one folder per experiment
```

## Setup

Python 3 with `numpy`, `scipy` and `matplotlib` is needed, plus two packages
from the authors' computation suite that are **not on PyPI and not included here**:

- `fea_engine` (1.0.1): finite-element assembly and the exact contact solvers;
- `rom_engine` (0.2.0): reduced contact models and the active-set LCP solver.

The reported results were produced with Python 3.14.3, NumPy 2.5.4 and SciPy 1.18.1. Every run writes these versions to
`results/<run>/config.json`.

```
python -m venv .venv
.venv\Scripts\pip install numpy scipy matplotlib                      # Windows; use .venv/bin/pip elsewhere
.venv\Scripts\pip install -e <path-to>/computation-suite/fea_engine -e <path-to>/computation-suite/rom_engine
```

## Running the experiments

Run all commands from the repository root. `python -m contact_rom` pins BLAS and OpenMP to one thread, which the timing
experiments require. Outputs go to `results/<benchmark>_<experiment>/` (or a study-specific folder). Timing experiments must run alone on the
machine, one process at a time. Runtimes are rough values from one Windows machine.

| Command | What it does | Runtime |
|---|---|---|
| `python -m contact_rom e1 --bench b1` | Silent-failure audit, signed SVD vs NMF | minutes |
| `python -m contact_rom e3 --bench b1` (and `b2`) | Five-basis comparison, bound and indicator tables | minutes |
| `python -m contact_rom e4 --bench b2` (and `b3`) | Break-even model against timings, mesh-size sweep | minutes (B3: 10 to 15 min) |
| `python -m contact_rom e5 --bench b3` | Circuit breaker, in-range and extrapolated queries | long |
| `python -m contact_rom summary` | Scores the acceptance criteria on B1 and B2 | about a minute |
| `python -m contact_rom robustness` | Criteria over five seeds, 40 and 200 validation queries | 15 to 30 min |
| `python -m contact_rom calibration --bench b1` (and `b2`) | Calibration study: eight seeds, five threshold rules | long |
| `python -m contact_rom basiscomp --bench b1` (and `b2`) | Matched-dimension and matched-accuracy basis comparison | long |
| `python -m contact_rom timingstudy --bench b2` (and `b3`) | Ten repeated timings with bootstrap intervals | long |
| `python -m contact_rom ablation --bench b1` (and `b2`) | Primal-only, no repair, repair, screened | long |
| `python -m contact_rom extrap --bench b1` (and `b2`) | Out-of-range and contact-patch translation | long |
| `python -m contact_rom mesh --bench b1` (and `b2`) | Mesh-size study | long |

B3 can be run on smaller meshes for development with `--b3-n-xy` and `--b3-nz`.

## Saved results

`results/` holds the outputs the manuscript analyses were computed from: `b1_e1`, `b1_e3`, `b1_e4`, `b2_all`, `b2_e3`, `b2_e4`,
`b3_e4`, `b3_e5` (audits and break-even sweeps), `calibration` (per-seed audits and analyses), `basis_comparison`,
`timing_study`, `ablation`, `scale_study` (extrapolation and mesh-size studies), `robustness` and `memory_study`
(stored entries of the sparse factorization against the reduced dense model). Each run also writes its configuration and
library versions to `config.json`. Re-running a command overwrites the matching folder.

## Notation and conventions

Paper sign convention: `K u + B^T lambda = f`, `B u <= c`, `lambda >= 0`. A **silent failure** is a query with reduced
KKT residual below `1e-12` and relative energy-norm error above 5%. The bound `bound1` uses the exact full-order
multiplier norm and serves as a validation quantity; `bound_Clam` replaces that norm by the offline constant
`C_lambda` (the largest multiplier norm among the training snapshots) and is the version usable online, rigorous only
while the query multiplier norm stays below `C_lambda`.

## Scope

Frictionless unilateral contact, direct factorization in the full-order solver, single-thread timings on one machine.
The break-even sizes are specific to this code and machine and should be read as estimates, not as predictors of
absolute speedup. Iterative full-order solvers and time-dependent dynamics are not covered.

## License and citation

No license has been chosen yet, so all rights are reserved until one is added. A citation entry will be added when the
manuscript is published.
