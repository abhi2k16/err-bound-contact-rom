# Reliability of reduced contact models: code, results and manuscript

This repository holds the `contact_rom` package (a posteriori admissibility estimator, five dual-basis constructions,
break-even cost model, circuit breaker), the scripts that generate the manuscript figures, and the manuscript sources.

The manuscript is `manuscript_master_v2.tex`. Section drafts and plans are the other `*.md` files;
`implementation_plan.md` documents the package layout and the acceptance criteria.

## Setup

The code needs Python 3 with `numpy`, `scipy` and `matplotlib`, plus two local packages from the computation suite
(`fea_engine` 1.0.1 for the finite-element parts and the exact contact solvers, `rom_engine` 0.2.0 for the reduced contact models).
`pyvista` is needed only for the geometry figures. The results in the manuscript were produced with Python 3.14.3,
numpy 2.5.4 and scipy 1.18.1; every run logs these versions to `results/<run>/config.json`.

```
python -m venv .venv
.venv\Scripts\pip install numpy scipy matplotlib pyvista
.venv\Scripts\pip install -e <path-to>/computation-suite/fea_engine -e <path-to>/computation-suite/rom_engine
```

LaTeX (`latexmk`, for example TinyTeX) is needed to build the PDF.

## Running the experiments

All commands are run from this directory. `python -m contact_rom` pins BLAS and OpenMP to one thread, which the
timing experiments (E4, E5) require. Outputs go to `results/<benchmark>_<experiment>/` (ignored by git).
Runtimes are rough values from one Windows machine.

| Command | What it does | Output | Runtime |
|---|---|---|---|
| `python -m contact_rom e1 --bench b1` | Silent-failure audit (signed SVD vs NMF) and estimator tables | `results/b1_e1/e1_b1.npz` | minutes |
| `python -m contact_rom all --bench b2` | E1, E3 and E4 on the membrane | `results/b2_all/` | minutes |
| `python -m contact_rom e3 --bench b1` and `--bench b2` | Five-basis comparison with calibrated thresholds | `results/b1_e3/`, `results/b2_e3/` | minutes |
| `python -m contact_rom e4 --bench b2` | Break-even model against timings, membrane size sweep | `results/b2_e4/e4_b2.npz`, `e4_sweep.npz` | minutes |
| `python -m contact_rom e4 --bench b3` | Break-even size sweep on the two-body problem, N = 24 to 16,653 | `results/b3_e4/e4_sweep_b3.npz` | about 10 to 15 minutes |
| `python -m contact_rom e5 --bench b3` | Circuit breaker on B3, in-range and extrapolated queries | `results/b3_e5/e5_circuit_breaker.npz` | long |
| `python -m contact_rom summary` | Scores the four acceptance criteria on B1 and B2 and reads the saved break-even sweeps | console | about a minute |
| `python -m contact_rom robustness` | Criteria 1 to 3 over five seeds (42 to 46) with 40 and 200 validation queries | `results/robustness/robustness.json` | about 15 to 30 minutes |

B3 can be run on smaller meshes for development with `--b3-n-xy` and `--b3-nz`. Experiment E6 (time-dependent dynamics) is
described in the protocol but is not implemented.

## Figures and manuscript

```
python make_figures.py             # figures 2, 5, 6, 8, 9 from results/ (vector PDFs in figures/)
python make_geometry_figures.py    # figures 3, 4, 7; needs pyvista; solves a few B2 and B3 queries
python compile_latex.py manuscript_master_v2.tex
```

`make_figures.py` reads these result folders: `b1_e1`, `b2_all` (or `b2_e1`), `b1_e3`, `b2_e3`, `b2_e4`, `b3_e4`, `b3_e5`.
The workflow diagram is the TikZ source `figures/fig_workflow.tex`. Intermediate renderings go to `figures/raster/` (ignored).
Close any PDF viewer before rebuilding, or the PDF cannot be overwritten.

## What the code computes

- The bound `bound1` in the audits uses the exact full-order multiplier norm and is a validation quantity (Theorem 1 holds
  for it by construction). `bound_Clam` replaces that norm by the offline constant `C_lambda`, the largest multiplier norm among the
  training snapshots, and is the version usable online; it is rigorous only while the query multiplier norm stays below `C_lambda`.
  Both are scored in `python -m contact_rom summary`.
- Acceptance criteria and their agreed scoping are described in `implementation_plan.md` and in Section 5 of the manuscript.

## Layout

```
contact_rom/    package: benchmarks, FOM adapter, bases, estimator, audit, calibration, cost model, timing, reports
figures/        generated figure PDFs and the TikZ diagram
results/        saved experiment outputs (not tracked)
*.md, *.tex     manuscript sources, section drafts, plans
```
