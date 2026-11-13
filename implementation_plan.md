# Plan: `contact_rom/` package — Implementation of the Proposed Research

Spec: engineering rules (modular, reproducible, verified) applied to the v2 sections 2–6 of the manuscript.

Built on the local **computation-suite** (`C:\Users\abhij\Desktop\Physical_ML\computation-suite`, repo `abhi2k16/fea-rom-computation-suite`): `fea_engine` 1.0.1 (full-order FE and exact contact solvers) and `rom_engine` 0.2.0 (reduction, contact ROMs). Their only dependencies are `numpy` and `scipy`; the script adds `matplotlib` for plots. No FEniCSx, GetFEM or PETSc.

## 1. Setup

- The project `.venv` is empty and neither package is importable yet. Install both in editable mode: `pip install -e <suite>/fea_engine -e <suite>/rom_engine` (plus numpy, scipy, matplotlib).
- Log `fea_engine.__version__` and `rom_engine.__version__` alongside python/numpy/scipy at runtime.
- `rom_engine` never imports `fea_engine`; the script is the only place the two meet, passing plain arrays `K, B, f, c`.

## 2. Reuse map: suite vs. new code

| Need (paper section) | From computation-suite | New code in the script |
|---|---|---|
| Full-order contact solve, FOM truth (§2.1) | `fea_engine.contact_solvers.solve_contact_static` / `solve_contact_static_sparse` (dense or sparse reduced/saddle routes, `ContactResult`), `contact_complementarity_residual` | thin wrapper that returns `(u, lambda_vec)` and asserts KKT ≤ 1e-12 |
| Mesh, K assembly, FEM benchmarks B2/B3 | `fea_engine` mesh, elements (solids, plates), `FESystem`, membrane obstacle problems from the gallery (361–3,481 unknowns, 39,601 sparse) | B1 spring chain builder; parameter samplers μ; B3 two-body assembly from suite solids |
| Reduced LCP solve (§2.2) | `rom_engine.complementarity.solve_reduced_contact`, `lcp_active_set` (warm start `active0`), `lcp_lemke` | iteration-count capture for k_ROM |
| Dual bases: SVD, NMF (§4.1) | `contact_dual_basis.svd_dual_basis`, `nonneg_dual_basis` / `NMFBasis`, `ContactROM`, `ContactROMResult` | rank-repair logic (`r_rep`, cond check ≤ 1e6) around SVD |
| Dual basis: cone-projected greedy | `contact_dual_basis.sparse_dictionary_fit` (NNLS) and `contact_greedy.greedy_contact_sampling` as reference | CPG selection loop per Benaceur et al. |
| Dual basis: gIS | none | greedy index selection (columns of I) |
| Dual basis: RID (§4.1) | `contact_rid.HybridHRContact`, `deim_rid_indices`, `rid_elements`, `contact_rows_in_rid` | adapter exposing RID under the common basis interface |
| Dynamics for E6 (optional) | `ContactROM.solve_dynamic(... scheme="bdf2")`, `solve_contact_dynamics` | η_adm(t_n) logging per step |
| Error metrics, 5-way comparison | `contact_validation.contact_errors`, `max_normalised_error`, `contact_comparison.compare_contact_models` | θ_eff, Spearman, recall metrics |
| **A posteriori estimator (§3)** | `contact_greedy.contact_condition_indicator` is related (violation of the three Hertz–Signorini–Moreau conditions) but is not Theorem 1 | **all new:** η_eq, η_sign, η_gap, η_comp, η_bound, η_adm |
| **Cost model (§4.2)** | none | **all new:** `flops_fom`, `flops_rom`, S_pred; measured FLOP and wall-clock timers |

Rule: wherever the suite provides a function, call it rather than reimplement it. Anything the paper contributes (estimator, rank-repair audit, cost model, acceptance checks) is written fresh.

## 3. Scope

| Item | Plan |
|---|---|
| B1 spring chain (N=60, n_c=30) | Built in the script (not in the suite); exactly as specified. |
| B2 membrane (18×18, N=324, n_c=144) | Build from suite mesh/assembly; matches the suite's membrane obstacle gallery problem. |
| B3 3D two-body | Hex8 element stiffness from the suite (`Hex8Solid3D`), structured grid assembled in scipy (all cells identical); FOM through own factored solve, verified against `solve_contact_static_sparse`. N ≈ 2·10⁴ default (configurable to 10⁵). Parametrized K(μ) via stiffness contrast. This replaces the earlier plain-scipy substitute. |
| E1–E5 | Core experiments, run by default. |
| E6 (BDF2 dynamics) | Optional, off by default (`CONFIG.run_e6 = False`). |

## 4. Package layout (`contact_rom/`, run with `python -m contact_rom`)

Split into small single-purpose modules (decided after step 3; originally planned as one file). Steps 1-3 are implemented; later steps add modules rather than growing existing ones.

| Module | Responsibility | Status |
|---|---|---|
| `config.py` | Frozen `Config` / `B1Config`: all hyperparameters, tolerances, paths | done |
| `environment.py` | Library/platform logging, `config.json` metadata | done |
| `benchmarks.py` | `ContactProblem`, B1 and B2 builders, well-posedness assertions (B3 pending, step 6) | B1, B2 done |
| `fom.py` | Full-order solve adapter over `fea_engine`, KKT checks, step-1 verification | done |
| `bases.py` | Snapshots, POD primal basis | done |
| `dual_bases.py` | SVD, NMF (suite), CPG, gIS, RID (suite DEIM/RID helpers) behind `make_dual_basis` | done (step 4) |
| `repair.py` | Supremizer rank repair to cond(B_r K_r^-1 B_r^T) <= 1e6, gives `r_rep` | done (step 4) |
| `estimator.py` | eta_eq/sign/gap/comp, eta_adm, Theorem 1 bounds, C_lambda bound, tau_adm | done |
| `calibration.py` | Validation-set calibration of tau_adm per (indicator, basis, r_lambda); test-set misses and false alarms | done (step 4b) |
| `audit.py` | Reduced solve via suite `ContactROM`, per-query metrics, E1/E2 in-sample audit | done |
| `metrics.py` | Acceptance-criteria metrics (effectivity, Spearman, detection, criterion 2) | done |
| `report.py` | Console tables and saved `.npz` | done |
| `cost.py`, `timing.py`, `break_even_experiment.py`, `report_break_even.py` | Eqs. (16)-(17) FLOP model, FOM/ROM wall-clock timers, E4 table + N sweep + crossover (criterion 4) | done (step 5) |
| `benchmark_3d.py` | B3: quarter-symmetric two-body hex8 contact, parametrized K(mu) = ratio*K_up + K_low, element stiffness from `fea_engine` `Hex8Solid3D` | done (step 6) |
| `query_context.py` | One LU per query shared by the exact solve and all estimators; cross-check against `fea_engine` sparse contact solver | done (step 6) |
| `circuit_breaker.py`, `report_breaker.py` | E5: calibrated tau_adm routes each query to ROM or FOM fallback; time/accuracy accounting | done (step 6) |
| `run.py` / `__main__.py` | CLI: `python -m contact_rom [e1|e3|all] [--bench b1|b2]` | E1, E3, E4, E5 done; `__main__` pins BLAS to 1 thread |

Rule: functions under 40 lines, no hidden globals, `Config` passed explicitly.

## 5. Reproducibility and verification

- One `rng = np.random.default_rng(cfg.seed)` passed explicitly; NMF seeded through the suite's `seed` argument.
- Log python, numpy, scipy, matplotlib, fea_engine and rom_engine versions; write `config.json` next to outputs.
- Outputs in `results/<run_id>/` as `*.npz`, `*.csv`, `*.png`.
- Assertions: K symmetry and SPD; λ_full ≥ 0 and non-negative bases ≥ 0; array shapes; FOM KKT ≤ 1e-12; **Theorem 1 holds numerically (bound ≥ true error on every query)**.
- Cross-check against the suite: reduced solution from my wrapper equals `ContactROM` output on the same bases (differences at machine precision).
- Functions under 40 lines, math names (`K_mat`, `eta_sign`, `r_rep`), NumPy docstrings, type hints.

## 6. Build order (each step runnable and checked first)

1. Install suite; config; B1; FOM wrapper; check KKT and complementarity.
2. Snapshots, SVD basis, reduced solver via suite; reproduce silent failure on B1 (negative lifted λ with KKT ≈ 1e-16).
3. Estimator components; check bound ≥ error → E1, E2.
4. NMF, CPG, gIS, RID bases; B2; E3 using `compare_contact_models` plus new Pareto axes.
5. Cost model, FLOP and wall-clock timers → E4.
6. B3 and circuit-breaker fallback → E5.
7. Reporting, summary table, acceptance checks. E6 optional.

## 7. Acceptance criteria (Section 5.3)

1. Zero false negatives: every query with ‖R_red‖∞ ≤ 1e-12 and true relative energy error > 5% has η_adm > τ_adm. *(Procedure agreed after step 4: τ_adm is calibrated on a held-out validation set as 0.9 × the smallest indicator value among validation silent failures, then scored on test queries; both Eq. (15) and the square-root-consistent `eta_adm_v2` are reported. Misses and false-alarm rate are reported together.)*
2. θ_eff ∈ [1, 10] on ≥ 90% of queries whose relative energy error exceeds the 5% tolerance (the errors the estimator must catch). The all-queries share is reported as a diagnostic only: Theorem 1's dominant term 2‖λ‖₂η_gap is first-order in the error, so θ_eff grows like err^(-1/2) as the error shrinks (fitted slope −0.54 on B1/NMF). *(Scoping agreed with the user after step 3; the manuscript's Section 5.3 wording has not been changed.)*
3. Spearman(η_bound, true error) > 0.90 across bases and dimensions.
4. Predicted break-even within a factor of 2 of measured, on B1–B3.

## 8. Risks

- Paper figures (`r_rep = 14`, λ ≈ −1.4·10², KKT ≈ 1e-16) may not reproduce; report what is obtained.
- Effectivity and Spearman targets are goals; the summary reports pass/fail honestly.
- The suite's contact ROM is "frictionless, constant gap map, small displacements"; B3 with parameter-varying K(μ) and c(μ) may need the adapters to handle a μ-dependent gap rather than `ContactROM` directly.
- The suite's `rel_tol` rule in greedy sampling is the author's reading of Algorithm 2 (noted as a probable typo in the original); I will not rely on it for CPG.
- B3 at N = 10⁵ is memory-heavy for dense compliance; the sparse route avoids it for the FOM, and the default size is capped.
- API signatures were read from the module listings, not run; step 1 verifies them against a small example before anything is built on them.

## 9. Open decisions (defaults if no objection)

1. Single file vs. small package → **small package** (changed from single file at user request).
2. B3 source → **suite solids + sparse contact solver**.
3. Spec source → **v2 sections**.
