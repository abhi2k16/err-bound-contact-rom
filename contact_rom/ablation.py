"""Ablation of the framework components on B1 and B2 (constant stiffness).

Variants, all with primal dimension ``r_primal`` and (where a dual basis is used) ``r_lambda = 16``:

``primal``       primal-only ROM: displacements reduced, all ``n_c`` multipliers kept (identity dual basis),
``pd_norepair``  primal-dual ROM with a dual basis (SVD or NMF) and no primal rank repair,
``pd_repair``    the same with supremizer rank repair,
``screened``     ``pd_repair`` plus the interface screen of Eq. (15) with the hybrid threshold calibrated on
                 validation queries; flagged queries are re-solved at full order.

Per variant we report the lifted error, the silent-failure share, the Schur condition number, the online FLOPs and
measured reduced-solve time, and for ``screened`` the fallback share, the errors of the accepted queries and the net
speedup against the full-order solve (constant-``K`` regime: factorisation and compliance are offline).
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from rom_engine.contact_dual_basis import ContactROM

from .audit import build_reduced_model, evaluate_query
from .bases import pod_primal_basis
from .config import Config
from .cost import ReducedSizes, flops_rom, problem_sizes
from .estimator import tau_adm_from_snapshots
from .query_context import build_query_context
from .repair import schur_condition
from .timing import time_fom_queries
from .calibration_study import threshold

R_DUAL = 16
KINDS = ("svd", "nmf", "cpg", "gis")
SEEDS = (42, 43, 44, 45, 46)


def stack(rows: list) -> dict:
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


def run_bench(cfg: Config, builder, prepare, bench: str) -> dict:
    """Pooled per-query arrays for every variant over the seeds: ``{variant_key: {metric: array}}``."""
    pooled: dict = {}
    for seed in SEEDS:
        run_cfg = replace(cfg, seed=seed)
        problem = builder(run_cfg)
        snaps, mu_test, mu_val = prepare(run_cfg, problem)
        V_pod = pod_primal_basis(snaps.U_snap, run_cfg.r_primal)
        C_lambda = float(np.linalg.norm(snaps.Lam_snap, axis=0).max())
        K0, B0 = problem.K_mat(snaps.mu_set[0]), problem.B_mat(snaps.mu_set[0])
        sizes = problem_sizes(K0, B0)
        tau_rel = tau_adm_from_snapshots(snaps.U_snap, K0, run_cfg.err_silent)
        t_fom = time_fom_queries(problem, mu_test, run_cfg.e4_repeats).t_affine
        models = {("primal", "-"): (ContactROM(K0, -B0, V_pod, np.eye(sizes.n_contact), method=run_cfg.lcp_method),
                                    {"r_dual_eff": float(sizes.n_contact), "r_rep": 0.0, "t_eta": 0.0}, V_pod,
                                    np.eye(sizes.n_contact))}
        for kind in KINDS:
            for repair, name in ((False, "pd_norepair"), (True, "pd_repair")):
                rom, meta = build_reduced_model(problem, run_cfg, snaps, V_pod, kind, R_DUAL, repair)
                models[(name, kind)] = (rom, meta, rom.U, rom.U_lam)
        ctxs_test = [build_query_context(problem, mu) for mu in mu_test]
        ctxs_val = [build_query_context(problem, mu) for mu in mu_val]
        for key, (rom, meta, V_mat, U_mat) in models.items():
            test = stack([evaluate_query(c, rom, C_lambda) for c in ctxs_test])
            cond = schur_condition(V_mat, U_mat, K0, B0)[0]
            r_u, r_l = int(run_cfg.r_primal + meta["r_rep"]), int(meta["r_dual_eff"])
            flops = flops_rom(sizes, ReducedSizes(r_u, r_l, float(np.median(test["k_rom"])),
                                                  float(np.median(test["m_r_act"]))), run_cfg.c_proj, False)
            slot = pooled.setdefault("|".join(key), {"cond": [], "flops": [], "r_u": [], "r_l": [], "arr": []})
            slot["cond"].append(float(cond)); slot["flops"].append(flops); slot["r_u"].append(r_u); slot["r_l"].append(r_l)
            slot["arr"].append({"err": test["err_energy"], "res": test["res_red"], "t_rom": test["t_online"],
                                "t_fom": t_fom, "ok": test["ok"]})
            if key[0] == "pd_repair":
                val = stack([evaluate_query(c, rom, C_lambda) for c in ctxs_val])
                tau = threshold("hybrid", val, "eta_adm", run_cfg, tau_rel)
                flagged = test["eta_adm"] > tau
                t_screen = test["t_online"] + meta["t_eta"] + flagged * t_fom
                s_slot = pooled.setdefault("|".join(("screened", key[1])), {"cond": [], "flops": [], "r_u": [], "r_l": [], "arr": []})
                s_slot["cond"].append(float(cond)); s_slot["flops"].append(flops); s_slot["r_u"].append(r_u); s_slot["r_l"].append(r_l)
                s_slot["arr"].append({"err": np.where(flagged, 0.0, test["err_energy"]), "res": test["res_red"],
                                      "t_rom": t_screen, "t_fom": t_fom, "ok": test["ok"], "flagged": flagged.astype(float),
                                      "err_raw": test["err_energy"]})
    return pooled


def summarize(cfg: Config, pooled: dict) -> dict:
    out = {}
    for key, slot in pooled.items():
        cat = {k: np.concatenate([a[k] for a in slot["arr"]]) for k in slot["arr"][0]}
        err, silent = cat["err"], (cat["res"] <= cfg.tol_kkt_red) & (cat["err"] > cfg.err_silent)
        row = {"err_med": float(np.median(err)), "err_p90": float(np.quantile(err, 0.9)), "err_max": float(err.max()),
               "share_above_5pct": float(np.mean(err > cfg.err_silent)), "silent_share": float(np.mean(silent)),
               "cond_med": float(np.median(slot["cond"])), "cond_max": float(np.max(slot["cond"])),
               "flops": float(np.median(slot["flops"])), "r_u": float(np.median(slot["r_u"])), "r_l": float(np.median(slot["r_l"])),
               "t_rom_mean": float(np.mean(cat["t_rom"])), "t_fom_mean": float(np.mean(cat["t_fom"])),
               "speedup": float(np.mean(cat["t_fom"]) / np.mean(cat["t_rom"])), "n": int(err.size),
               "lcp_failures": int(np.sum(cat["ok"] < 0.5))}
        if "flagged" in cat:
            acc = cat["flagged"] < 0.5
            row.update({"fallback_share": float(1 - acc.mean()),
                        "accepted_err_max": float(cat["err_raw"][acc].max()) if acc.any() else float("nan"),
                        "accepted_above_5pct": int(np.sum(cat["err_raw"][acc] > cfg.err_silent))})
        out[key] = row
    return out


def run(cfg: Config, builders: dict, prepare, bench: str) -> Path:
    summary = summarize(cfg, run_bench(cfg, builders[bench], prepare, bench))
    path = Path(cfg.out_dir) / "ablation" / f"ablation_{bench}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"bench": bench, "seeds": list(SEEDS), "r_dual": R_DUAL, "rows": summary}, indent=1), encoding="utf-8")
    return path
