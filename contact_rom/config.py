"""Central configuration: every hyperparameter, tolerance and path lives here."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class B1Config:
    """Benchmark B1: 1D spring chain above a rigid parametrized obstacle."""

    n_dof: int = 60                     # N, chain nodes
    n_contact: int = 30                 # n_c, every (N/n_c)-th node can touch
    k_spring: float = 1.0               # spring stiffness
    k_ground: float = 0.05              # weak grounding spring (keeps K SPD)
    clearance0: float = 1.0             # base obstacle clearance
    load_range: tuple = (0.05, 0.3)      # mu_1: downward load amplitude
    curvature_range: tuple = (1.0, 6.0)  # mu_2: obstacle curvature


@dataclass(frozen=True)
class B2Config:
    """Benchmark B2: clamped membrane (n_grid x n_grid interior nodes) over a curved obstacle."""

    n_grid: int = 18                    # interior grid -> N = n_grid**2 = 324
    n_sub: int = 12                     # contact sub-grid -> n_c = n_sub**2 = 144
    clearance0: float = 0.5             # base obstacle clearance
    curvature: float = 4.0              # obstacle curvature kappa
    load_range: tuple = (0.05, 0.25)    # mu_1: peak load F0
    center_range: tuple = (0.35, 0.65)  # mu_2, mu_3: indenter centre (xc, yc)


@dataclass(frozen=True)
class B3Config:
    """Benchmark B3: quarter-symmetric two-body 3D contact (hex8 blocks), parametrized K(mu)."""

    n_xy: int = 20                      # elements per side in x and y (n_c = (n_xy + 1)^2)
    nz_up: int = 6                      # element layers, upper block
    nz_low: int = 6                     # element layers, lower block
    nu: float = 0.3                     # Poisson ratio (both blocks)
    k_ground: float = 0.01              # weak ground spring on the upper top face (keeps K SPD)
    kappa: float = 0.5                  # curvature of the gap profile delta0 + kappa r^2
    load_sigma: float = 0.5             # width of the Gaussian pressure profile (fraction of the side)
    ratio_range: tuple = (0.5, 2.0)     # mu_1 = E_up / E_low
    load_range: tuple = (0.01, 0.06)    # mu_2 = peak pressure P
    delta_range: tuple = (0.05, 0.3)    # mu_3 = initial clearance delta0
    extrap_load: tuple = (1.3, 2.0)     # extrapolation: P in this multiple of the training maximum


@dataclass(frozen=True)
class B4Config:
    """Benchmark B4: 3D elastic block, clamped on its lateral faces, resting on a rigid curved foundation (constant K)."""

    n_xy: int = 16                      # elements per side in x and y; n_c = (n_xy - 1)^2 candidate bottom nodes
    nz: int = 4                         # element layers
    nu: float = 0.3                     # Poisson ratio (E = 1)
    clearance0: float = 0.05            # base clearance between the block bottom and the foundation
    curvature: float = 1.0              # foundation curvature kappa (clearance grows with the squared distance from the centre)
    load_sigma: float = 0.5             # width of the Gaussian pressure on the top face (fraction of the side)
    load_range: tuple = (0.02, 0.10)    # mu_1: peak pressure
    center_range: tuple = (0.35, 0.65)  # mu_2, mu_3: foundation centre (xc, yc) as fractions of the side


@dataclass(frozen=True)
class Config:
    """Single source of truth for every hyperparameter and path."""

    seed: int = 42
    n_train: int = 60
    n_test: int = 40
    tol_kkt_fom: float = 1e-10          # FOM KKT verification tolerance
    tol_sym: float = 1e-12              # matrix-symmetry tolerance
    lcp_method: str = "active_set"
    r_primal: int = 20                  # POD primal dimension r_u
    r_dual_list: tuple = (4, 8, 12, 16, 20)   # dual dimensions r_lambda
    tol_kkt_red: float = 1e-12          # 'machine-precision' reduced KKT threshold
    err_silent: float = 0.05            # relative energy error counted as failure
    n_val: int = 40                     # held-out validation queries for tau_adm calibration
    tau_margin: float = 0.9             # tau = margin * min indicator over validation silent failures
    c_fact: float = 1.0                 # Eq. (16) factorisation constant
    c_solve: float = 1.0                # Eq. (16) triangular-solve constant
    c_proj: float = 2.0                 # Eq. (17) projection constant (2 FLOPs per multiply-add)
    e4_repeats: int = 5                 # best-of repeats for wall-clock timings
    e4_grids: tuple = (5, 6, 8, 12, 16, 20, 24, 32, 40, 48)   # membrane n_grid values for the N sweep
    e4_n_train: int = 40                # snapshots per sweep size
    e4_n_test: int = 10                 # test queries per sweep size
    e4_grids_b3: tuple = (2, 3, 4, 6, 8, 12, 16, 20)   # B3 n_xy values for the N sweep
    e4_codes_b3_small: tuple = (11, 12, 21)             # B3 very small meshes, code = 10 n_xy + layers (N = 24, 40, 63)
    e4_kind: str = "nmf"                # dual basis used in the sweep
    e4_r_dual: int = 12                 # r_lambda used in the sweep
    e5_kinds: tuple = ("nmf", "cpg")      # dual bases run through the circuit breaker
    e5_r_dual_list: tuple = (8, 16)     # their dual dimensions
    e5_n_train: int = 40                # E5 training snapshots
    e5_n_val: int = 25                  # E5 validation queries (threshold calibration)
    e5_n_test: int = 25                 # E5 in-distribution test queries
    e5_n_extrap: int = 25               # E5 extrapolation queries
    e5_timing_queries: int = 6          # queries used to time the full-order solve
    e5_repeats: int = 2                 # best-of repeats for FOM timing (B3 solves are slow)
    cond_max: float = 1e6               # rank-repair target for cond(B_r K_r^-1 B_r^T)
    r_dual_list_e3: tuple = (4, 8, 12, 16, 24)  # E3 dual dimensions
    out_dir: Path = Path("results")
    run_id: str = "step1"
    b1: B1Config = field(default_factory=B1Config)
    b2: B2Config = field(default_factory=B2Config)
    b3: B3Config = field(default_factory=B3Config)
    b4: B4Config = field(default_factory=B4Config)
