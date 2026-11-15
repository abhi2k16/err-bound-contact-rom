"""Benchmark B3: quarter-symmetric 3D two-body contact with parameter-varying stiffness.

Two stacked hexahedral elastic blocks (upper over lower) on a structured grid. The
interface is node-to-node and frictionless; the gap is ``delta0 + kappa r^2`` (a curved
indenter profile). Parameters ``mu = (E_up/E_low, P, delta0)``.

Element stiffness comes from ``fea_engine``'s ``Hex8Solid3D`` (one reference element:
every cell of the uniform grid is identical). Symmetry planes ``x = 0`` and ``y = 0`` are
eliminated; the lower block is clamped at its base; a weak ground spring on the upper
top face keeps ``K`` SPD before contact (a stand-in for flange bolts).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from fea_engine.elements.solids import Hex8Solid3D
from fea_engine.material import D_solid3d, Material

from .benchmarks import ContactProblem
from .config import B3Config, B4Config

HEX_CORNERS = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
CORNER_OFFSETS = HEX_CORNERS.astype(int)


@dataclass
class GridLayout:
    """Node numbering for the two stacked blocks (all nodes, before eliminating constrained dofs)."""

    n_xy: int
    nz_up: int
    nz_low: int

    @property
    def nodes_per_body(self) -> tuple:
        plane = (self.n_xy + 1) ** 2
        return plane * (self.nz_up + 1), plane * (self.nz_low + 1)

    def node(self, body: int, i: np.ndarray, j: np.ndarray, k: np.ndarray) -> np.ndarray:
        """Global node id of grid index ``(i, j, k)`` in ``body`` (0 = upper, 1 = lower)."""
        n1 = self.n_xy + 1
        offset = 0 if body == 0 else self.nodes_per_body[0]
        return offset + (k * n1 + j) * n1 + i

    @property
    def n_nodes(self) -> int:
        return sum(self.nodes_per_body)


def reference_stiffness(nu: float) -> np.ndarray:
    """Unit-cube Hex8 stiffness for ``E = 1`` from the suite, shape (24, 24)."""
    D_mat = D_solid3d(Material(E=1.0, nu=nu))
    return Hex8Solid3D().stiffness(HEX_CORNERS, D_mat)


def body_connectivity(layout: GridLayout, body: int, nz: int) -> np.ndarray:
    """Hex node connectivity of one body, shape (n_elem, 8), suite node order."""
    ex, ey, ez = np.meshgrid(np.arange(layout.n_xy), np.arange(layout.n_xy), np.arange(nz), indexing="ij")
    ex, ey, ez = ex.ravel(), ey.ravel(), ez.ravel()
    cols = [layout.node(body, ex + o[0], ey + o[1], ez + o[2]) for o in CORNER_OFFSETS]
    return np.column_stack(cols)


def assemble_body(Ke: np.ndarray, conn: np.ndarray, n_nodes: int) -> sp.csr_matrix:
    """Scatter the identical element matrix over all cells of one body into an (3 n_nodes)^2 matrix."""
    dofs = (3 * conn[:, :, None] + np.arange(3)[None, None, :]).reshape(len(conn), 24)
    rows = np.repeat(dofs, 24, axis=1).ravel()
    cols = np.tile(dofs, (1, 24)).ravel()
    vals = np.tile(Ke.ravel(), len(conn))
    return sp.coo_matrix((vals, (rows, cols)), shape=(3 * n_nodes, 3 * n_nodes)).tocsr()


def constrained_dofs(layout: GridLayout) -> np.ndarray:
    """Dofs eliminated: symmetry planes (u_x on i = 0, u_y on j = 0) and the clamped lower base."""
    n1 = layout.n_xy + 1
    fixed = []
    for body, nz in ((0, layout.nz_up), (1, layout.nz_low)):
        kk, jj, ii = np.meshgrid(np.arange(nz + 1), np.arange(n1), np.arange(n1), indexing="ij")
        ids = layout.node(body, ii.ravel(), jj.ravel(), kk.ravel())
        fixed.append(3 * ids[ii.ravel() == 0])
        fixed.append(3 * ids[jj.ravel() == 0] + 1)
    base = layout.node(1, *[a.ravel() for a in np.meshgrid(np.arange(n1), np.arange(n1), [0], indexing="ij")])
    fixed.append(np.concatenate([3 * base, 3 * base + 1, 3 * base + 2]))
    return np.unique(np.concatenate(fixed))


def top_face_nodes(layout: GridLayout) -> tuple:
    """Upper-body top-face node ids and their tributary weights (quarter model), plus (x, y) coordinates."""
    n1 = layout.n_xy + 1
    ii, jj = np.meshgrid(np.arange(n1), np.arange(n1), indexing="ij")
    ids = layout.node(0, ii.ravel(), jj.ravel(), np.full(ii.size, layout.nz_up))
    edge = lambda a: np.where((a == 0) | (a == layout.n_xy), 0.5, 1.0)
    return ids, edge(ii.ravel()) * edge(jj.ravel()), ii.ravel().astype(float), jj.ravel().astype(float)


def contact_pairs(layout: GridLayout) -> tuple:
    """Upper-bottom and lower-top node ids of the matched interface and their (x, y) positions."""
    n1 = layout.n_xy + 1
    ii, jj = np.meshgrid(np.arange(n1), np.arange(n1), indexing="ij")
    ii, jj = ii.ravel(), jj.ravel()
    up = layout.node(0, ii, jj, np.zeros_like(ii))
    low = layout.node(1, ii, jj, np.full_like(ii, layout.nz_low))
    return up, low, ii.astype(float), jj.astype(float)


def build_b3(cfg: B3Config) -> ContactProblem:
    """Benchmark B3 as a :class:`ContactProblem` (K varies with the stiffness ratio)."""
    layout = GridLayout(cfg.n_xy, cfg.nz_up, cfg.nz_low)
    Ke = reference_stiffness(cfg.nu)
    n_all = 3 * layout.n_nodes
    K_up = assemble_body(Ke, body_connectivity(layout, 0, cfg.nz_up), layout.n_nodes)
    K_low = assemble_body(Ke, body_connectivity(layout, 1, cfg.nz_low), layout.n_nodes)
    top_ids, top_w, top_x, top_y = top_face_nodes(layout)
    ground = sp.coo_matrix((np.full(top_ids.size, cfg.k_ground), (3 * top_ids + 2, 3 * top_ids + 2)),
                           shape=(n_all, n_all)).tocsr()
    free = np.setdiff1d(np.arange(n_all), constrained_dofs(layout))
    pos = -np.ones(n_all, dtype=int)
    pos[free] = np.arange(free.size)
    K_up_f, K_low_f, ground_f = (m[free][:, free].tocsr() for m in (K_up, K_low, ground))
    up_ids, low_ids, cx, cy = contact_pairs(layout)
    rows = np.arange(up_ids.size)
    B_const = sp.csr_matrix((np.concatenate([-np.ones(rows.size), np.ones(rows.size)]),
                             (np.concatenate([rows, rows]), np.concatenate([pos[3 * up_ids + 2], pos[3 * low_ids + 2]]))),
                            shape=(rows.size, free.size))
    radius2 = (cx**2 + cy**2) / (cfg.n_xy**2)                      # normalised distance from the symmetry axis
    bump = np.exp(-((top_x**2 + top_y**2) / cfg.n_xy**2) / (2 * cfg.load_sigma**2))
    area = top_w * bump

    def K_mat(mu: np.ndarray) -> sp.csr_matrix:
        return (float(mu[0]) * K_up_f + K_low_f + ground_f).tocsr()

    def f_vec(mu: np.ndarray) -> np.ndarray:
        vec = np.zeros(free.size)
        vec[pos[3 * top_ids + 2]] = -float(mu[1]) * area
        return vec

    def c_vec(mu: np.ndarray) -> np.ndarray:
        return float(mu[2]) + cfg.kappa * radius2

    def sampler(ranges: tuple):
        def sample(rng: np.random.Generator, n_samples: int) -> np.ndarray:
            lo, hi = np.array([r[0] for r in ranges]), np.array([r[1] for r in ranges])
            return lo + rng.random((n_samples, 3)) * (hi - lo)
        return sample

    # element connectivity in the free-dof numbering (constrained dofs dropped, so rows have different lengths);
    # only the RID dual basis uses it
    elem_dofs = []
    for body, nz in ((0, cfg.nz_up), (1, cfg.nz_low)):
        conn = body_connectivity(layout, body, nz)
        dofs = pos[(3 * conn[:, :, None] + np.arange(3)[None, None, :]).reshape(len(conn), 24)]
        elem_dofs.extend(row[row >= 0] for row in dofs)

    in_ranges = (cfg.ratio_range, cfg.load_range, cfg.delta_range)
    load_hi = cfg.load_range[1]
    out_ranges = (cfg.ratio_range, (cfg.extrap_load[0] * load_hi, cfg.extrap_load[1] * load_hi), cfg.delta_range)
    return ContactProblem("B3", free.size, rows.size, K_mat, f_vec, lambda mu: B_const, c_vec, sampler(in_ranges),
                          elem_dofs=elem_dofs, k_varies=True, sample_mu_extrap=sampler(out_ranges))


def build_b4(cfg: B4Config) -> ContactProblem:
    """Benchmark B4: a single 3D elastic block, clamped on its four lateral faces, resting on a rigid curved foundation.

    Constant stiffness. The load is a Gaussian pressure on the top face; the foundation centre ``(xc, yc)`` moves with the
    parameters, so that the contact patch translates over the candidate nodes of the bottom face (no symmetry plane is
    used). ``mu = (P, xc, yc)``. Constraint (paper convention): ``-u_z <= c`` at each candidate bottom node, with
    ``c = clearance0 + kappa ((x - xc)^2 + (y - yc)^2)`` and coordinates as fractions of the side.
    """
    n, nz = cfg.n_xy, cfg.nz
    layout = GridLayout(n, nz, 1)
    n_nodes = layout.nodes_per_body[0]
    n_all = 3 * n_nodes
    K_all = assemble_body(reference_stiffness(cfg.nu), body_connectivity(layout, 0, nz), n_nodes)
    n1 = n + 1
    kk, jj, ii = np.meshgrid(np.arange(nz + 1), np.arange(n1), np.arange(n1), indexing="ij")
    kk, jj, ii = kk.ravel(), jj.ravel(), ii.ravel()
    ids = layout.node(0, ii, jj, kk)
    lateral = ids[(ii == 0) | (ii == n) | (jj == 0) | (jj == n)]
    fixed = np.unique(np.concatenate([3 * lateral, 3 * lateral + 1, 3 * lateral + 2]))
    free = np.setdiff1d(np.arange(n_all), fixed)
    pos = -np.ones(n_all, dtype=int)
    pos[free] = np.arange(free.size)
    K_const = K_all[free][:, free].tocsr()
    # contact candidates: interior nodes of the bottom face (k = 0)
    bi, bj = np.meshgrid(np.arange(1, n), np.arange(1, n), indexing="ij")
    bi, bj = bi.ravel(), bj.ravel()
    bottom = layout.node(0, bi, bj, np.zeros_like(bi))
    rows = np.arange(bottom.size)
    B_const = sp.csr_matrix((-np.ones(rows.size), (rows, pos[3 * bottom + 2])), shape=(rows.size, free.size))
    cx, cy = bi / n, bj / n
    # top-face Gaussian pressure, centred on the block
    ti, tj = np.meshgrid(np.arange(1, n), np.arange(1, n), indexing="ij")
    ti, tj = ti.ravel(), tj.ravel()
    top = layout.node(0, ti, tj, np.full_like(ti, nz))
    bump = np.exp(-(((ti / n - 0.5) ** 2 + (tj / n - 0.5) ** 2)) / (2 * cfg.load_sigma ** 2 * 0.25))

    def f_vec(mu: np.ndarray) -> np.ndarray:
        vec = np.zeros(free.size)
        vec[pos[3 * top + 2]] = -float(mu[0]) * bump
        return vec

    def c_vec(mu: np.ndarray) -> np.ndarray:
        return cfg.clearance0 + cfg.curvature * ((cx - mu[1]) ** 2 + (cy - mu[2]) ** 2)

    def sample_mu(rng: np.random.Generator, n_samples: int) -> np.ndarray:
        lo = np.array([cfg.load_range[0], cfg.center_range[0], cfg.center_range[0]])
        hi = np.array([cfg.load_range[1], cfg.center_range[1], cfg.center_range[1]])
        return lo + rng.random((n_samples, 3)) * (hi - lo)

    conn = body_connectivity(layout, 0, nz)
    dofs = pos[(3 * conn[:, :, None] + np.arange(3)[None, None, :]).reshape(len(conn), 24)]
    elem_dofs = [row[row >= 0] for row in dofs]
    return ContactProblem("B4", free.size, rows.size, lambda mu: K_const, f_vec, lambda mu: B_const, c_vec, sample_mu,
                          elem_dofs=elem_dofs)
