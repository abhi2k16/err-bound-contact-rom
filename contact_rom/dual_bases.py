"""Five dual-basis constructions (manuscript Section 4.1) behind one interface."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy.optimize import nnls

from rom_engine.contact_dual_basis import nonneg_dual_basis, svd_dual_basis
from rom_engine.contact_rid import contact_rows_in_rid, deim_rid_indices, rid_dofs, rid_elements

DUAL_KINDS = ("svd", "nmf", "cpg", "gis", "rid")


@dataclass
class BasisContext:
    """Problem data some constructions need (RID uses geometry; others ignore it)."""

    V_mat: np.ndarray                 # (N, r_u) primal POD basis
    B_mat: sp.spmatrix                # (n_c, N) constraint Jacobian
    elem_dofs: np.ndarray | None      # (n_elem, dofs_per_elem)
    rid_layers: int = 1               # extra element rings around the DEIM points


def gis_dual_basis(Lam_snap: np.ndarray, r_dual: int) -> np.ndarray:
    """Greedy index selection: columns of ``I_{n_c}`` with the largest snapshot energy.

    Columns are orthonormal, so the greedy choice that minimises the lifted
    error ``||lambda_outside||`` is the top-``r_dual`` rows by ``sum_s lambda_s[j]^2``.
    """
    energy = np.sum(Lam_snap**2, axis=1)
    rows = np.sort(np.argsort(energy)[::-1][:r_dual])
    return np.eye(Lam_snap.shape[0])[:, rows]


def cone_residual(U_cols: np.ndarray, target: np.ndarray) -> float:
    """Distance from ``target`` to ``cone(U_cols)`` (non-negative least squares)."""
    if U_cols.shape[1] == 0:
        return float(np.linalg.norm(target))
    coeffs, resid = nnls(U_cols, target)
    return float(resid)


def cpg_dual_basis(Lam_snap: np.ndarray, r_dual: int) -> np.ndarray:
    """Cone-projected greedy: add the snapshot worst approximated by ``cone(U)`` (Benaceur et al.).

    Start from the largest-norm snapshot; columns are normalised snapshots, so
    the basis is non-negative by construction.
    """
    n_c, n_s = Lam_snap.shape
    chosen: list[int] = [int(np.argmax(np.linalg.norm(Lam_snap, axis=0)))]
    for _ in range(min(r_dual, n_s) - 1):
        U_cols = Lam_snap[:, chosen]
        resid = np.array([cone_residual(U_cols, Lam_snap[:, s]) if s not in chosen else -1.0 for s in range(n_s)])
        if resid.max() <= 1e-12:
            break
        chosen.append(int(np.argmax(resid)))
    U_mat = Lam_snap[:, chosen]
    return U_mat / np.maximum(np.linalg.norm(U_mat, axis=0), 1e-300)


def rid_dual_basis(Lam_snap: np.ndarray, r_dual: int, ctx: BasisContext) -> np.ndarray:
    """RID restriction: identity columns on the contact rows inside the reduced integration domain.

    The RID is built from DEIM points of the primal POD basis and of the
    leading multiplier modes (Fauque et al.), grown by ``ctx.rid_layers`` rings.
    """
    if ctx.elem_dofs is None:
        raise ValueError("RID needs element connectivity (ContactProblem.elem_dofs)")
    lam_modes = svd_dual_basis(Lam_snap, r_dual)
    dofs_primal = deim_rid_indices(ctx.V_mat[:, :r_dual])
    rows_dual = deim_rid_indices(lam_modes)
    dofs_dual = sp.csr_matrix(ctx.B_mat)[rows_dual].indices
    dof_idx = np.unique(np.concatenate([dofs_primal, dofs_dual]))
    rid_elems = rid_elements(ctx.elem_dofs, ctx.V_mat.shape[0], dof_idx, layers=ctx.rid_layers)
    P_dofs, _ = rid_dofs(ctx.elem_dofs, ctx.V_mat.shape[0], rid_elems)
    rows = contact_rows_in_rid(ctx.B_mat, P_dofs)
    if rows.size == 0:
        raise ValueError("RID contains no contact rows; increase rid_layers or r_dual")
    return np.eye(Lam_snap.shape[0])[:, rows]


def make_dual_basis(kind: str, Lam_snap: np.ndarray, r_dual: int, ctx: BasisContext | None = None) -> np.ndarray:
    """Dual basis ``U_mat`` of shape ``(n_c, r)`` for ``kind`` in :data:`DUAL_KINDS`."""
    if kind == "svd":
        return svd_dual_basis(Lam_snap, r_dual)
    if kind == "nmf":
        U_mat = nonneg_dual_basis(Lam_snap, r_dual).U
    elif kind == "cpg":
        U_mat = cpg_dual_basis(Lam_snap, r_dual)
    elif kind == "gis":
        U_mat = gis_dual_basis(Lam_snap, r_dual)
    elif kind == "rid":
        U_mat = rid_dual_basis(Lam_snap, r_dual, ctx)
    else:
        raise ValueError(f"unknown dual basis kind: {kind}")
    assert U_mat.min() >= 0.0, f"{kind} basis must be non-negative"
    return U_mat
