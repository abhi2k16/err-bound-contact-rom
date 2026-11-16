"""Geometry and field figures built with the suite's visualization tools (PyVista via ``fea_engine.visualization``).

Run ``python -m scripts.figures.make_geometry_figures``. Outputs go to ``figures/``:
  fig_geometry.pdf    B1 chain, B2 membrane, B3 two-body mesh with loads, supports and active contact
  fig_b3_rom.pdf      B3 full-order vs reduced displacement and error, NMF and SVD dual bases (``compare_fields``)
  fig_multiplier.pdf  B2 contact tractions: full order vs lifted reduced multipliers (adhesive values visible)
Intermediate renderings are written to ``figures/raster/``.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib.colors import TwoSlopeNorm

from contact_rom.audit import build_reduced_model, rom_for_query, solve_reduced_timed
from contact_rom.bases import build_snapshots, pod_primal_basis
from contact_rom.benchmark_3d import GridLayout, build_b3, constrained_dofs
from contact_rom.benchmarks import build_b1, build_b2
from contact_rom.config import Config
from contact_rom.query_context import build_query_context
from fea_engine.mesh import Mesh, box_mesh, rectangle_mesh
from fea_engine.visualization import to_pyvista

pv.OFF_SCREEN = True
OUT, RASTER = Path("figures"), Path("figures/raster")
BLUE, ORANGE, INK, MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "font.size": 8,
                     "mathtext.fontset": "stix", "axes.linewidth": 0.6, "pdf.fonttype": 42})
CFG = Config()


# ----------------------------------------------------------------------------------------------------- B3 helpers
def b3_small() -> tuple:
    """B3 on a small mesh (illustration): ``(problem, layout, b3 config)``."""
    b3 = replace(CFG.b3, n_xy=8, nz_up=3, nz_low=3)
    return build_b3(b3), GridLayout(b3.n_xy, b3.nz_up, b3.nz_low), b3


def two_body_mesh(b3) -> tuple:
    """Upper and lower hex8 blocks (unit elements, interface at z = 0) and their union as one ``Mesh``."""
    upper = box_mesh(b3.n_xy, b3.n_xy, b3.nz_up, b3.n_xy, b3.n_xy, b3.nz_up)
    lower = box_mesh(b3.n_xy, b3.n_xy, b3.nz_low, b3.n_xy, b3.n_xy, b3.nz_low)
    lower.nodes[:, 2] -= b3.nz_low
    nodes = np.vstack([upper.nodes, lower.nodes])
    elems = np.vstack([upper.elements, lower.elements + len(upper.nodes)])
    return upper, lower, Mesh(nodes, elems, dim=3)


def embed_b3(u_free: np.ndarray, layout: GridLayout, b3, n_up: int, mesh: Mesh) -> np.ndarray:
    """Scatter the free-dof vector into a full ``3 n_nodes`` displacement of the union mesh."""
    free = np.setdiff1d(np.arange(3 * layout.n_nodes), constrained_dofs(layout))
    pos = -np.ones(3 * layout.n_nodes, dtype=int)
    pos[free] = np.arange(free.size)
    full = np.zeros(3 * len(mesh.nodes))
    for n, (x, y, z) in enumerate(mesh.nodes):
        body, k = (0, int(round(z))) if n < n_up else (1, int(round(z)) + b3.nz_low)
        gid = layout.node(body, np.array(int(round(x))), np.array(int(round(y))), np.array(k))
        for d in range(3):
            p = pos[3 * int(gid) + d]
            full[3 * n + d] = u_free[p] if p >= 0 else 0.0
    return full


def pick_query(problem, cfg, snaps, V_pod, kinds=("nmf", "svd"), r_dual=8) -> tuple:
    """A test query on which SVD fails silently (error > 5%) and NMF stays accurate; returns lifted solutions."""
    rng = np.random.default_rng(7)
    roms = {k: build_reduced_model(problem, cfg, snaps, V_pod, k, r_dual, repair=True)[0] for k in kinds}
    for mu in problem.sample_mu(rng, 30):
        ctx = build_query_context(problem, mu)
        sols = {k: solve_reduced_timed(rom_for_query(problem, cfg, roms[k], ctx), ctx.f_vec, -ctx.c_vec) for k in kinds}
        err = {k: np.linalg.norm(sols[k].u - ctx.fom.u_vec) / np.linalg.norm(ctx.fom.u_vec) for k in kinds}
        if err["svd"] > 0.2 and err["nmf"] < 0.05:
            return mu, ctx, sols, err
    raise SystemExit("no suitable query found")


# ----------------------------------------------------------------------------------------------------- geometry
def crop_white(path: Path, pad: int = 12) -> None:
    """Trim the white margin of a rendered PNG in place."""
    img = mpimg.imread(path)
    mask = (img[..., :3] < 0.98).any(axis=2)
    rows, cols = np.where(mask.any(axis=1))[0], np.where(mask.any(axis=0))[0]
    r0, r1, c0, c1 = max(rows[0] - pad, 0), rows[-1] + pad, max(cols[0] - pad, 0), cols[-1] + pad
    mpimg.imsave(path, img[r0:r1, c0:c1])


def render_b3_geometry(path: Path) -> None:
    """Two-body model: upper block, lower block, load arrows, clamped base, symmetry planes."""
    problem, layout, b3 = b3_small()
    upper, lower, _ = two_body_mesh(b3)
    pl = pv.Plotter(off_screen=True, window_size=[1100, 900])
    pl.set_background("white")
    pl.add_mesh(to_pyvista(upper), color="#9ec3ee", show_edges=True, edge_color="#2d4f7c", line_width=0.6)
    pl.add_mesh(to_pyvista(lower), color="#e3ded3", show_edges=True, edge_color="#6b665c", line_width=0.6)
    for x in np.linspace(1, b3.n_xy - 1, 4):
        for y in np.linspace(1, b3.n_xy - 1, 4):
            h = 1.7 * np.exp(-((x**2 + y**2) / b3.n_xy**2) / (2 * b3.load_sigma**2))
            pl.add_mesh(pv.Arrow(start=(x, y, b3.nz_up + 0.2 + h), direction=(0, 0, -1), scale=h, tip_length=0.4,
                                 tip_radius=0.12, shaft_radius=0.04), color=ORANGE)
    base = pv.Plane(center=(b3.n_xy / 2, b3.n_xy / 2, -b3.nz_low - 0.05), direction=(0, 0, 1), i_size=b3.n_xy * 1.12,
                    j_size=b3.n_xy * 1.12)
    pl.add_mesh(base, color="#444444", opacity=0.9)
    pl.camera_position = [(2.9 * b3.n_xy, -1.5 * b3.n_xy, 1.9 * b3.n_xy), (b3.n_xy / 2, b3.n_xy / 2, 0.3), (0, 0, 1)]
    pl.screenshot(str(path))
    pl.close()
    crop_white(path)


def render_b2_geometry(path: Path) -> dict:
    """Membrane deflection (FOM) under the curved obstacle, active contact nodes marked."""
    prob = build_b2(CFG.b2)
    n = CFG.b2.n_grid
    mu = np.array([0.18, 0.5, 0.5])
    ctx = build_query_context(prob, mu)
    mesh = rectangle_mesh(1.0, 1.0, n + 1, n + 1)
    z = np.zeros(len(mesh.nodes))
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            z[i * (n + 2) + j] = ctx.fom.u_vec[(i - 1) * n + (j - 1)]
    grid = to_pyvista(mesh)
    pts = grid.points.copy()
    pts[:, 2] = z
    grid.points = pts
    grid.point_data["deflection u"] = z
    xy = np.array([((i + 1) / (n + 1), (j + 1) / (n + 1)) for i in range(n) for j in range(n)])
    start = (n - CFG.b2.n_sub) // 2
    cols = [(i * n + j) for i in range(start, start + CFG.b2.n_sub) for j in range(start, start + CFG.b2.n_sub)]
    active = ctx.fom.lam_vec > 1e-12
    obstacle = rectangle_mesh(0.62, 0.62, 24, 24, x0=0.19, y0=0.19)
    og = to_pyvista(obstacle)
    op = og.points.copy()
    op[:, 2] = CFG.b2.clearance0 + CFG.b2.curvature * ((op[:, 0] - mu[1]) ** 2 + (op[:, 1] - mu[2]) ** 2)
    og.points = op
    pl = pv.Plotter(off_screen=True, window_size=[1100, 900])
    pl.set_background("white")
    pl.add_mesh(grid, scalars="deflection u", cmap="Blues", show_edges=True, edge_color="#7f98b8", line_width=0.4,
                show_scalar_bar=False)
    pl.add_mesh(og, color="#9a9a96", opacity=0.45, show_edges=False)
    act_pts = np.column_stack([xy[np.array(cols)[active]], ctx.fom.u_vec[np.array(cols)[active]]])
    pl.add_mesh(pv.PolyData(act_pts), color=ORANGE, point_size=11, render_points_as_spheres=True)
    pl.camera_position = [(2.6, -1.9, 1.9), (0.5, 0.5, 0.45), (0, 0, 1)]
    pl.screenshot(str(path))
    crop_white(path)
    pl.close()
    return {"n_active": int(active.sum()), "n_c": int(prob.n_contact)}


def fig_geometry() -> None:
    RASTER.mkdir(parents=True, exist_ok=True)
    render_b2_geometry(RASTER / "b2_geometry.png")
    render_b3_geometry(RASTER / "b3_geometry.png")
    prob = build_b1(CFG.b1)
    mu = np.array([0.2, 3.0])
    ctx = build_query_context(prob, mu)
    n, n_c = prob.n_dof, prob.n_contact
    x = np.linspace(-1, 1, n)
    xc = np.linspace(-1, 1, n_c)
    cols = np.asarray(prob.B_mat(mu).argmax(axis=1)).ravel()
    fig = plt.figure(figsize=(6.3, 2.4))
    ax = fig.add_axes([0.04, 0.12, 0.30, 0.74])
    ax.plot(xc, prob.c_vec(mu), color=MUTED, lw=1.2)
    ax.fill_between(xc, prob.c_vec(mu), prob.c_vec(mu).max() + 0.4, color="#d9d8d2")
    ax.plot(x, ctx.fom.u_vec, "-", color=BLUE, lw=1.0)
    ax.plot(x[::2], ctx.fom.u_vec[::2], "o", color=BLUE, ms=2.2)
    act = ctx.fom.lam_vec > 1e-12
    ax.plot(xc[act], ctx.fom.u_vec[cols][act], "o", color=ORANGE, ms=4, mfc="none", mew=1.0)
    for xi in np.linspace(-0.8, 0.8, 6):
        ax.annotate("", xy=(xi, mu[0] * 0.0 + ctx.fom.u_vec.max() * 0.35 + 0.05), xytext=(xi, 0.0),
                    arrowprops={"arrowstyle": "-|>", "color": ORANGE, "lw": 0.7})
    ax.set(xlim=(-1.05, 1.05), ylim=(-0.05, prob.c_vec(mu).max() + 0.25), xticks=[], yticks=[])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title(f"B1: spring chain ($N={n}$, $n_c={n_c}$)", fontsize=8)
    ax.text(0, prob.c_vec(mu).max() + 0.08, "parabolic obstacle", ha="center", fontsize=6.5, color=MUTED)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    legends = {
        "b2_geometry.png": [Patch(fc="#6ea0d4", label="membrane (deflection)"), Patch(fc="#b8b8b4", label="rigid obstacle"),
                            Line2D([], [], marker="o", ls="", color=ORANGE, label="active contact nodes")],
        "b3_geometry.png": [Patch(fc="#9ec3ee", label="upper body"), Patch(fc="#e3ded3", label="lower body"),
                            Patch(fc="#444444", label="clamped base"), Line2D([], [], color=ORANGE, lw=2, label="pressure load")]}
    for pos, name, img in ((0.36, "B2: membrane", "b2_geometry.png"), (0.68, "B3: two-body hex8 contact", "b3_geometry.png")):
        a = fig.add_axes([pos, 0.2, 0.31, 0.68])
        a.imshow(mpimg.imread(RASTER / img))
        a.axis("off")
        a.set_title(name + (" ($N=324$, $n_c=144$)" if name.startswith("B2") else " (quarter model)"), fontsize=8)
        a.legend(handles=legends[img], loc="upper center", bbox_to_anchor=(0.5, 0.0), fontsize=6.5, frameon=False,
                 ncol=2, handlelength=1.2, columnspacing=1.0)
    fig.savefig(OUT / "fig_geometry.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / "fig_geometry.pdf")


# ----------------------------------------------------------------------------------------------------- B3 ROM
def render_b3_panels(union: Mesh, U_full: np.ndarray, U_roms: dict, top: float, path: Path) -> dict:
    """2 x 3 grid of 3-D renderings (full | reduced | error) with the suite's field helpers; labels are added later.

    Fields are drawn on the undeformed mesh: the initial gap is part of the constraint data, not of the mesh, so
    deformed blocks would interpenetrate on screen."""
    from fea_engine.visualization.pyvista_comparison import field_error

    pl = pv.Plotter(shape=(len(U_roms), 3), off_screen=True, window_size=[2400, 750 * len(U_roms)], border=False)
    pl.set_background("white")
    err_top, metrics = {}, {}
    for row, (kind, U_rom) in enumerate(U_roms.items()):
        _, fv, rv, err, met = field_error(union, U_full, U_rom, "displacement", None)
        metrics[kind], err_top[kind] = met, float(err.max())
        for col, (U, vals, clim, cmap) in enumerate(((U_full, fv, (0.0, top), "viridis"), (U_rom, rv, (0.0, top), "viridis"),
                                                     (U_full, err, (0.0, err_top[kind]), "magma"))):
            pl.subplot(row, col)
            grid = to_pyvista(union)
            grid.point_data["f"] = vals
            pl.add_mesh(grid, scalars="f", cmap=cmap, clim=clim, show_scalar_bar=False, show_edges=False)
            pl.camera_position = [(2.2 * union.nodes[:, 0].max(), -1.3 * union.nodes[:, 0].max(), 1.5 * union.nodes[:, 0].max()),
                                  (union.nodes[:, 0].max() / 2, union.nodes[:, 0].max() / 2, 0.0), (0, 0, 1)]
            pl.reset_camera_clipping_range()
    pl.screenshot(str(path))
    pl.close()
    crop_white(path)
    return {"err_top": err_top, "rel_l2": {k: float(m["rel_l2"]) for k, m in metrics.items()}}


def fig_b3_rom() -> dict:
    """Full-order | reduced | error displacement for the NMF and SVD ROM of one B3 query."""
    problem, layout, b3 = b3_small()
    cfg = replace(CFG, b3=b3)
    rng = np.random.default_rng(cfg.seed + 1)
    snaps = build_snapshots(problem, problem.sample_mu(rng, cfg.e5_n_train), cfg.lcp_method)
    V_pod = pod_primal_basis(snaps.U_snap, cfg.r_primal)
    mu, ctx, sols, err = pick_query(problem, cfg, snaps, V_pod)
    upper, _, union = two_body_mesh(b3)
    n_up = len(upper.nodes)
    U_full = embed_b3(ctx.fom.u_vec, layout, b3, n_up, union)
    U_roms = {k: embed_b3(sols[k].u, layout, b3, n_up, union) for k in ("nmf", "svd")}
    mag = lambda U: np.linalg.norm(U.reshape(-1, 3), axis=1)
    top = float(max(mag(U_full).max(), *(mag(U).max() for U in U_roms.values())))
    RASTER.mkdir(parents=True, exist_ok=True)
    info = render_b3_panels(union, U_full, U_roms, top, RASTER / "b3_rom_grid.png")
    img = mpimg.imread(RASTER / "b3_rom_grid.png")
    fig = plt.figure(figsize=(6.3, 4.1))
    gs = fig.add_gridspec(2, 3, left=0.07, right=0.99, top=0.93, bottom=0.12, wspace=0.02, hspace=0.25)
    h, w = img.shape[0] // 2, img.shape[1] // 3
    heads = ("full-order", "reduced", "absolute error")
    for row, kind in enumerate(("nmf", "svd")):
        for col in range(3):
            ax = fig.add_subplot(gs[row, col])
            ax.imshow(img[row * h:(row + 1) * h, col * w:(col + 1) * w])
            ax.axis("off")
            if row == 0:
                ax.set_title(heads[col], fontsize=8)
            if col == 2:
                sm = plt.cm.ScalarMappable(norm=plt.Normalize(0, info["err_top"][kind]), cmap="magma")
                cb = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.05, pad=0.02, aspect=18)
                cb.ax.tick_params(labelsize=6)
                cb.set_ticks([0, info["err_top"][kind]])
                cb.set_ticklabels(["0", f"{info['err_top'][kind]:.2g}"])
        label = {"nmf": "NMF dual basis", "svd": "SVD dual basis"}[kind]
        fig.text(0.015, gs[row, 0].get_position(fig).y0 + 0.5 * gs[row, 0].get_position(fig).height, label, rotation=90,
                 va="center", fontsize=8)
    sm = plt.cm.ScalarMappable(norm=plt.Normalize(0, top), cmap="viridis")
    cax = fig.add_axes([0.20, 0.04, 0.38, 0.018])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=6)
    cb.set_label("displacement magnitude (first two columns, both rows)", fontsize=6.5)
    fig.savefig(OUT / "fig_b3_rom.pdf", dpi=300)
    plt.close(fig)
    stats = {k: {"rel_l2": info["rel_l2"][k], "energy_err": float(err[k]), "lam_min": float(sols[k].lam.min())}
             for k in ("nmf", "svd")}
    print("wrote", OUT / "fig_b3_rom.pdf", stats)
    return stats


# ----------------------------------------------------------------------------------------------------- multipliers
def fig_multiplier() -> dict:
    """B2 contact traction on the 12x12 contact patch: full order vs lifted NMF and SVD multipliers."""
    problem = build_b2(CFG.b2)
    rng = np.random.default_rng(CFG.seed + 1)
    snaps = build_snapshots(problem, problem.sample_mu(rng, CFG.n_train), CFG.lcp_method)
    V_pod = pod_primal_basis(snaps.U_snap, CFG.r_primal)
    roms = {k: build_reduced_model(problem, CFG, snaps, V_pod, k, 8, repair=True)[0] for k in ("nmf", "svd")}
    mu = problem.sample_mu(rng, 1)[0]
    ctx = build_query_context(problem, mu)
    fields = {"full-order": ctx.fom.lam_vec}
    for kind, label in (("nmf", "NMF dual basis"), ("svd", "SVD dual basis")):
        res = solve_reduced_timed(rom_for_query(problem, CFG, roms[kind], ctx), ctx.f_vec, -ctx.c_vec)
        fields[label] = res.lam
    n_sub = CFG.b2.n_sub
    lim = max(abs(v).max() for v in fields.values())
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
    fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.2))
    for ax, (name, lam) in zip(axes, fields.items()):
        im = ax.imshow(lam.reshape(n_sub, n_sub).T, origin="lower", cmap="RdBu_r", norm=norm, extent=(0, 1, 0, 1))
        ax.set_title(name, fontsize=8)
        ax.text(0.03, 0.04, f"min = {lam.min():.2g}", transform=ax.transAxes, fontsize=6.5,
                color="white" if lam.min() < -0.25 * lim else INK)
        ax.set_xticks([])
        ax.set_yticks([])
    cb = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02)
    cb.set_label(r"contact traction $\lambda$ (blue: negative, non-physical tensile)", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    fig.savefig(OUT / "fig_multiplier.pdf", bbox_inches="tight")
    plt.close(fig)
    out = {k: float(v.min()) for k, v in fields.items()}
    print("wrote", OUT / "fig_multiplier.pdf", out)
    return out


if __name__ == "__main__":
    fig_geometry()
    fig_b3_rom()
    fig_multiplier()
