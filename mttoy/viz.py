"""
viz.py  —  Visualisation utilities for the MT-on-condensate model.

Public API
──────────
render_state(st, out_png, ...)        Render a single State to PNG (3D matplotlib)
trim_white(im, ...)                   Crop near-white borders from a PIL image
make_gif(png_files, out_gif, ...)     Assemble a GIF from a list of PNG paths
gif_from_states(states, out_gif, ...) Render a sequence of States into a GIF
save_snapshots(snap_dict, ...)        Save early/late PNGs and return a path LUT
draw_schematic(out_path, ...)         Reaction-scheme diagram with rate equations
"""

from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageChops, ImageOps

from mttoy import State


# ════════════════════════════════════════════════════════════════
# Colour palette
# ════════════════════════════════════════════════════════════════
_C_BOUND = "#2a7fff"    # surface-bound
_C_ESCAPED = "#ff6b1a"  # escaped / free


# ════════════════════════════════════════════════════════════════
# Image helpers
# ════════════════════════════════════════════════════════════════

def trim_white(
    im: Image.Image,
    pad: int = 8,
    bg: Tuple[int, int, int] = (255, 255, 255),
    thr: int = 6,
) -> Image.Image:
    if im.mode != "RGB":
        im = im.convert("RGB")

    bg_im = Image.new("RGB", im.size, bg)
    diff = ImageChops.difference(im, bg_im)
    mask = ImageOps.grayscale(diff).point(lambda p: 255 if p > thr else 0)

    bbox = mask.getbbox()
    if bbox is None:
        return im

    x0, y0, x1, y1 = bbox
    return im.crop((
        max(0, x0 - pad),
        max(0, y0 - pad),
        min(im.size[0], x1 + pad),
        min(im.size[1], y1 + pad),
    ))


# ════════════════════════════════════════════════════════════════
# Core renderer
# ════════════════════════════════════════════════════════════════

def render_state(
    st: State,
    out_png: str,
    lim: Optional[float] = None,
    title: bool = False,
    dpi: int = 180,
    elev: float = 20,
    azim: float = 35,
) -> None:
    """
    Render a single State to PNG.

    Important behavior:
    - Physics is untouched.
    - Every filament segment is drawn individually.
    - Any segment touching a free/escaped node is colored orange.
      This guarantees the peeled-off suffix is visible.
    """
    d = st.droplet
    p = d.params
    R = float(d.R)

    if lim is not None:
        ax_lim = float(lim)
    elif st.mts and any(len(mt.nodes) >= 2 for mt in st.mts):
        all_pts = np.vstack([np.vstack(mt.nodes) for mt in st.mts if len(mt.nodes) >= 2])
        max_extent = float(np.abs(all_pts).max())
        ax_lim = max(1.5 * R, min(1.08 * max_extent, 3.0 * R))
    else:
        ax_lim = 1.5 * R

    fig = plt.figure(figsize=(4.0, 4.0), dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_box_aspect((1, 1, 1))

    # droplet wireframe
    u = np.linspace(0, 2 * np.pi, 56)
    v = np.linspace(0, np.pi, 32)
    xs = R * np.outer(np.cos(u), np.sin(v))
    ys = R * np.outer(np.sin(u), np.sin(v))
    zs = R * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(
        xs, ys, zs,
        rstride=2, cstride=2,
        linewidth=0.35, alpha=0.35, color="steelblue"
    )

    # filaments: draw segment-by-segment
    for mt in st.mts:
        if len(mt.nodes) < 2:
            continue

        pts = np.vstack(mt.nodes)
        nb = list(mt.node_bound)

        # guard against malformed states
        if len(nb) < len(pts):
            nb = nb + [False] * (len(pts) - len(nb))
        elif len(nb) > len(pts):
            nb = nb[:len(pts)]

        for i in range(len(pts) - 1):
            p0 = pts[i]
            p1 = pts[i + 1]
            b0 = bool(nb[i])
            b1 = bool(nb[i + 1])

            # If either endpoint is free, show the segment as escaped.
            is_bound_seg = b0 and b1
            color = _C_BOUND if is_bound_seg else _C_ESCAPED
            lw = 1.2

            ax.plot(
                [p0[0], p1[0]],
                [p0[1], p1[1]],
                [p0[2], p1[2]],
                color=color,
                linewidth=lw,
                alpha=0.9,
                solid_capstyle="round",
            )

        # tip marker + radial spike for peeled (free) tips
        tip = pts[-1]
        tip_bound = bool(nb[-1])
        tip_color = _C_BOUND if tip_bound else _C_ESCAPED

        ax.scatter(
            [tip[0]], [tip[1]], [tip[2]],
            color=tip_color, s=6,
            alpha=0.95, depthshade=False,
        )

    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_zlim(-ax_lim, ax_lim)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.view_init(elev=elev, azim=azim)

    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False
        pane.set_edgecolor("white")
    ax.grid(False)

    if title:
        n_bound = sum(1 for mt in st.mts if mt.bound)
        n_escaped = sum(1 for mt in st.mts if mt.bound and not mt.tip_bound)
        mean_surf = float(np.mean(st.surface.n)) if st.surface is not None else float("nan")
        ax.set_title(
            f"{p.name}\n"
            f"t={st.t/60:.1f} min  MTs={len(st.mts)}  "
            f"bound={n_bound}  esc={n_escaped}\n"
            f"η={p.eta_eff:.2f}  c={p.c_bulk:.1f}  "
            f"n_surf={mean_surf:.1f}",
            fontsize=8,
        )

    fig.savefig(out_png, facecolor="white", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


# ════════════════════════════════════════════════════════════════
# GIF helpers
# ════════════════════════════════════════════════════════════════

def make_gif(
    png_files: Sequence[str],
    out_gif: str,
    duration_ms: int = 140,
) -> None:
    frames = [Image.open(p) for p in png_files]
    if not frames:
        raise ValueError("png_files is empty.")
    frames[0].save(
        out_gif,
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0,
    )


def gif_from_states(
    states: Sequence[State],
    out_gif: str = "./sim.gif",
    every: int = 1,
    duration_ms: int = 120,
    dpi: int = 180,
    elev: float = 20,
    azim: float = 35,
    trim: bool = True,
    trim_pad: int = 8,
) -> str:
    if not states:
        raise ValueError("states is empty.")

    tmpdir = Path(tempfile.mkdtemp(prefix="mtgif_"))
    try:
        seq = sorted(states, key=lambda s: float(getattr(s, "t", 0.0)))
        seq = seq[::max(1, int(every))]
        frames = []

        for i, st in enumerate(seq):
            png = tmpdir / f"frame_{i:04d}.png"
            render_state(st, str(png), title=False, dpi=dpi, elev=elev, azim=azim)

            if trim:
                im = Image.open(png).convert("RGB")
                trim_white(im, pad=trim_pad).save(png)

            frames.append(Image.open(png))

        frames[0].save(
            out_gif,
            save_all=True,
            append_images=frames[1:],
            duration=duration_ms,
            loop=0,
        )
    finally:
        shutil.rmtree(tmpdir)

    return out_gif


# ════════════════════════════════════════════════════════════════
# Snapshot batch saver
# ════════════════════════════════════════════════════════════════

def save_snapshots(
    snap_dict: Dict[str, Sequence[State]],
    render_fn,
    out_dir: str = "./sample_data",
    times_min: Tuple[int, ...] = (10, 50),
    prefix_map: Optional[Dict[str, str]] = None,
    tol_sec: float = 180.0,
    trim: bool = True,
    trim_pad: int = 10,
) -> Dict[Tuple[str, int], str]:
    os.makedirs(out_dir, exist_ok=True)
    if prefix_map is None:
        prefix_map = {k: k for k in snap_dict}

    def _snap_at(snaps: Sequence[State], t_min: int) -> State:
        target = float(t_min) * 60.0
        best = min(snaps, key=lambda s: abs(s.t - target))
        if abs(best.t - target) > tol_sec:
            print(f"[warn] closest snapshot is at {best.t/60:.1f} min (target {t_min} min)")
        return best

    records: list = []
    lut: Dict[Tuple[str, int], str] = {}

    for key, snaps in snap_dict.items():
        tag = prefix_map[key]
        for tmin in times_min:
            tmin = int(tmin)
            state = _snap_at(snaps, tmin)
            fname = f"{tag}_{tmin:02d}min.png"
            fpath = os.path.join(out_dir, fname)

            render_fn(state, fpath)

            if trim:
                im = Image.open(fpath).convert("RGB")
                trim_white(im, pad=trim_pad).save(fpath)

            lut[(key, tmin)] = fpath
            records.append(dict(
                condition=key,
                tag=tag,
                time_min=tmin,
                time_sec=float(state.t),
                filename=fname,
                path=fpath,
            ))

    pd.DataFrame(records).to_csv(os.path.join(out_dir, "snapshot_index.csv"), index=False)
    return lut


# ════════════════════════════════════════════════════════════════
# Reaction-scheme schematic
# ════════════════════════════════════════════════════════════════

def draw_schematic(
    out_path: str = "figures/schematic.pdf",
    dpi: int = 200,
    also_png: bool = True,
) -> str:
    """
    Draw the condensate-MT reaction schematic with rate equations and save to file.

    Parameters
    ----------
    out_path  : output file path (PDF recommended; PNG also saved if also_png=True)
    dpi       : raster resolution for PNG output
    also_png  : if True and out_path ends in .pdf, also save a .png version

    Returns
    -------
    out_path  : the primary saved file path
    """
    import matplotlib.patches as mpatches

    C_DROP = "#dbeafe"       # soft sky-blue — droplet interior
    C_BULK = "#fafaf8"       # near-white warm — bulk region
    C_SURF = "#1f78b4"       # strong blue — surface-bound MT / density
    C_FREE = "#e6550d"       # burnt orange — free (escaped) MT
    C_NUC  = "#31a354"       # forest green — nucleation
    C_ADS  = "#984ea3"       # deep purple — adsorption
    C_DIFF = "#08888a"       # teal — surface diffusion
    C_PEEL = "#b2182b"       # deep red — tip peeling

    FS  = 16
    FSL = 18
    ARROW_KW = dict(arrowstyle="-|>", mutation_scale=18, lw=2.2)

    plt.rcParams["font.family"] = "Arial"
    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.55, 1.65)
    ax.set_ylim(-1.45, 1.48)

    R = 0.50

    def pt(deg, r=R):
        th = np.radians(deg)
        return np.array([r * np.cos(th), r * np.sin(th)])

    # backgrounds
    ax.add_patch(plt.Circle((0, 0), 1.5, color=C_BULK, zorder=0))
    ax.add_patch(plt.Circle((0, 0), R,   color=C_DROP, zorder=1, lw=2.2, ec="#3a6ea8"))

    # region labels
    ax.text(0, -0.05, "Condensate\ndroplet", ha="center", va="center",
            fontsize=FSL, color="#3a5a80", style="italic", zorder=5)
    ax.text(0, 1.38, r"Dilute phase ($c_\mathrm{dilute},\;D_\mathrm{3D}$)",
            ha="center", va="center", fontsize=FSL, color="#7a6020", zorder=5)

    # 1. surface tubulin dots
    for th in np.concatenate([np.linspace(95, 172, 18), np.linspace(188, 265, 18)]):
        ax.plot(*pt(th, R), 'o', ms=4.5, color=C_SURF, alpha=0.55, zorder=3)
    ps = pt(225, R + 0.12)
    ax.text(ps[0] - 0.08, ps[1] + 0.10, r"$n(\mathbf{u},t)$",
            ha="center", va="top", fontsize=FS, color=C_SURF,
            fontweight="bold", zorder=6)

    # 2. adsorption arrow (110°)
    p_out = pt(110, 1.38)
    p_in  = pt(110, R + 0.02)
    ax.annotate("", xy=p_in, xytext=p_out,
                arrowprops=dict(color=C_ADS, **ARROW_KW), zorder=6)
    mid_ads = (p_out + p_in) / 2
    ax.text(mid_ads[0] + 0.09, mid_ads[1],
            r"$k_\mathrm{ads} = \dfrac{k_0^\mathrm{ads}\,D_\mathrm{3D}}{R}$",
            ha="left", va="center", fontsize=FS, color=C_ADS, zorder=7)

    # 3. surface diffusion arc (150–188°)
    r_arc = R + 0.10
    theta_arc = np.linspace(np.radians(150), np.radians(188), 60)
    xarc = r_arc * np.cos(theta_arc)
    yarc = r_arc * np.sin(theta_arc)
    ax.plot(xarc[3:-3], yarc[3:-3], color=C_DIFF, lw=2.4, zorder=6)
    for end_idx, nbr_idx in [(3, 5), (-4, -6)]:
        ep   = np.array([xarc[end_idx], yarc[end_idx]])
        nbr  = np.array([xarc[nbr_idx], yarc[nbr_idx]])
        tang = ep - nbr
        tang /= np.linalg.norm(tang) + 1e-9
        ax.annotate("", xy=ep + 0.045 * tang, xytext=ep - 0.01 * tang,
                    arrowprops=dict(color=C_DIFF, arrowstyle="-|>",
                                    mutation_scale=15, lw=1.8), zorder=9)
    mid_diff = np.array([xarc[30], yarc[30]])
    ax.text(mid_diff[0] - 0.13, mid_diff[1],
            r"$D_\mathrm{surf} = 1/\eta_\mathrm{eff}$",
            ha="right", va="center", fontsize=FS, color=C_DIFF, zorder=7)

    # 4. nucleation star (268°)
    p_nuc = pt(268, R)
    ax.plot(*p_nuc, '*', ms=17, color=C_NUC, zorder=7)
    ax.text(p_nuc[0] + 0.02, p_nuc[1] - 0.10,
            r"$k_\mathrm{nuc} \propto D_\mathrm{surf}"
            r"\langle (n - n_\mathrm{thr})^{n_H} \rangle_+$",
            ha="center", va="top", fontsize=FS, color=C_NUC, zorder=7)

    # 5. bound MT (300–342°)
    bound_angles = np.linspace(300, 342, 8)
    bound_pts = np.array([pt(a) for a in bound_angles])
    ax.plot(bound_pts[:, 0], bound_pts[:, 1],
            color=C_SURF, lw=5.0, solid_capstyle="round", zorder=5)
    ax.plot(*bound_pts[0], '*', ms=14, color=C_NUC, zorder=6)

    tip_b  = bound_pts[-1]
    tang_b = pt(347) - pt(337)
    tang_b /= np.linalg.norm(tang_b)
    ax.annotate("", xy=tip_b + 0.19 * tang_b, xytext=tip_b + 0.02 * tang_b,
                arrowprops=dict(color=C_SURF, **ARROW_KW), zorder=8)
    ax.text(tip_b[0] + 0.04, tip_b[1] - 0.10,
            r"$k_\mathrm{grow}^\mathrm{bound} = k_+^\mathrm{surf}"
            r"\cdot\dfrac{D_\mathrm{surf}}{D_\mathrm{surf}+D_0}$",
            ha="left", va="top", fontsize=FS, color=C_SURF, zorder=7)

    # 6. tip peeling (dashed, to upper-right)
    peel_end = pt(18, 1.30)
    ax.annotate("", xy=peel_end, xytext=tip_b,
                arrowprops=dict(color=C_PEEL, linestyle="dashed",
                                arrowstyle="-|>", mutation_scale=15, lw=2.0),
                zorder=7)
    ax.text(peel_end[0] + 0.06, peel_end[1],
            r"$k_\mathrm{off} = e^{-\alpha\eta_\mathrm{eff}}\,\sigma(L)$",
            ha="left", va="center", fontsize=FS, color=C_PEEL, zorder=7)

    # 7. free MT (38°, growing toward 50°)
    free_anchor = pt(38, R)
    free_dir    = np.array([np.cos(np.radians(50)), np.sin(np.radians(50))])
    free_pts    = np.array([free_anchor + i * 0.14 * free_dir for i in range(7)])
    ax.plot(free_pts[:, 0], free_pts[:, 1],
            color=C_FREE, lw=5.0, solid_capstyle="round", zorder=5)
    ax.plot(*free_pts[0], '*', ms=14, color=C_NUC, zorder=6)

    tip_f = free_pts[-1]
    ax.annotate("", xy=tip_f + 0.19 * free_dir, xytext=tip_f + 0.02 * free_dir,
                arrowprops=dict(color=C_FREE, **ARROW_KW), zorder=8)
    ax.text(tip_f[0] + 0.12, tip_f[1] + 0.08,
            r"$k_\mathrm{grow}^\mathrm{free} = k_0^\mathrm{on}\,D_\mathrm{3D}\,c_\mathrm{dilute}$",
            ha="left", va="bottom", fontsize=FS, color=C_FREE, zorder=7)

    # legend
    legend_elements = [
        mpatches.Patch(color=C_SURF,  label="Surface-bound MT"),
        mpatches.Patch(color=C_FREE,  label="Free growing MT"),
        mpatches.Patch(color=C_NUC,   label="Nucleation site"),
        mpatches.Patch(color=C_ADS,   label="Adsorption"),
        mpatches.Patch(color=C_DIFF,  label="Surface diffusion"),
        mpatches.Patch(color=C_PEEL,  label="Tip detachment"),
    ]
    ax.legend(handles=legend_elements, loc="lower right",
              fontsize=FS, framealpha=0.93, edgecolor="#aaaaaa",
              bbox_to_anchor=(1.18, -0.04))

    fig.savefig(out_path, bbox_inches="tight", dpi=dpi)

    if also_png and out_path.endswith(".pdf"):
        png_path = out_path[:-4] + ".png"
        fig.savefig(png_path, bbox_inches="tight", dpi=dpi)

    plt.close(fig)
    return out_path
