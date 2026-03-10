"""
pyvista_viz.py  —  PyVista-based 3D visualisation for the MT-on-condensate model.

Public API
──────────
render_state_pv(st, ...)          Render a single State to an off-screen PNG
live_viewer(states, ...)          Animate a list of States in an interactive window
snapshot_pv(st, out_path, ...)    Quick single-frame export
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np

try:
    import pyvista as pv
except ImportError as e:
    raise ImportError("pyvista is required: pip install pyvista") from e

from mttoy.core import State

# ── Colour palette ────────────────────────────────────────────────────────────
_C_BOUND   = (0.16, 0.44, 1.0)   # blue  — surface-bound tip
_C_ESCAPED = (1.0,  0.42, 0.10)  # orange — escaped tip
_CMAP_DENSITY = "plasma"          # surface density colourmap


# ════════════════════════════════════════════════════════════════
# Internal helpers
# ════════════════════════════════════════════════════════════════

def _sphere_mesh(center: np.ndarray, R: float, resolution: int = 80) -> pv.PolyData:
    sphere = pv.Sphere(radius=R, center=center, theta_resolution=resolution,
                       phi_resolution=resolution)
    return sphere


def _density_on_sphere(sphere: pv.PolyData, n: np.ndarray,
                        center: np.ndarray, R: float) -> pv.PolyData:
    """
    Map the flat (Ny, Nx) density grid onto sphere vertices via spherical
    coordinates.  The grid is treated as a lat/lon texture:
        x-axis → longitude (0..2π)
        y-axis → latitude  (-π/2..π/2)
    """
    pts = np.array(sphere.points) - center          # shift to origin
    r   = np.linalg.norm(pts, axis=1, keepdims=True).clip(1e-12)
    pts_n = pts / r                                  # unit vectors

    lon = np.arctan2(pts_n[:, 1], pts_n[:, 0])      # -π .. π
    lat = np.arcsin(pts_n[:, 2].clip(-1, 1))         # -π/2 .. π/2

    Ny, Nx = n.shape
    ix = ((lon + np.pi) / (2 * np.pi) * Nx).astype(int).clip(0, Nx - 1)
    iy = ((lat + np.pi / 2) / np.pi * Ny).astype(int).clip(0, Ny - 1)

    sphere["density"] = n[iy, ix].astype(float)
    return sphere


def _add_mts(plotter: pv.Plotter, st: State, tube_radius: float = 0.08):
    """Draw each MT as a tube; colour by bound/escaped status."""
    for mt in st.mts:
        if len(mt.nodes) < 2:
            continue
        pts = np.vstack(mt.nodes)
        spline = pv.Spline(pts, n_points=max(len(pts) * 4, 20))
        color = _C_BOUND if mt.tip_bound else _C_ESCAPED
        plotter.add_mesh(
            spline.tube(radius=tube_radius),
            color=color,
            smooth_shading=True,
        )


def _base_plotter(off_screen: bool = False, window_size=(900, 800)) -> pv.Plotter:
    pl = pv.Plotter(off_screen=off_screen, window_size=window_size)
    pl.set_background("black")
    pl.enable_anti_aliasing()
    return pl


def _populate_plotter(pl: pv.Plotter, st: State,
                      show_density: bool = True,
                      show_mts: bool = True,
                      sphere_opacity: float = 0.55,
                      sphere_resolution: int = 80,
                      tube_radius: float = 0.08):
    center = st.droplet.center
    R      = st.droplet.R

    sphere = _sphere_mesh(center, R, resolution=sphere_resolution)

    if show_density and st.surface is not None:
        sphere = _density_on_sphere(sphere, st.surface.n, center, R)
        pl.add_mesh(sphere, scalars="density", cmap=_CMAP_DENSITY,
                    opacity=sphere_opacity, smooth_shading=True,
                    show_scalar_bar=True, scalar_bar_args=dict(
                        title="Surface density n", title_font_size=12,
                        label_font_size=10, color="white"))
    else:
        pl.add_mesh(sphere, color="#4488cc", opacity=sphere_opacity,
                    smooth_shading=True)

    if show_mts:
        _add_mts(pl, st, tube_radius=tube_radius)

    # time label
    pl.add_text(f"t = {st.t:.1f} s", position="upper_left",
                font_size=10, color="white")


# ════════════════════════════════════════════════════════════════
# Public API
# ════════════════════════════════════════════════════════════════

def snapshot_pv(
    st: State,
    out_path: str | Path,
    *,
    show_density: bool = True,
    show_mts: bool = True,
    sphere_opacity: float = 0.55,
    sphere_resolution: int = 80,
    tube_radius: float = 0.08,
    window_size: tuple = (900, 800),
    camera_position: Optional[tuple] = None,
) -> Path:
    """Render a single State off-screen and save to PNG."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    pl = _base_plotter(off_screen=True, window_size=window_size)
    _populate_plotter(pl, st, show_density=show_density, show_mts=show_mts,
                      sphere_opacity=sphere_opacity,
                      sphere_resolution=sphere_resolution,
                      tube_radius=tube_radius)

    if camera_position is not None:
        pl.camera_position = camera_position

    pl.screenshot(str(out_path))
    pl.close()
    return out_path


def gif_pv(
    states: Sequence[State],
    out_gif: str | Path,
    *,
    show_density: bool = True,
    show_mts: bool = True,
    sphere_opacity: float = 0.55,
    tube_radius: float = 0.08,
    window_size: tuple = (900, 800),
    frame_rate: int = 8,
) -> Path:
    """Render a sequence of States into an animated GIF using PyVista."""
    out_gif = Path(out_gif)
    out_gif.parent.mkdir(parents=True, exist_ok=True)

    pl = _base_plotter(off_screen=True, window_size=window_size)
    pl.open_gif(str(out_gif), fps=frame_rate)

    for st in states:
        pl.clear()
        _populate_plotter(pl, st, show_density=show_density, show_mts=show_mts,
                          sphere_opacity=sphere_opacity, tube_radius=tube_radius)
        pl.write_frame()

    pl.close()
    return out_gif


def live_viewer(
    states: Sequence[State],
    *,
    show_density: bool = True,
    show_mts: bool = True,
    sphere_opacity: float = 0.55,
    tube_radius: float = 0.08,
    window_size: tuple = (1000, 900),
    frame_rate: int = 10,
):
    """
    Step through a list of States in an interactive PyVista window.
    Press Q to quit, Space to pause/resume (basic loop).

    Usage in notebook:
        from mttoy.pyvista_viz import live_viewer
        live_viewer(out["snapshots"])
    """
    pl = _base_plotter(off_screen=False, window_size=window_size)
    pl.show(interactive_update=True, auto_close=False)

    import time
    delay = 1.0 / max(frame_rate, 1)

    for st in states:
        pl.clear()
        _populate_plotter(pl, st, show_density=show_density, show_mts=show_mts,
                          sphere_opacity=sphere_opacity, tube_radius=tube_radius)
        pl.update()
        time.sleep(delay)

    pl.close()
