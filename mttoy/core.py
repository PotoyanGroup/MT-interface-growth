from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Tuple

import numpy as np
import pandas as pd


# ════════════════════════════════════════════════════════════════
# Parameters
# ════════════════════════════════════════════════════════════════

@dataclass
class Params:
    name: str

    # ── Primary knobs ────────────────────────────────────────────
    eta_eff: float = 1.0    # interface viscosity (RGRGG ~ 8, KGKGG ~ 0.6)
    # Dual role: (1) sets D_surf = 1/eta_eff (Stokes-Einstein);
    #            (2) enters k_off as the dimensionless adhesion energy E_adh/kT
    #                (Kramers barrier). Both reflect the same molecular quantity —
    #                the density/strength of IDP-tubulin contacts at the interface.
    R:       float = 8.0    # droplet radius
    c_dilute:  float = 5.0    # dilute-phase tubulin concentration

    # ── Surface tubulin field ────────────────────────────────────
    n_max:        float = 80.0   # surface saturation density (Langmuir)
    Kd_bulk:      float = 5.0    # half-saturation concentration
    Nx:           int   = 40     # slightly finer grid helps local depletion
    Ny:           int   = 40
    Lx:           float = 20.0   # grid physical size
    Ly:           float = 20.0
    init_surface_frac: float = 0.15
    # initial surface loading as a fraction of n_eq
    # < 1 is important: otherwise surface starts already fully fueled everywhere

    # ── Bulk 3D diffusion — unifies adsorption and free-tip growth ───────────
    D_3D:         float = 1.0
    k0_ads:       float = 0.0125  # k_ads = k0_ads · D_3D / R  (Smoluchowski per-area flux)
    k0_on:        float = 2e-5   # k_on_3D = k0_on · D_3D

    # ── Nucleation ───────────────────────────────────────────────
    k0_nuc:       float = 7e-7
    n_hill:       float = 3.0
    n_nuc_thresh: float = 18.0
    dn_nuc:       float = 30.0   # stronger local depletion at nucleation

    # ── Growth ───────────────────────────────────────────────────
    k_plus_surf:  float = 0.07
    D0_grow:      float = 3.5
    # higher D0_grow strengthens the viscosity penalty on bound growth
    dn_grow:      float = 5.0
    # bound growth now truly consumes local tubulin
    ds:           float = 0.18

    # smaller step -> tip stays in a local patch longer -> easier to deplete
    L_p:          float = 450.0
    # MT persistence length (same units as ds, R).
    # Sets bound-tip angular noise via WLC: noise = sqrt(ds/L_p).
    # Default 450 ~ 0.45 mm, consistent with dynamic MT stiffness.
    noise_free:   float = 0.05 #75 microm
    # Angular noise for free (escaped) tips. Larger than bound noise —
    # not MT stiffness but aster-spread: tips grow radially with some
    # directional scatter set by this parameter.

    # ── Tip peeling ──────────────────────────────────────────────
    k0_off:       float = 0.005
    n_thresh:     float = 2.5
    w_thresh:     float = 0.5

    # ── Numerical limits ─────────────────────────────────────────
    max_mts:      int   = 60


# ════════════════════════════════════════════════════════════════
# Derived quantities
# ════════════════════════════════════════════════════════════════

def D_surf(p: Params) -> float:
    """2D surface diffusion: D_surf = 1 / eta_eff."""
    return 1.0 / max(p.eta_eff, 1e-12)


def k_ads(p: Params) -> float:
    """Surface replenishment rate from bulk (Smoluchowski per-area flux: j ~ D·c/R)."""
    return p.k0_ads * p.D_3D / max(p.R, 1e-12)


def k_on_3D(p: Params) -> float:
    """Free-tip bulk polymerization rate constant."""
    return p.k0_on * p.D_3D


def n_eq(p: Params) -> float:
    """Langmuir equilibrium surface density."""
    return p.n_max * p.c_dilute / (p.Kd_bulk + p.c_dilute)


def k_off(p: Params, L: float) -> float:
    """
    Kramers escape rate for MT tip detachment:

        k_off(L) = k0_off · exp(-E_adh/kT) · σ(L)

    where E_adh/kT ≡ eta_eff  (kT = 1 in simulation units).

    Physical basis:
      exp(-eta_eff): Kramers factor. The condensate interfacial viscosity is
        proportional to the tip adhesion free energy — denser / more cohesive
        interfaces (RGRGG, high eta_eff) present a higher barrier to tip
        detachment, exponentially suppressing peeling.

      σ(L) = sigmoid((L/R − n_thresh) / w_thresh): length gate. A long MT
        accumulates elastic bending stress at the contact zone, progressively
        lowering the effective barrier (Bell/Evans force-dependent unbinding).
        Short MTs (L/R < n_thresh) remain pinned; long ones peel spontaneously.
    """
    x = (L / max(p.R, 1e-12) - p.n_thresh) / max(p.w_thresh, 1e-12)
    x = max(-60.0, min(60.0, x))
    gate = 1.0 / (1.0 + math.exp(-x))
    return p.k0_off * math.exp(-p.eta_eff) * gate


def bound_growth_sat(p: Params) -> float:
    """
    Diffusive feeding penalty for surface-bound growth.
    High eta -> small D_surf -> small sat.
    """
    dsurf = D_surf(p)
    return dsurf / (dsurf + p.D0_grow)


def tip_noise(p: Params) -> float:
    """
    Angular noise per growth step for a surface-bound tip (WLC statistics):
        sigma = sqrt(ds / L_p)
    Recovers <theta^2> = ds/L_p per step, consistent with MT semiflexibility.
    """
    return math.sqrt(p.ds / max(p.L_p, 1e-12))



# ════════════════════════════════════════════════════════════════
# State
# ════════════════════════════════════════════════════════════════

@dataclass
class Droplet:
    center: np.ndarray
    R: float
    params: Params


@dataclass
class SurfaceField:
    n: np.ndarray   # shape (Ny, Nx)


@dataclass
class MT:
    nodes:      List[np.ndarray]
    node_bound: List[bool] = field(default_factory=list)
    uv:         Optional[np.ndarray] = None  # 2D grid coordinate; None after tip peels

    def __post_init__(self):
        if not self.node_bound:
            self.node_bound = [False] * len(self.nodes)

    @property
    def tip_bound(self) -> bool:
        return bool(self.node_bound[-1]) if self.node_bound else False

    @property
    def bound(self) -> bool:
        return any(self.node_bound)

    def length(self) -> float:
        if len(self.nodes) < 2:
            return 0.0
        pts = np.vstack(self.nodes)
        return float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))


@dataclass
class State:
    t:       float
    droplet: Droplet
    mts:     List[MT] = field(default_factory=list)
    surface: Optional[SurfaceField] = None
    rng:     np.random.Generator = field(default_factory=lambda: np.random.default_rng(0))


# ════════════════════════════════════════════════════════════════
# Helpers
# ════════════════════════════════════════════════════════════════

def unit(v: np.ndarray) -> np.ndarray:
    return v / (np.linalg.norm(v) + 1e-12)


def sample_tangent(rng: np.random.Generator, normal: np.ndarray) -> np.ndarray:
    r = rng.normal(size=3)
    return unit(r - np.dot(r, normal) * normal)


def prob(k: float, dt: float) -> float:
    if k <= 0 or dt <= 0:
        return 0.0
    return 1.0 - math.exp(-k * dt)


def project_to_sphere(center: np.ndarray, R: float, x: np.ndarray) -> np.ndarray:
    return center + R * unit(x - center)


# ════════════════════════════════════════════════════════════════
# 2D surface grid utilities
# ════════════════════════════════════════════════════════════════

def _uv_wrap(p: Params, uv: np.ndarray) -> np.ndarray:
    x = ((uv[0] + 0.5 * p.Lx) % p.Lx) - 0.5 * p.Lx
    y = ((uv[1] + 0.5 * p.Ly) % p.Ly) - 0.5 * p.Ly
    return np.array([x, y], dtype=float)


def _uv_to_ij(p: Params, uv: np.ndarray) -> Tuple[int, int]:
    uv = _uv_wrap(p, uv)
    j = int(np.floor((uv[0] + 0.5 * p.Lx) / p.Lx * p.Nx)) % p.Nx
    i = int(np.floor((uv[1] + 0.5 * p.Ly) / p.Ly * p.Ny)) % p.Ny
    return i, j


def _ij_to_uv(p: Params, i: int, j: int) -> np.ndarray:
    return np.array([
        (j + 0.5) / p.Nx * p.Lx - 0.5 * p.Lx,
        (i + 0.5) / p.Ny * p.Ly - 0.5 * p.Ly,
    ])


def _sphere_to_uv(p: Params, pos: np.ndarray, center: np.ndarray) -> np.ndarray:
    """
    Equirectangular (plate carrée) projection: 3D sphere position → UV grid coordinate.
        phi   in [-π, π]  →  u in [-Lx/2, Lx/2]
        theta in [0,  π]  →  v in [-Ly/2, Ly/2]
    Exact for any point on the sphere; replaces incremental tang[:2] update.
    """
    r_hat = unit(pos - center)
    theta = np.arccos(np.clip(r_hat[2], -1.0, 1.0))
    phi   = np.arctan2(r_hat[1], r_hat[0])
    u = phi   / np.pi * (0.5 * p.Lx)
    v = (theta / np.pi - 0.5) * p.Ly
    return _uv_wrap(p, np.array([u, v]))


def _uv_to_sphere(p: Params, uv: np.ndarray, center: np.ndarray, R: float) -> np.ndarray:
    """
    Inverse equirectangular projection: UV grid coordinate → 3D sphere position.
    Used in spawn so the nucleation anchor sits on the selected surface patch.
    """
    uv    = _uv_wrap(p, uv)
    phi   = uv[0] / (0.5 * p.Lx) * np.pi
    theta = (uv[1] / p.Ly + 0.5) * np.pi
    r_hat = np.array([
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta),
    ])
    return center + R * r_hat


def surf_consume(st: State, uv: np.ndarray, amount: float) -> bool:
    i, j = _uv_to_ij(st.droplet.params, uv)
    if st.surface.n[i, j] < amount:
        return False
    st.surface.n[i, j] -= amount
    return True


def surf_step(st: State, dt: float) -> None:
    """
    Surface tubulin dynamics:
        dn/dt = k_ads (n_eq - n) + D_surf ∇² n

    Implemented with:
      - exact local relaxation toward n_eq over dt
      - explicit diffusion substeps for numerical stability
    """
    p = st.droplet.params
    n = st.surface.n

    # adsorption / replenishment toward Langmuir equilibrium
    n += (n_eq(p) - n) * (1.0 - np.exp(-k_ads(p) * dt))

    # explicit diffusion with stable substeps
    dsurf = D_surf(p)
    if dsurf > 0 and dt > 0:
        dx = p.Lx / p.Nx
        dy = p.Ly / p.Ny
        dmin2 = min(dx * dx, dy * dy)

        # conservative explicit stability condition in 2D
        dt_stable = 0.20 * dmin2 / max(dsurf, 1e-12)
        n_sub = max(1, int(math.ceil(dt / max(dt_stable, 1e-12))))
        dt_sub = dt / n_sub

        ax = dsurf * dt_sub / (dx * dx)
        ay = dsurf * dt_sub / (dy * dy)

        for _ in range(n_sub):
            lap = (
                ay * (np.roll(n, 1, axis=0) + np.roll(n, -1, axis=0) - 2.0 * n) +
                ax * (np.roll(n, 1, axis=1) + np.roll(n, -1, axis=1) - 2.0 * n)
            )
            n += lap
            np.clip(n, 0.0, p.n_max, out=n)

    np.clip(n, 0.0, p.n_max, out=n)


# ════════════════════════════════════════════════════════════════
# Rates
# ════════════════════════════════════════════════════════════════

def nuc_rate(st: State) -> float:
    """
    Productive surface nucleation:
        k_nuc = k0_nuc · D_surf · <excess^n_hill> / R²
    """
    p = st.droplet.params
    if len(st.mts) >= p.max_mts or st.surface is None:
        return 0.0

    n = st.surface.n
    y = np.maximum(0.0, (n - p.n_nuc_thresh) / max(p.n_nuc_thresh, 1e-12))
    hill = float(np.mean(y ** p.n_hill)) * (p.Nx * p.Ny)
    return p.k0_nuc * D_surf(p) * hill / max(p.R ** 2, 1e-12)


def grow_rate(st: State, mt: MT) -> float:
    """
    Bound tip:
        k_plus_surf * sat
    where sat = D_surf / (D_surf + D0_grow)

    Free tip:
        k_on_3D * c_dilute
    Free tips are no longer throttled by surface viscosity.
    Once detached, they polymerize from bulk.
    """
    p = st.droplet.params
    if mt.tip_bound:
        return p.k_plus_surf * bound_growth_sat(p)
    return k_on_3D(p) * p.c_dilute


# ════════════════════════════════════════════════════════════════
# Events
# ════════════════════════════════════════════════════════════════

def spawn(st: State, dt: float) -> None:
    """
    Nucleate at a surface patch weighted by local tubulin excess.
    Anchor and tip are born bound.
    """
    p = st.droplet.params
    rng = st.rng
    if st.surface is None:
        return

    n = st.surface.n
    y = np.maximum(0.0, (n - p.n_nuc_thresh) / max(p.n_nuc_thresh, 1e-12))
    w = y ** p.n_hill
    wsum = float(np.sum(w))
    if wsum <= 0.0:
        return

    flat = w.ravel()
    idx = int(rng.choice(flat.size, p=flat / wsum))
    uv = _ij_to_uv(p, idx // p.Nx, idx % p.Nx)

    if not surf_consume(st, uv, p.dn_nuc):
        return

    # Place anchor at the 3D point corresponding to the selected UV patch.
    anchor = _uv_to_sphere(p, uv, st.droplet.center, p.R)
    nh = unit(anchor - st.droplet.center)
    dirv = unit(sample_tangent(rng, nh) + tip_noise(p) * rng.normal(size=3))
    tip = project_to_sphere(st.droplet.center, p.R, anchor + p.ds * dirv)

    st.mts.append(MT(nodes=[anchor, tip], node_bound=[True, True],
                     uv=_sphere_to_uv(p, tip, st.droplet.center)))


def grow(st: State, idx: int, dt: float) -> None:
    """
    Extend MT tip by one step.

    Important physical change:
      bound growth now consumes local surface tubulin.
      If the local patch is depleted, the MT stalls.
    """
    p = st.droplet.params
    rng = st.rng
    mt = st.mts[idx]
    d = st.droplet

    # need at least a tangent estimate
    if len(mt.nodes) < 2:
        return
    t_hat = unit(mt.nodes[-1] - mt.nodes[-2])

    # ── bound branch: local consumption + constrained surface walk ──────────
    if mt.tip_bound:
        if mt.uv is None:
            return

        # true kinetic arrest: if local patch lacks tubulin, growth does not occur
        if not surf_consume(st, mt.uv, p.dn_grow):
            return

        dirv = unit(t_hat + tip_noise(p) * rng.normal(size=3))

        new_pos = project_to_sphere(d.center, d.R, mt.nodes[-1] + p.ds * dirv)

        mt.nodes.append(new_pos)
        mt.node_bound.append(True)

        # Recompute UV directly from the new 3D tip position (equirectangular projection).
        # Exact for all sphere positions; eliminates the tang[:2] approximation error.
        mt.uv = _sphere_to_uv(p, new_pos, d.center)

    # ── free branch: persistent walk into bulk ──────────────────────────────
    else:
        dirv = unit(t_hat + p.noise_free * rng.normal(size=3))
        mt.nodes.append(mt.nodes[-1] + p.ds * dirv)
        mt.node_bound.append(False)


def peel(st: State, idx: int) -> None:
    """Detach tip from surface (Kramers escape). Anchor stays; tip goes free."""
    mt = st.mts[idx]
    mt.node_bound[-1] = False
    mt.uv = None


# ════════════════════════════════════════════════════════════════
# Morphometrics
# ════════════════════════════════════════════════════════════════

def compute_morphometrics(st: State) -> Tuple[float, float]:
    """
    f_bound  : fraction of MTs with tip still on surface
    aster_S  : mean radial alignment of tip tangent (0 = tangent, 1 = radial)
    """
    if not st.mts:
        return 0.0, 0.0

    center = st.droplet.center
    n_tip_bound = 0
    aligns = []

    for mt in st.mts:
        if len(mt.nodes) < 2:
            continue
        if mt.tip_bound:
            n_tip_bound += 1
        t_hat = unit(mt.nodes[-1] - mt.nodes[-2])
        n_hat = unit(mt.nodes[-1] - center)
        aligns.append(float(max(0.0, np.dot(t_hat, n_hat))))

    f_bound = n_tip_bound / max(len(st.mts), 1)
    aster_S = float(np.mean(aligns)) if aligns else 0.0
    return f_bound, aster_S


# ════════════════════════════════════════════════════════════════
# Gillespie loop
# ════════════════════════════════════════════════════════════════

def _run_loop(
    st: State,
    T_end: float,
    sample_times: Optional[np.ndarray] = None,
    snap_dt: Optional[float] = None,
) -> Tuple[List[dict], List[State]]:
    rng = st.rng
    records: List[dict] = []
    snaps: List[State] = []
    si = 0
    next_snap = 0.0 if snap_dt is not None else None

    while st.t < T_end:
        rates: List[float] = []
        events: List[Tuple] = []

        rn = nuc_rate(st)
        if rn > 0:
            rates.append(rn)
            events.append(("nuc",))

        for i, mt in enumerate(st.mts):
            rg = grow_rate(st, mt)
            if rg > 0:
                rates.append(rg)
                events.append(("grow", i))
            if mt.tip_bound:
                rp = k_off(st.droplet.params, mt.length())
                if rp > 0:
                    rates.append(rp)
                    events.append(("peel", i))

        total = float(np.sum(rates)) if rates else 0.0

        dt = min(
            rng.exponential(1.0 / total) if total > 0 else 1.0,
            T_end - st.t,
        )

        if st.surface is not None:
            surf_step(st, dt)

        if total > 0:
            r = rng.random() * total
            cum = 0.0
            chosen = None
            for rate, ev in zip(rates, events):
                cum += rate
                if r <= cum:
                    chosen = ev
                    break

            if chosen is not None:
                if chosen[0] == "nuc":
                    spawn(st, dt)
                elif chosen[0] == "grow":
                    grow(st, chosen[1], dt)
                elif chosen[0] == "peel":
                    peel(st, chosen[1])

        st.t += dt

        if sample_times is not None:
            while si < len(sample_times) and st.t >= sample_times[si]:
                _record(st, records, sample_times[si])
                si += 1

        if snap_dt is not None and st.t >= next_snap:
            snaps.append(_copy_state(st))
            next_snap += snap_dt

    return records, snaps


def _record(st: State, records: List[dict], t_sample: float) -> None:
    p = st.droplet.params
    mts = st.mts
    L = [mt.length() for mt in mts]
    n_tip_bound = sum(1 for mt in mts if mt.tip_bound)
    n_escaped = sum(1 for mt in mts if mt.bound and not mt.tip_bound)
    f_bound, aster_S = compute_morphometrics(st)
    mean_surf = float(np.mean(st.surface.n)) if st.surface is not None else 0.0

    records.append(dict(
        condition      = p.name,
        t              = float(t_sample),
        eta_eff        = float(p.eta_eff),
        R              = float(p.R),
        c_dilute       = float(p.c_dilute),
        D_surf         = float(D_surf(p)),
        k_ads          = float(k_ads(p)),
        k_off_at_L0    = float(k_off(p, 0.0)),
        n_mt           = int(len(mts)),
        n_tip_bound    = int(n_tip_bound),
        n_escaped      = int(n_escaped),
        frac_tip_bound = float(n_tip_bound / max(1, len(mts))),
        mean_L         = float(np.mean(L)) if L else 0.0,
        max_L          = float(np.max(L)) if L else 0.0,
        mean_surf      = mean_surf,
        min_surf       = float(np.min(st.surface.n)) if st.surface is not None else 0.0,
        escape_index   = float(1.0 - f_bound),
        aster_S        = float(aster_S),
    ))


def _copy_state(st: State) -> State:
    d = st.droplet
    dd = Droplet(center=d.center.copy(), R=d.R, params=d.params)
    mts = [
        MT(
            nodes=[nd.copy() for nd in mt.nodes],
            node_bound=list(mt.node_bound),
            uv=None if mt.uv is None else mt.uv.copy(),
        )
        for mt in st.mts
    ]
    surf = SurfaceField(n=st.surface.n.copy()) if st.surface is not None else None
    return State(t=st.t, droplet=dd, mts=mts, rng=np.random.default_rng(0), surface=surf)


def _make_state(p: Params, seed: int) -> State:
    rng = np.random.default_rng(seed)
    d = Droplet(center=np.zeros(3), R=p.R, params=p)

    # Start below equilibrium so replenishment matters dynamically.
    n0 = np.full((p.Ny, p.Nx), p.init_surface_frac * n_eq(p), dtype=float)

    return State(t=0.0, droplet=d, rng=rng, surface=SurfaceField(n=n0))


# ════════════════════════════════════════════════════════════════
# Public API
# ════════════════════════════════════════════════════════════════

def simulate(
    params: Params,
    T_end: float = 1800.0,
    seed: int = 0,
    sample_n: int = 46,
) -> Tuple[State, pd.DataFrame]:
    st = _make_state(params, seed)
    records, _ = _run_loop(st, T_end, sample_times=np.linspace(0, T_end, sample_n))
    return st, pd.DataFrame(records)


def snapshot(
    params: Params,
    T_end: float = 1800.0,
    snapshot_dt: float = 120.0,
    seed: int = 0,
) -> List[State]:
    st = _make_state(params, seed)
    _, snaps = _run_loop(st, T_end, snap_dt=snapshot_dt)
    return snaps


def run_simulation(
    params: Params,
    T_end: float = 3600.0,
    seed: int = 0,
    sample_n: int = 61,
    snapshot_dt: Optional[float] = None,
) -> Dict[str, Any]:
    st = _make_state(params, seed)
    records, snaps = _run_loop(
        st,
        T_end,
        sample_times=np.linspace(0, T_end, sample_n),
        snap_dt=snapshot_dt if (snapshot_dt and snapshot_dt > 0) else None,
    )
    return {"final_state": st, "metrics": pd.DataFrame(records), "snapshots": snaps}


