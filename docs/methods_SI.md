## Supporting Information: Simulation Model and Methods

### S1. Overview of stochastic model of microtubule nucleation and growth at condensate interfaces

We developed a coarse-grained stochastic model of microtubule (MT) nucleation, growth, and detachment at the surface of a spherical condensate droplet of radius $R$. The model treats tubulin as a continuous surface density field $n(\mathbf{u}, t)$ defined on a 2D equirectangular grid mapped onto the droplet surface, and represents each MT as a chain of discrete nodes in 3D space. The condensate interface viscosity $\eta_\mathrm{eff}$ is the primary control parameter, encoding the density and cohesive strength of IDP–tubulin contacts at the droplet surface.

![](/Users/potoyan/Dropbox/minima_MT2/schematic.png)


**Figure S1. Reaction scheme for the condensate-surface MT growth model.** A spherical condensate droplet (blue, radius $R$) is surrounded by bulk solution containing free tubulin at concentration $c_\mathrm{bulk}$ with 3D diffusivity $D_\mathrm{3D}$. Tubulin adsorbs from bulk onto the droplet surface to form a 2D density field $n(\mathbf{u},t)$ (blue dots) at rate $k_\mathrm{ads} = k_0^\mathrm{ads} D_\mathrm{3D}/R$ (purple arrow), and redistributes laterally by surface diffusion with coefficient $D_\mathrm{surf} = 1/\eta_\mathrm{eff}$ (cyan double-headed arc), where $\eta_\mathrm{eff}$ is the interface viscosity set by the IDP sequence. When local surface density exceeds a threshold, MT nucleation occurs (green star) at a rate proportional to $D_\mathrm{surf}$. Surface-bound MTs (blue filament) grow along the droplet surface at rate $k_\mathrm{grow}^\mathrm{bound} = k_+^\mathrm{surf} \cdot D_\mathrm{surf}/(D_\mathrm{surf}+D_0)$, which is directly throttled by interface viscosity through the diffusive-saturation factor. Tip detachment (peeling) occurs via Kramers escape at rate $k_\mathrm{off}(L) = k_0^\mathrm{off},e^{-\eta_\mathrm{eff}},\sigma!\left((L/R - \ell_\mathrm{thr})/w_\mathrm{thr}\right)$ (dashed red arrow), where the exponential Kramers factor exponentially suppresses peeling at high $\eta_\mathrm{eff}$ and the sigmoidal gate $\sigma$ captures length-dependent bending stress at the contact zone. Once detached, free tips (orange filament) polymerise from bulk at rate $k_\mathrm{grow}^\mathrm{free} = k_0^\mathrm{on} D_\mathrm{3D} c_\mathrm{bulk}$, independent of surface viscosity. The interplay between diffusion-limited surface growth and Kramers-controlled detachment gives rise to viscosity-dependent MT morphologies ranging from stalled surface-confined filaments (high $\eta_\mathrm{eff}$, RGRGG-like) to radial asters (low $\eta_\mathrm{eff}$, KGKGG-like). Filled circles mark MT anchor points fixed to the droplet surface.

**Surface Tubulin adosprtion.** Surface tubulin is assumed to adsorb from bulk according to Langmuir isotherm kinetics. The equilibrium surface density is

$$
n_\mathrm{eq} = n_\mathrm{max}\,\frac{c_\mathrm{bulk}}{K_d + c_\mathrm{bulk}},
$$

where $n_\mathrm{max}$ is the maximum surface saturation density, $c_\mathrm{bulk}$ is the bulk tubulin concentration, and $K_d$ is the half-saturation constant.

**Surface dynamics** The surface field evolves according to

$$
\frac{\partial n}{\partial t} = k_\mathrm{ads}\!\left(n_\mathrm{eq} - n\right) + D_\mathrm{surf}\,\nabla^2 n,
$$

where the adsorption rate from bulk follows a Smoluchowski flux,

$$
k_\mathrm{ads} = \frac{k_0^\mathrm{ads}\,D_\mathrm{3D}}{R},
$$

and the 2D surface diffusivity obeys a Stokes–Einstein-like relation,

$$
D_\mathrm{surf} = \frac{1}{\eta_\mathrm{eff}}.
$$

High interface viscosity therefore simultaneously slows lateral redistribution of surface tubulin and reduces the adsorption flux from bulk. The PDE is integrated numerically using an exact exponential relaxation for the adsorption term and an explicit finite-difference scheme for the diffusion term, with adaptive substeps to satisfy the stability condition $D_\mathrm{surf}\,\Delta t / \Delta x^2 \leq 0.2$. MT dynamics are simulated as a continuous-time Gillespie (next-reaction) process with three classes of events: nucleation, growth, and tip detachment (peeling).

**Nucleation.** The global nucleation rate is

$$
k_\mathrm{nuc} = \frac{k_0^\mathrm{nuc}\,D_\mathrm{surf}}{R^2}\sum_{i,j}\left\langle\left(\frac{n_{ij} - n_\mathrm{thresh}}{n_\mathrm{thresh}}\right)^{n_\mathrm{Hill}}\right\rangle_+,
$$

where $\langle\cdot\rangle_+$ denotes the positive part, $n_\mathrm{thresh}$ is the local density threshold required for productive nucleation, and $n_\mathrm{Hill}$ is the cooperativity exponent. The $D_\mathrm{surf}$ prefactor reflects that nucleation requires lateral assembly of a tubulin cluster, which is diffusion-limited. When a nucleation event fires, the anchor and initial tip are placed at the selected surface patch and a fixed amount $\Delta n_\mathrm{nuc}$ of tubulin is consumed locally.

**zgrowth for surface-bound tips.** The growth rate for a tip whose last node lies on the droplet surface is

$$
k_\mathrm{grow}^\mathrm{bound} = k_+^\mathrm{surf}\cdot\frac{D_\mathrm{surf}}{D_\mathrm{surf} + D_0},
$$

where $k_+^\mathrm{surf}$ is the intrinsic on-rate and $D_0$ sets the crossover between diffusion-limited and reaction-limited regimes. This saturation form is the steady-state solution of a 2D diffusion-limited capture problem: the flux of tubulin delivered to a moving tip on the surface scales as $D_\mathrm{surf}/(D_\mathrm{surf}+D_0)$. High $\eta_\mathrm{eff}$ reduces $D_\mathrm{surf}$, directly throttling the growth rate. Each successful growth event additionally consumes $\Delta n_\mathrm{grow}$ from the local surface patch; if the patch is depleted below this amount, the event is rejected (kinetic arrest).

**Growth for free tips.** Once a tip has detached from the surface it polymerises directly from bulk:

$$
k_\mathrm{grow}^\mathrm{free} = k_0^\mathrm{on}\,D_\mathrm{3D}\,c_\mathrm{bulk}.
$$

Free tips are not throttled by surface viscosity.

**Simplified microtubule mechanics for free and bound tips** At each growth step the tip advances by a fixed arc length $\delta s$ along the surface (bound) or in 3D (free). The new direction is

$$
\hat{d} = \hat{t} + \sigma\,\boldsymbol{\xi}, \qquad \boldsymbol{\xi}\sim\mathcal{N}(0,\mathbf{I}),
$$

normalised and (for bound tips) projected back onto the sphere. The angular noise amplitude for bound tips is set by worm-like chain statistics,

$$
\sigma_\mathrm{bound} = \sqrt{\delta s / L_p},
$$

where $L_p$ is the MT persistence length. Free tips use a fixed noise $\sigma_\mathrm{free}$ that controls the angular spread of the aster.

**Tip detachment**

The rate at which a surface-bound tip detaches follows Kramers' escape theory modified by a length-dependent force gate (Bell–Evans model):

$$
k_\mathrm{off}(L) = k_0^\mathrm{off}\,e^{-\eta_\mathrm{eff}}\,\sigma\!\left(\frac{L/R - \ell_\mathrm{thresh}}{w_\mathrm{thresh}}\right),
$$

where $\sigma(\cdot)$ is the logistic function, $L$ is the current MT contour length, and $\ell_\mathrm{thresh}$, $w_\mathrm{thresh}$ set the length scale and width of the detachment gate. The factor $e^{-\eta_\mathrm{eff}}$ is the Kramers barrier: denser, more cohesive interfaces (high $\eta_\mathrm{eff}$, corresponding to RGRGG-type sequences) exponentially suppress tip peeling, while low-viscosity interfaces (KGKGG-type) allow spontaneous detachment. Long MTs accumulate bending stress at the contact zone, progressively lowering the effective barrier and triggering peeling.



### S2. Parameter Table

All simulations in this work used the parameter values listed below unless otherwise noted. Varied parameters ($\eta_\mathrm{eff}$, $c_\mathrm{bulk}$, $R$) are indicated in each figure panel.

| Symbol | Description | Value |
|---|---|---|
| **Primary** | | |
| $\eta_\mathrm{eff}$ | Interface viscosity (dimensionless) | varied (0.1–10) |
| $c_\mathrm{bulk}$ | Bulk tubulin concentration (a.u.) | varied (2–20) |
| $R$ | Droplet radius (μm) | 8.0 |
| **Surface field** | | |
| $n_\mathrm{max}$ | Maximum surface density (a.u.) | 80 |
| $K_d$ | Langmuir half-saturation (a.u.) | 5.0 |
| $n_0 / n_\mathrm{eq}$ | Initial surface loading fraction | 0.45 |
| $N_x \times N_y$ | Surface grid resolution | 40 × 40 |
| $L_x \times L_y$ | Grid physical size (μm) | 20 × 20 |
| **Bulk transport** | | |
| $D_\mathrm{3D}$ | Bulk diffusivity (a.u.) | 1.0 |
| $k_0^\mathrm{ads}$ | Adsorption prefactor | 0.14 / 8.0 |
| $k_0^\mathrm{on}$ | Free-tip on-rate prefactor | 3 × 10⁻³ |
| **Nucleation** | | |
| $k_0^\mathrm{nuc}$ | Nucleation prefactor | 2 × 10⁻⁴ |
| $n_\mathrm{thresh}$ | Surface density nucleation threshold | 6.0 |
| $n_\mathrm{Hill}$ | Hill cooperativity exponent | 3.0 |
| $\Delta n_\mathrm{nuc}$ | Tubulin consumed per nucleation event | 8.0 |
| **Growth** | | |
| $k_+^\mathrm{surf}$ | Surface on-rate (s⁻¹) | 0.11 |
| $D_0$ | Diffusion-saturation crossover scale | 4.0 |
| $\Delta n_\mathrm{grow}$ | Tubulin consumed per growth step | 5.0 |
| $\delta s$ | Growth step length (μm) | 0.18 |
| $L_p$ | MT persistence length (μm) | 1000 |
| $\sigma_\mathrm{free}$ | Angular noise for free tips (rad) | 0.05 |
| **Tip detachment** | | |
| $k_0^\mathrm{off}$ | Basal detachment rate (s⁻¹) | 0.005 |
| $\ell_\mathrm{thresh}$ | Length gate threshold ($L/R$) | 0.30 |
| $w_\mathrm{thresh}$ | Length gate width ($L/R$) | 0.30 |
| **Numerics** | | |
| $N_\mathrm{MT}^\mathrm{max}$ | Maximum MTs per droplet | 28 |
| $T_\mathrm{end}$ | Simulation duration (s) | 2000 |
| Seeds | Independent replicates per condition | 3 (42, 43, 44) |

---

## S3. Simulation Algorithm

At each step of the Gillespie loop:

1. Compute all event rates: one nucleation rate $k_\mathrm{nuc}$; per-MT growth rates $k_\mathrm{grow}^{(i)}$ and detachment rates $k_\mathrm{off}^{(i)}$ for each MT $i$.
2. Draw the waiting time $\Delta t \sim \mathrm{Exp}(1/k_\mathrm{tot})$, where $k_\mathrm{tot} = \sum_i$ (all rates).
3. Advance the surface field by $\Delta t$ (adsorption + diffusion substeps).
4. Select and execute one event with probability proportional to its rate.
5. Record observables at pre-specified sample times.

Morphological observables reported are: mean MT contour length $\langle L \rangle$, escape fraction $1 - f_\mathrm{bound}$ (fraction of MTs whose tip has detached from the surface), and the aster order parameter

$$
S = \left\langle \hat{t}_\mathrm{tip} \cdot \hat{n}_\mathrm{tip} \right\rangle,
$$

where $\hat{t}_\mathrm{tip}$ is the tip growth direction and $\hat{n}_\mathrm{tip}$ is the outward radial unit vector. $S \to 1$ indicates a radial aster; $S \to 0$ indicates tangential (tactoid-like) organisation.
