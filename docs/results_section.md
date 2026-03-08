# Results

## A minimal coarse-grained model recapitulates condensate-directed MT assembly

To interpret how condensate material properties encode distinct MT assembly outcomes,
we developed a coarse-grained stochastic model of MT nucleation and growth at a spherical
condensate surface (Fig. S1; full description in SI Section S1–S5).
The model represents tubulin as a 2D surface density field $n(\mathbf{u},t)$ that
replenishes from bulk via adsorption and redistributes laterally by surface diffusion.
Microtubule dynamics — nucleation, growth, and tip detachment — are governed by
Gillespie kinetics with rates coupled to the local surface field.
The condensate interface viscosity $\eta_\mathrm{eff}$ is the single parameter encoding
sequence identity: it simultaneously sets the 2D surface diffusivity
$D_\mathrm{surf} = 1/\eta_\mathrm{eff}$ (lateral tubulin mobility) and the Kramers
barrier $e^{-\alpha\eta_\mathrm{eff}}$ controlling tip detachment.
High $\eta_\mathrm{eff}$ therefore corresponds to dense, cohesive interfaces
characteristic of RGRGG-type sequences, while low $\eta_\mathrm{eff}$ represents
fluid, weakly adhesive interfaces characteristic of KGKGG-type sequences.



Figure 1. Condensate interface viscosity and tubulin concentration jointly control microtubule assembly morphology and kinetics.

(A) Simulation snapshots at $t = 2000,\mathrm{s}$ illustrating four assembly regimes. Blue filaments have tips in contact with the condensate surface; orange filaments have escaped and grow radially into the bulk. From left to right: Stalled (high $\eta_\mathrm{eff}$, high $c$)  short, surface-arrested filaments with no escape, representative of dense RGRGG-type interfaces; Tactoid (low $\eta_\mathrm{eff}$, low $c$) — tangentially ordered surface filaments; Mixed (low $\eta_\mathrm{eff}$, intermediate $c$) — coexistence of surface-bound and escaped filaments; Aster (low $\eta_\mathrm{eff}$, high $c$) — radially organised, predominantly free MTs characteristic of fluid KGKGG-type interfaces.

(B) Mean MT length over time across four interface viscosities at fixed $c_\mathrm{bulk}$ ($n = 3$ replicates, shading $\pm$ SEM). High $\eta_\mathrm{eff}$ stalls elongation; low $\eta_\mathrm{eff}$ supports rapid growth.

(C) Surface escape fraction ($1 - f_\mathrm{bound}$) over time for the same conditions. High-viscosity interfaces pin tips to the surface while low-viscosity interfaces permit progressive filament release into the bulk.

(D) Total MT polymer length across $N$ condensate droplets at fixed total condensate volume ($n = 3$, error bars $\pm$ SEM). Fragmentation into smaller droplets increases polymer output by raising the surface-to-volume ratio, predicting that mechanical condensate remodelling by growing MTs amplifies net polymerization.

(E) Aster order parameter $S$ (green) and mean MT length (purple, dashed) versus bulk tubulin concentration at low $\eta_\mathrm{eff}$ ($n = 3$, error bars $\pm$ SEM). Increasing $c_\mathrm{bulk}$ drives a cooperative transition from tangential tactoid organisation to radial aster geometry, with filament elongation and reorientation occurring together.

## Interface viscosity governs MT growth kinetics and surface escape

Simulations across a four-decade range of $\eta_\mathrm{eff}$ at fixed bulk tubulin
concentration ($c_\mathrm{bulk} = 10$ a.u.) reveal a clear viscosity-dependent
divergence in MT assembly kinetics (Fig. 1B–C).
At low viscosity ($\eta_\mathrm{eff} = 0.1$), surface-bound tips grow rapidly —
the diffusive-saturation factor $D_\mathrm{surf}/(D_\mathrm{surf}+D_0)$ approaches
unity — and long MTs accumulate quickly before peeling from the surface and escaping
radially into the bulk, producing a high escape fraction (Fig. 1C).
At high viscosity ($\eta_\mathrm{eff} = 10$), the same factor is suppressed by two
orders of magnitude, stalling surface growth, while the Kramers barrier simultaneously
pins tips to the interface and prevents detachment.
The result is a population of short, surface-arrested filaments that neither grow nor
escape — a stalled state.
These two regimes are consistent with the experimentally observed contrast between
KGKGG condensates, which support vigorous MT nucleation and elongation, and denser
RGRGG-like condensates, where filament dynamics are markedly reduced
(Srinivasan et al., this work).
The intermediate viscosity conditions ($\eta_\mathrm{eff} = 0.5$–$3$) produce mixed
phenotypes with moderate growth rates and partial escape, underscoring that the
condensate sequence tunes assembly outcomes continuously rather than through a binary switch.

## Tubulin concentration drives a tactoid-to-aster morphology transition

At low interface viscosity, increasing bulk tubulin concentration drives a sharp
transition in MT organization from tangentially ordered tactoids to radial asters
(Fig. 1A, E).
At low $c_\mathrm{bulk}$, the sparse surface field supports few nucleation events
and MTs grow tangentially along the droplet surface, producing a low aster order
parameter $S \approx 0$.
As $c_\mathrm{bulk}$ increases, both the nucleation rate and the surface replenishment
flux grow, enabling rapid filament elongation and eventual escape; MT tips detach and
grow radially outward, raising $S$ toward 1 and mean MT length simultaneously
(Fig. 1E).
The co-occurrence of increasing length and order mirrors the experimental observation
that enhanced tubulin self-assembly produces long filaments with progressively more
pronounced radial organisation, and that condensate-directed polymerization can switch
between bundled tactoid-like and star-like aster morphologies depending on
tubulin availability.

## Condensate fragmentation non-monotonically regulates total polymer mass

A key prediction of the transport-limited growth mechanism is that droplet size and
number matter: at fixed total condensate volume, fragmenting a single droplet into
$N$ smaller ones of radius $R = R_0 / N^{1/3}$ modulates the effective adsorption
flux $k_\mathrm{ads} \propto 1/R$ and nucleation geometry.
Simulations show that total MT polymer length (summed across all $N$ droplets)
increases sharply as $N$ rises from 1 to $\sim$8 and then saturates or slightly
decreases at large $N$ (Fig. 1D).
Smaller droplets present a larger surface-to-volume ratio, enhancing per-droplet tubulin
recruitment and nucleation, but at very large $N$ each droplet is too small to sustain
long filaments before they geometrically outgrow the surface.
This non-monotonic dependence on fragmentation provides a testable experimental
prediction: pelleting or turbidity assays on mechanically fragmented condensate
emulsions should show a peak in total polymer mass at intermediate droplet number,
consistent with the observation that growing MTs mechanically destabilize and
rearrange the condensate network (Srinivasan et al., this work), which
would effectively increase $N$ over time.

## Sequence-encoded material properties as tunable regulatory knobs

Taken together, the simulation results establish a mechanistic picture in which
condensate sequence identity controls MT assembly through two orthogonal channels
encoded in $\eta_\mathrm{eff}$: diffusion-limited surface feeding sets the growth
rate, and Kramers-barrier-controlled adhesion sets the lifetime of surface contact.
Low-viscosity, fluid condensates (KGKGG-type) are efficient MT catalysts — they
rapidly deliver tubulin to growing tips and release filaments into the bulk to form
asters — while high-viscosity, cohesive condensates (RGRGG-type) trap and stall MTs
at the interface.
The concentration of free tubulin acts as an independent axis, pushing the system
from tactoid to aster organisation above a threshold set by nucleation cooperativity.
These predictions are in quantitative agreement with the experimental finding that
condensate material properties, achieved through modular peptide and nucleic acid
design, act as programmable regulatory forces over MT assembly, and that the same
condensate platform can be tuned to produce qualitatively distinct cytoskeletal
architectures.
