"""mttoy: minimal microtubule-on-condensate growth model.

Gillespie/PDMP simulation of MT polymerization at a condensate interface.
Captures tubulin starvation-limited stalling (high-η condensates, e.g. RGRGG)
vs rapid surface-hugging or aster formation (low-η condensates, e.g. KGKGG).
"""

from .core import Params, Droplet, MT, State, run_simulation
from .analysis import phase_map, run_point
from .viz import render_state, make_gif, gif_from_states, save_snapshots, draw_schematic

__all__ = [
    "Params", "Droplet", "MT", "State", "run_simulation",
    "phase_map", "run_point",
    "render_state", "make_gif", "gif_from_states", "save_snapshots", "draw_schematic",
]
