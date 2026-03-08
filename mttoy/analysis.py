from __future__ import annotations

from dataclasses import asdict
from typing import Iterable, Optional, Tuple, Dict, Any, Sequence

import numpy as np
import pandas as pd

from .core import Params, run_simulation


def run_point(
    base_params: Params,
    c_bulk: float,
    eta_eff: float,
    seeds: Sequence[int] = (1, 2, 3),
    T_end: float = 3600.0,
    sample_n: int = 61,
    tail_frac: float = 0.25,   # average over last 25% timepoints
) -> Dict[str, float]:
    vals_escape = []
    vals_aster = []
    vals_maxL = []
    vals_meanL = []

    base = asdict(base_params)
    for sd in seeds:
        d = dict(base)
        d["c_bulk"] = float(c_bulk)
        d["eta_eff"] = float(eta_eff)
        d["name"] = base.get("name", "cond")

        p = Params(**d)
        out = run_simulation(p, T_end=T_end, seed=int(sd), sample_n=sample_n, snapshot_dt=None)
        df = out["metrics"].sort_values("t")

        # use tail-average to reduce all-or-none noise
        k0 = int(np.floor((1.0 - tail_frac) * len(df)))
        tail = df.iloc[max(0, k0):]

        # escape
        if "frac_tip_bound" not in tail.columns:
            raise KeyError("Expected column 'frac_tip_bound' in metrics.")
        vals_escape.append(1.0 - float(tail["frac_tip_bound"].mean()))

        # aster
        if "aster_S" in tail.columns:
            vals_aster.append(float(tail["aster_S"].mean()))

        # length
        if "max_L" in tail.columns:
            vals_maxL.append(float(tail["max_L"].mean()))
        if "mean_L" in tail.columns:
            vals_meanL.append(float(tail["mean_L"].mean()))

    def _sem(x):
        x = np.asarray(x, float)
        if len(x) <= 1:
            return 0.0
        return float(np.std(x, ddof=1) / np.sqrt(len(x)))

    out = {
        "escape_mean": float(np.mean(vals_escape)),
        "escape_sem": _sem(vals_escape),
    }
    if len(vals_aster) > 0:
        out["aster_mean"] = float(np.mean(vals_aster))
        out["aster_sem"] = _sem(vals_aster)
    if len(vals_maxL) > 0:
        out["maxL_mean"] = float(np.mean(vals_maxL))
        out["maxL_sem"] = _sem(vals_maxL)
    if len(vals_meanL) > 0:
        out["meanL_mean"] = float(np.mean(vals_meanL))
        out["meanL_sem"] = _sem(vals_meanL)

    return out


def phase_map(
    base_params: Params,
    c_grid: Iterable[float],
    eta_grid: Iterable[float],
    seeds: Sequence[int] = (1, 2, 3),
    T_end: float = 3600.0,
    sample_n: int = 61,
) -> Dict[str, np.ndarray]:
    c_grid = list(c_grid)
    eta_grid = list(eta_grid)

    shape = (len(eta_grid), len(c_grid))
    escape = np.zeros(shape, float)
    aster  = np.full(shape, np.nan, float)
    maxL   = np.full(shape, np.nan, float)
    meanL  = np.full(shape, np.nan, float)

    for iy, eta in enumerate(eta_grid):
        for ix, c in enumerate(c_grid):
            res = run_point(base_params, c, eta, seeds=seeds, T_end=T_end, sample_n=sample_n)
            escape[iy, ix] = res["escape_mean"]
            if "aster_mean" in res: aster[iy, ix] = res["aster_mean"]
            if "maxL_mean"  in res: maxL[iy, ix]  = res["maxL_mean"]
            if "meanL_mean" in res: meanL[iy, ix] = res["meanL_mean"]

    out = {"escape": escape}
    if not np.all(np.isnan(aster)): out["aster"] = aster
    if not np.all(np.isnan(maxL)):  out["maxL"]  = maxL
    if not np.all(np.isnan(meanL)): out["meanL"] = meanL
    return out
