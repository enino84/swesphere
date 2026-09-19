# -*- coding: utf-8 -*-
"""Global invariants of the shallow-water equations and error norms.

All integrals are area weighted (cos φ) over the whole grid, or over the rows
outside the polar sponge with ``interior=True``. The continuous equations
conserve mass, total energy and potential enstrophy; the discrete model loses a
little of each through the spectral filter and the polar sponge, and the forced
presets exchange them with the relaxation. EXP-00 reports the measured drifts.
"""
from __future__ import annotations

import numpy as np

from .dynamics import relative_vorticity


def area_weights(model, interior: bool = False) -> np.ndarray:
    w = np.clip(np.asarray(model.grid.cosL), 0.0, None).copy()
    if interior:
        k = model.n_pole_rows; w[:k] = 0.0; w[model.Nlat - k:] = 0.0
    return w / w.sum()


def invariants(model, x, interior: bool = False) -> dict:
    """Mean depth (m), total energy per unit area (J/kg m) and potential enstrophy."""
    u, v, h = model.unpack(np.asarray(x, dtype=float)); G = model.grid; w = area_weights(model, interior)
    eta = relative_vorticity(u, v, G) + G.f_cor
    return dict(mass=float((w * h).sum()),
                energy=float((w * (0.5 * h * (u ** 2 + v ** 2) + 0.5 * G.g * h ** 2)).sum()),
                enstrophy=float((w * 0.5 * eta ** 2 / h).sum()))


def error_norms(model, x, x_ref, field: str = "h", interior: bool = True) -> dict:
    """Normalized l1, l2 and l-infinity errors of Williamson et al. (1992)."""
    i = "uvh".index(field); q = model.unpack(np.asarray(x, dtype=float))[i]; r = model.unpack(np.asarray(x_ref, dtype=float))[i]
    w = area_weights(model, interior); m = w > 0
    return dict(l1=float((w * np.abs(q - r)).sum() / (w * np.abs(r)).sum()),
                l2=float(np.sqrt((w * (q - r) ** 2).sum() / (w * r ** 2).sum())),
                linf=float(np.abs(q - r)[m].max() / np.abs(r)[m].max()))
