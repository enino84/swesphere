# -*- coding: utf-8 -*-
"""Named configurations. The numbers that characterize each one are measured by the experiments
(``results/EXP-01/stability.csv``, ``results/EXP-02/regimes.csv``, ``results/EXP-03/footprint.csv``)
for the default ``dt = 120 s`` with the filter every 3 steps; they change with the filter cadence.

``waves``
    Solid-body rotation (Williamson test case 2) plus a wavenumber-4 wave packet in the north and a
    wavenumber-3 packet in the south, RK4, no forcing. The initial eddies decay in ~40 days and zonal
    flow with weak Rossby waves remains. With no source of energy the waves keep decaying slowly (the
    spectral filter is the only sink): the regime is quasi-linear and slowly decaying, not stationary.
    Spin-up 40 days; snapshots every 3 days.

``one_jet``
    Galewsky jet at 45 N, RK4, relaxation of the zonal mean (tau = 5 d) toward the jet. Sustained and
    fairly regular eddies in the north, the south quiescent. Spin-up 40 days; snapshots every 2 days.

``two_jets``
    Galewsky jets at +-45 deg, same relaxation. Barotropic instability sustains irregular eddies in both
    hemispheres: a stationary turbulent regime in which balanced perturbations travel with the mean wind
    (EXP-03). Spin-up 40 days; snapshots every 2 days.
"""
from __future__ import annotations

from .model import SWEModel
from .forcing import ZonalRelaxation

DAY = 86400.0


def waves(**overrides) -> SWEModel:
    kw = dict(dt=120.0, scheme="rk4"); kw.update(overrides)
    return SWEModel(**kw)


def one_jet(tau_days: float = 5.0, h0: float = 10000.0, **overrides):
    """Galewsky jet in the north only: sustained, fairly regular eddies; the southern hemisphere quiescent.
    Returns ``(model, x0)``."""
    return _jets(False, tau_days, h0, **overrides)


def two_jets(tau_days: float = 5.0, h0: float = 10000.0, **overrides):
    """Galewsky jets at +-45 deg: irregular sustained turbulence in both hemispheres. Returns ``(model, x0)``."""
    return _jets(True, tau_days, h0, **overrides)


def _jets(both, tau_days, h0, **overrides):
    kw = dict(dt=120.0, scheme="rk4"); kw.update(overrides)
    model = SWEModel(**kw)
    x0 = model.initial_condition("galewsky", h0=h0, both_hemispheres=both)
    u, v, h = model.unpack(x0)
    model.forcing = ZonalRelaxation(u, v, h, tau_days * DAY)
    return model, x0


PRESETS = {"waves": lambda **kw: (waves(**kw), waves(**kw).initial_condition("rossby")), "one_jet": one_jet, "two_jets": two_jets}
SPINUP_DAYS = {"waves": 40.0, "one_jet": 40.0, "two_jets": 40.0}
SNAPSHOT_DAYS = {"waves": 3.0, "one_jet": 2.0, "two_jets": 2.0}


def make(name: str, **overrides):
    """``(model, x0)`` for a preset name."""
    if name == "waves":
        m = waves(**overrides); return m, m.initial_condition("rossby")
    return PRESETS[name](**overrides)
