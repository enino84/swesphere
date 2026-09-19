# -*- coding: utf-8 -*-
"""Named configurations whose behaviour is documented and tested (Sept 2026).

``waves``
    Rossby-Haurwitz initial condition, RK4, no forcing. The initial eddies
    decay in ~40 days and a stationary zonal state with weak Rossby waves
    remains (h std ~380 m, eddy std ~20 m, |U| ~36 m/s). Advection strong,
    error growth slow. Spin-up 40 days; climatology snapshots every 3 days
    are decorrelated (corr 0.13).

``two_jets``
    Galewsky jets at +-45 deg, RK4, zonal relaxation (tau = 5 d) toward the
    jets. Barotropic instability sustains irregular eddies: h anomaly std
    ~130-160 m, |U|max ~70 m/s, stationary from day 30. Advection strong,
    error large and advected (a balanced perturbation travels with the mean
    wind: 8.7 cells in 24 h against 9.1 by pure advection). Spin-up 40 days;
    snapshots every 2 days (corr 0.58).
"""
from __future__ import annotations

from .model import SWEModel
from .forcing import ZonalRelaxation

DAY = 86400.0


def waves(**overrides) -> SWEModel:
    kw = dict(dt=120.0, scheme="rk4"); kw.update(overrides)
    return SWEModel(**kw)


def one_jet(tau_days: float = 5.0, h0: float = 10000.0, **overrides):
    """Galewsky jet in the north only: sustained but fairly regular eddies (wavenumber ~6,
    eddy std ~110 m), the southern hemisphere quiescent. Returns ``(model, x0)``."""
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
