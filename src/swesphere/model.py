# -*- coding: utf-8 -*-
"""The shallow-water model on the sphere with a state-vector API.

State layout: ``x = [u.ravel(), v.ravel(), h.ravel()]`` on the ``(Nlat, Nlon)``
grid, row-major, rows from north to south. ``field_size = Nlat * Nlon``.
"""
from __future__ import annotations

import numpy as np

from .grid import SphereGrid, sh_filter, make_sponge
from .dynamics import rhs as _rhs
from .integrators import INTEGRATORS
from . import initial_conditions as ics


class SWEModel:
    def __init__(self, LMAX: int = 32, dt: float = 120.0, scheme: str = "rk4",
                 filter_every: int = 3, filter_alpha: float = 36.0, filter_p: int = 16,
                 n_pole_rows: int = 10, forcing=None, a: float = 6.371e6, Omega: float = 7.292e-5, g: float = 9.81):
        self.grid = SphereGrid.make(LMAX=LMAX, a=a, Omega=Omega, g=g)
        self.dt, self.scheme = float(dt), scheme
        self.filter_every, self.filter_alpha, self.filter_p = int(filter_every), float(filter_alpha), int(filter_p)
        self.n_pole_rows = int(n_pole_rows)
        self._sponge = make_sponge(self.grid.Nlat, n_pole_rows=n_pole_rows)
        self.forcing = forcing
        self._step = INTEGRATORS[scheme]

    # ---- layout
    @property
    def Nlat(self): return self.grid.Nlat
    @property
    def Nlon(self): return self.grid.Nlon
    @property
    def field_size(self): return self.grid.Nlat * self.grid.Nlon
    @property
    def dim(self): return 3 * self.field_size
    @property
    def var_blocks(self):
        fs = self.field_size
        return {"u": slice(0, fs), "v": slice(fs, 2 * fs), "h": slice(2 * fs, 3 * fs)}

    def pack(self, u, v, h):
        return np.concatenate([np.asarray(u).ravel(), np.asarray(v).ravel(), np.asarray(h).ravel()])

    def unpack(self, x):
        fs = self.field_size; s = (self.Nlat, self.Nlon)
        return x[:fs].reshape(s), x[fs:2 * fs].reshape(s), x[2 * fs:3 * fs].reshape(s)

    def interior_mask(self):
        """Boolean over one field: True outside the polar sponge (the estimated points)."""
        rows = np.repeat(np.arange(self.Nlat), self.Nlon)
        return (rows >= self.n_pole_rows) & (rows < self.Nlat - self.n_pole_rows)

    # ---- dynamics
    def rhs(self, u, v, h):
        du, dv, dh = _rhs(u, v, h, self.grid)
        if self.forcing is not None:
            du, dv, dh = self.forcing(u, v, h, du, dv, dh)
        return du, dv, dh

    def propagate(self, x0, T, just_final_state=True):
        if not just_final_state:
            raise NotImplementedError
        u, v, h = (f.copy() for f in self.unpack(np.asarray(x0, dtype=float)))
        n_steps = max(1, int(round((float(T[-1]) - float(T[0])) / self.dt)))
        for step in range(1, n_steps + 1):
            u, v, h = self._step(self.rhs, u, v, h, self.dt)
            u = u * self._sponge; v = v * self._sponge
            if step % self.filter_every == 0:
                u = sh_filter(u, self.grid, alpha=self.filter_alpha, p=self.filter_p)
                v = sh_filter(v, self.grid, alpha=self.filter_alpha, p=self.filter_p)
                h = sh_filter(h, self.grid, alpha=self.filter_alpha, p=self.filter_p)
            if not np.isfinite(h).all():
                raise RuntimeError(f"SWE integration blew up at step {step}/{n_steps} (dt={self.dt}); reduce dt or filter more often")
        return self.pack(u, v, h)

    # ---- initial conditions
    def initial_condition(self, kind: str = "rossby", seed: int = 0, **kw):
        if kind == "rossby":
            u, v, h = ics.default_ic(self.grid, **kw) if not kw.get("perturbed") else ics.perturbed_ic(self.grid, seed=seed, **kw)
        elif kind == "galewsky":
            u, v, h = ics.galewsky_jet(self.grid, **kw)
        else:
            raise ValueError(kind)
        return self.pack(u, v, h)
