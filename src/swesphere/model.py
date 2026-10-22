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
                 n_pole_rows: int = 0, forcing=None, a: float = 6.371e6, Omega: float = 7.292e-5, g: float = 9.81):
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

    def set_coriolis(self, f_cor):
        """Replace the Coriolis field, in place, and return ``self``.

        A thin wrapper over :meth:`SphereGrid.with_coriolis`, there because
        the test cases that need another ``f`` are reached through the model.
        Passing ``None`` restores the Earth's ``f = 2 Omega sin(phi)``.
        """
        if f_cor is None:
            f_cor = 2.0 * self.grid.Omega * self.grid.sinL
        self.grid = self.grid.with_coriolis(f_cor)
        return self

    def interior_mask(self):
        """Boolean over one field: True outside the polar sponge (the estimated points)."""
        rows = np.repeat(np.arange(self.Nlat), self.Nlon)
        return (rows >= self.n_pole_rows) & (rows < self.Nlat - self.n_pole_rows)

    def interior(self, field):
        """The rows of ``field`` outside the polar sponge.

        Use this rather than ``field[m.n_pole_rows:-m.n_pole_rows]``. That
        idiom was correct while the sponge was always on, and it silently
        returns an EMPTY array when ``n_pole_rows`` is 0, because ``-0`` is
        ``0`` and ``field[0:0]`` selects nothing. Since 0 is the default from
        version 0.6.0, the idiom is now a trap, and one that fails quietly:
        a mean over an empty array is nan, and a standard deviation is nan,
        so the first sign of it is usually a plot that is blank.
        """
        k = self.n_pole_rows
        return field if k == 0 else field[k:self.Nlat - k]

    # ---- dynamics
    def rhs(self, u, v, h):
        du, dv, dh = _rhs(u, v, h, self.grid)
        if self.forcing is not None:
            du, dv, dh = self.forcing(u, v, h, du, dv, dh)
        return du, dv, dh

    def propagate(self, x0, T, just_final_state=True):
        """Integrate the state vector from ``T[0]`` to ``T[-1]``, in seconds.

        The interval must be a whole number of time steps; anything else
        raises, rather than silently returning a state at a different time.
        A zero-length interval returns the state unchanged.

        The spherical-harmonic filter fires on **absolute model time**, every
        ``filter_every`` steps counted from ``t = 0``, not every
        ``filter_every`` steps of this call. Splitting a forecast into
        several calls therefore gives the same answer as one call, which is
        what an external driver (a data assimilation loop, say) needs: the
        length of its forecast window is its own choice and must not change
        the model. Anchoring the schedule to time rather than to a counter on
        the model also keeps successive calls independent, so propagating the
        members of an ensemble through one model object is safe.

        With ``just_final_state=False`` the whole trajectory is returned, as
        an array of shape ``(n_steps + 1, dim)`` whose first row is ``x0``.
        That array costs ``8 * (n_steps + 1) * dim`` bytes, which at the
        default truncation is about 200 kB per step.
        """
        t0, t1 = float(T[0]), float(T[-1])
        span = t1 - t0
        if span < 0.0:
            raise ValueError(
                f"propagate cannot integrate backwards: T[0]={t0:g} is after "
                f"T[-1]={t1:g}."
            )
        n_steps = int(round(span / self.dt))
        if abs(n_steps * self.dt - span) > 1e-6 * self.dt:
            raise ValueError(
                f"the interval {span:g} s is not a whole number of time steps "
                f"(dt={self.dt:g} s). Ask for {n_steps * self.dt:g} s or "
                f"{(n_steps + 1) * self.dt:g} s, or build the model with a dt "
                f"that divides the interval you need."
            )

        u, v, h = (f.copy() for f in self.unpack(np.asarray(x0, dtype=float)))
        traj = None
        if not just_final_state:
            traj = np.empty((n_steps + 1, self.dim), dtype=float)
            traj[0] = self.pack(u, v, h)

        # Index of each step on the absolute time grid, so the filter
        # schedule does not depend on where this call starts.
        k0 = int(round(t0 / self.dt))
        for step in range(1, n_steps + 1):
            u, v, h = self._step(self.rhs, u, v, h, self.dt)
            u = u * self._sponge; v = v * self._sponge
            if (k0 + step) % self.filter_every == 0:
                u = sh_filter(u, self.grid, alpha=self.filter_alpha, p=self.filter_p)
                v = sh_filter(v, self.grid, alpha=self.filter_alpha, p=self.filter_p)
                h = sh_filter(h, self.grid, alpha=self.filter_alpha, p=self.filter_p)
            if not np.isfinite(h).all():
                raise RuntimeError(f"SWE integration blew up at step {step}/{n_steps} (dt={self.dt}); reduce dt or filter more often")
            if traj is not None:
                traj[step] = self.pack(u, v, h)
        return self.pack(u, v, h) if traj is None else traj

    # ---- initial conditions
    def initial_condition(self, kind: str = "rossby", seed: int = 0, **kw):
        if kind == "rossby":
            perturbed = kw.pop("perturbed", False)
            u, v, h = ics.perturbed_ic(self.grid, seed=seed, **kw) if perturbed else ics.default_ic(self.grid, **kw)
        elif kind == "tc2":
            u, v, h = ics.williamson_tc2(self.grid, **kw)
        elif kind in ("rossby_haurwitz", "tc6"):
            u, v, h = ics.rossby_haurwitz(self.grid, **kw)
        elif kind == "galewsky":
            u, v, h = ics.galewsky_jet(self.grid, **kw)
        else:
            raise ValueError(kind)
        return self.pack(u, v, h)
