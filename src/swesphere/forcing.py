# -*- coding: utf-8 -*-
"""Forcing that gives the model a stationary turbulent regime.

Without forcing the shallow-water equations on the sphere do not sustain eddies: the filter dissipates
whatever the initial condition put in and the flow relaxes toward a zonal state with weak, slowly
decaying waves (the ``waves`` preset). ``ZonalRelaxation`` relaxes only the zonal mean of (u, v, h)
toward a reference jet and leaves the eddies alone; with a barotropically unstable reference
(Galewsky) the instability regenerates eddies against the dissipation and the flow reaches a
stationary turbulent state (EXP-02 measures its statistics). Relaxing the full fields toward the
zonal state kills the eddies instead.
"""
from __future__ import annotations

import numpy as np


class ZonalRelaxation:
    def __init__(self, u_ref, v_ref, h_ref, tau):
        self.k = 1.0 / float(tau)
        self.ref = tuple(np.repeat(np.asarray(f).mean(axis=1, keepdims=True), np.asarray(f).shape[1], axis=1) for f in (u_ref, v_ref, h_ref))

    def __call__(self, u, v, h, du, dv, dh):
        for f, r_, d_ in ((u, self.ref[0], du), (v, self.ref[1], dv), (h, self.ref[2], dh)):
            zm = np.repeat(f.mean(axis=1, keepdims=True), f.shape[1], axis=1)
            d_ -= self.k * (zm - r_)
        return du, dv, dh
