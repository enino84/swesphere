# -*- coding: utf-8 -*-
"""Time stepper for the shallow-water right-hand side: fourth-order Runge-Kutta.

RK4 with the spectral filter and the polar sponge keeps every preset stationary
for hundreds of days at dt = 120 s (see EXP-01 in ``experiments/`` for the
stable range of dt and filter cadence).
"""
from __future__ import annotations


def rk4(rhs, u, v, h, dt):
    k1 = rhs(u, v, h)
    k2 = rhs(u + 0.5 * dt * k1[0], v + 0.5 * dt * k1[1], h + 0.5 * dt * k1[2])
    k3 = rhs(u + 0.5 * dt * k2[0], v + 0.5 * dt * k2[1], h + 0.5 * dt * k2[2])
    k4 = rhs(u + dt * k3[0], v + dt * k3[1], h + dt * k3[2])
    return (u + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
            v + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]),
            h + dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]))


INTEGRATORS = {"rk4": rk4}
