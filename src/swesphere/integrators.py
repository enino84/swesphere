# -*- coding: utf-8 -*-
"""Time steppers for the shallow-water right-hand side.

``euler`` is the original scheme of the pyteda implementation. It is
unconditionally unstable for the wave part of the system: the spectral
filter absorbs the growth for a few weeks and the run blows up (measured:
day 25 at dt = 120 s, day 178 at dt = 60 s). ``rk4`` with the same
right-hand side, filter and sponge is stable (100+ days at dt = 120 s with
a stationary energy). Use ``rk4``.
"""
from __future__ import annotations


def euler(rhs, u, v, h, dt):
    du, dv, dh = rhs(u, v, h)
    return u + dt * du, v + dt * dv, h + dt * dh


def rk4(rhs, u, v, h, dt):
    k1 = rhs(u, v, h)
    k2 = rhs(u + 0.5 * dt * k1[0], v + 0.5 * dt * k1[1], h + 0.5 * dt * k1[2])
    k3 = rhs(u + 0.5 * dt * k2[0], v + 0.5 * dt * k2[1], h + 0.5 * dt * k2[2])
    k4 = rhs(u + dt * k3[0], v + dt * k3[1], h + dt * k3[2])
    return (u + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
            v + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]),
            h + dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]))


INTEGRATORS = {"euler": euler, "rk4": rk4}
