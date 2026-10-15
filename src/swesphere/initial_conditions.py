# -*- coding: utf-8 -*-
"""
Initial conditions for the SWE on the sphere.

Two ingredients:

* **Solid-body rotation** (Williamson et al. 1992, Test Case 2). An
  exact stationary solution of the SWE on the sphere::

       u = U0 · cos φ
       v = 0
       h = H0 − (a Ω U0 + U0² / 2) / g · sin² φ

  The default base state of every IC.

* **Rossby-wave perturbations** in each hemisphere, gaussian-enveloped
  in latitude and zonally sinusoidal::

       δh = A_h · cos(k λ) · exp(−(φ − φ₀)² / (2 σ²))
       δu = A_uv · cos(k λ) · exp(...)
       δv = A_uv · sin(k λ) · exp(...)

  Default amplitudes and shapes match the reference setup from the
  upstream visualization script.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple

import numpy as np

from .grid import SphereGrid


@dataclass
class RossbyPert:
    """Rossby-wave-like perturbation envelope and amplitudes."""
    amp_h: float = 220.0      # Amplitude in h (m)
    amp_uv: float = 6.0       # Amplitude in (u, v) (m/s)
    wavenumber: int = 4       # Zonal wavenumber
    lat_center: float = 42.0  # Latitude of the gaussian center (deg)
    width: float = 12.0       # Gaussian width (deg)
    phase: float = 0.0        # Phase offset (rad)


def tc2_sin_rotated_latitude(grid: SphereGrid,
                             alpha_deg: float = 0.0) -> np.ndarray:
    """``sin`` of the latitude measured in the frame tilted by ``alpha_deg``."""
    a = np.radians(float(alpha_deg))
    return (-np.cos(grid.LON_ns) * grid.cosL * np.sin(a)
            + grid.sinL * np.cos(a))


def tc2_coriolis(grid: SphereGrid, alpha_deg: float = 0.0) -> np.ndarray:
    """The Coriolis field test case 2 prescribes at a rotation angle.

    ``f = 2 Omega sin(phi')``, the planetary rotation tilted along with
    everything else, which is what the test case specifies. **This is not
    optional.** Run the tilted flow against the grid's own
    ``f = 2 Omega sin(phi)`` and it is not a steady solution of anything: the
    balance

        dB/dlam = a cos(phi) (f + zeta) v,   dB/dphi = -a (f + zeta) u

    then fails its own integrability condition by
    ``-2 Omega U0 a sin(alpha) sin(lam) cos^2(phi)``, so no height field at
    all holds it steady, and the drift you would measure would be the drift
    of a problem that was never in balance.

    Pass it through :meth:`SphereGrid.with_coriolis`, or let
    :func:`swesphere.presets.williamson_tc2` wire it up for you.
    """
    return 2.0 * grid.Omega * tc2_sin_rotated_latitude(grid, alpha_deg)


def williamson_tc2(grid: SphereGrid, U0: float = 38.0,
                   H0: float = 2800.0, alpha_deg: float = 0.0
                   ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Williamson 1992 Test Case 2: solid-body rotation at an angle.

    ``alpha_deg`` is the angle between the axis of the flow and the grid's
    polar axis, and the default of 0 reproduces the zonally symmetric case
    bit for bit. At 90 degrees the flow runs pole to pole.

    The angle matters for verification. At ``alpha_deg = 0`` the solution has
    ``v == 0`` everywhere, so the meridional momentum equation is never
    evaluated on anything but zero, and the meridional derivatives are all
    taken along grid lines. A scheme can look second-order accurate on that
    case and be considerably worse on any flow that crosses a pole.

    Exact and steady for every angle, **provided the model uses**
    :func:`tc2_coriolis`.
    """
    a = np.radians(float(alpha_deg))
    u = U0 * (grid.cosL * np.cos(a)
              + grid.sinL * np.cos(grid.LON_ns) * np.sin(a))
    v = -U0 * np.sin(grid.LON_ns) * np.sin(a)
    rot = tc2_sin_rotated_latitude(grid, alpha_deg)
    h = H0 - (grid.a * grid.Omega * U0 + 0.5 * U0 ** 2) / grid.g * rot ** 2
    return u, v, h


def rossby_haurwitz(grid: SphereGrid, R: int = 4, omega: float = 7.848e-6,
                    K: float = 7.848e-6, h0: float = 8000.0
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Williamson 1992 Test Case 6: a Rossby-Haurwitz wave.

    A wavenumber-``R`` pattern on a zonal flow, neither zonal nor steady. In
    the non-divergent barotropic limit it travels east without changing
    shape, at the rate :func:`rh_angular_phase_speed` gives; in the
    shallow-water system it holds its form for weeks. Unlike test case 2 it
    checks that the model carries a structure at the right speed, not only
    that it leaves one alone.
    """
    lam = grid.LON_ns
    c, s = grid.cosL, grid.sinL
    a, Om, g = grid.a, grid.Omega, grid.g

    u = a * omega * c + a * K * c ** (R - 1) * (R * s ** 2 - c ** 2) * np.cos(R * lam)
    v = -a * K * R * c ** (R - 1) * s * np.sin(R * lam)

    # The cos^(-2) term of A is folded into the power, so nothing divides by
    # cos(phi) at the poles.
    A = (0.5 * omega * (2.0 * Om + omega) * c ** 2
         + 0.25 * K ** 2 * c ** (2 * R) * ((R + 1) * c ** 2 + (2 * R ** 2 - R - 2))
         - 0.5 * K ** 2 * R ** 2 * c ** (2 * R - 2))
    B = (2.0 * (Om + omega) * K / ((R + 1) * (R + 2)) * c ** R
         * ((R ** 2 + 2 * R + 2) - (R + 1) ** 2 * c ** 2))
    C = 0.25 * K ** 2 * c ** (2 * R) * ((R + 1) * c ** 2 - (R + 2))

    h = h0 + (a ** 2 / g) * (A + B * np.cos(R * lam) + C * np.cos(2 * R * lam))
    return u, v, h


def rh_angular_phase_speed(grid: SphereGrid, R: int = 4,
                           omega: float = 7.848e-6) -> float:
    """Eastward angular phase speed of the Rossby-Haurwitz wave, rad/s.

        nu = (R (3 + R) omega - 2 Omega) / ((1 + R) (2 + R))

    About 2.46e-6 for the standard R = 4, a pattern that circles the globe in
    roughly 29.5 days. Compare a measured phase against this.
    """
    return (R * (3.0 + R) * omega - 2.0 * grid.Omega) / ((1.0 + R) * (2.0 + R))


def add_rossby_perturbation(
    u: np.ndarray, v: np.ndarray, h: np.ndarray,
    grid: SphereGrid, pert: RossbyPert,
):
    """Add a single Rossby-wave perturbation in-place."""
    phi0 = np.radians(pert.lat_center)
    sigma = np.radians(pert.width)
    env = np.exp(-((grid.LAT_ns - phi0) ** 2) / (2 * sigma ** 2))
    cos_arg = pert.wavenumber * grid.LON_ns + pert.phase
    sin_arg = pert.wavenumber * grid.LON_ns + pert.phase
    h += pert.amp_h * np.cos(cos_arg) * env
    u += pert.amp_uv * np.cos(cos_arg) * env
    v += pert.amp_uv * np.sin(sin_arg) * env


def default_ic(grid: SphereGrid, U0: float = 38.0,
               H0: float = 2800.0
               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Default IC: TC2 + NH wavenumber-4 + SH wavenumber-3 perturbations."""
    u, v, h = williamson_tc2(grid, U0=U0, H0=H0)
    add_rossby_perturbation(u, v, h, grid, RossbyPert(
        amp_h=220.0, amp_uv=6.0, wavenumber=4,
        lat_center=42.0, width=12.0, phase=0.0,
    ))
    add_rossby_perturbation(u, v, h, grid, RossbyPert(
        amp_h=150.0, amp_uv=0.0, wavenumber=3,
        lat_center=-38.0, width=11.0, phase=1.1,
    ))
    # Small extra meridional component on the SH band, matching the
    # reference setup (proportional to amp_h / 50 in the original).
    phi2 = np.radians(-38.0)
    env2 = np.exp(-((grid.LAT_ns - phi2) ** 2) / (2 * np.radians(11.0) ** 2))
    v += (150.0 / 50) * np.sin(3 * grid.LON_ns + 1.1) * env2
    return u, v, h


def perturbed_ic(
    grid: SphereGrid,
    seed: int = 0,
    U0: float = 38.0,
    H0: float = 2800.0,
    pert_strength: float = 0.05,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Default IC plus a small random perturbation of ``u`` and ``h`` at zonal wavenumbers 1-8,
    confined to low and middle latitudes, to start runs from different states.

    The perturbations of ``u`` and ``h`` have independent phases, so they are not in geostrophic
    balance: they excite weak gravity waves that the model radiates and filters within a few days.
    Use a spin-up before sampling.
    """
    u, v, h = default_ic(grid, U0=U0, H0=H0)
    rng = np.random.default_rng(seed)

    # Low zonal wavenumbers, gaussian envelope in latitude.
    Nlon = grid.Nlon
    n_modes = 8
    for k in range(1, n_modes + 1):
        amp_u = pert_strength * U0 * rng.standard_normal()
        amp_h = pert_strength * H0 * 0.05 * rng.standard_normal()
        phase_u = 2 * np.pi * rng.random()
        phase_h = 2 * np.pi * rng.random()
        # Latitude-banded modulation
        lat_mod = np.exp(-(grid.LAT_ns / np.radians(45.0)) ** 2)
        u += amp_u * np.cos(k * grid.LON_ns + phase_u) * lat_mod
        h += amp_h * np.cos(k * grid.LON_ns + phase_h) * lat_mod
    return u, v, h


# ---------------------------------------------------------------------------
# Galewsky, Scott & Polvani (2004): a balanced, barotropically unstable jet
# ---------------------------------------------------------------------------
def galewsky_jet(grid: SphereGrid, umax: float = 80.0, phi0_deg: float = 25.7, phi1_deg: float = 64.6,
                 h0: float = 10000.0, bump: bool = True, both_hemispheres: bool = True):
    """Balanced zonal jet(s) with the small height bump that triggers the instability.

    Returns ``(u, v, h)`` on the grid. With ``both_hemispheres`` a mirror jet
    (also westerly) is placed in the south so the whole sphere is active; the
    two hemispheres perturb each other across the equator and the resulting
    turbulence is irregular (with the zonal relaxation of :mod:`swesphere.forcing`).
    """
    Nlat, Nlon = grid.Nlat, grid.Nlon
    lats = np.radians(90 - np.arange(Nlat) * 180 / Nlat); lons = np.radians(np.arange(Nlon) * 360 / Nlon)
    phi0, phi1 = np.radians(phi0_deg), np.radians(phi1_deg); en = np.exp(-4.0 / (phi1 - phi0) ** 2)

    def uprof(phi):
        out = np.zeros_like(phi); m = (phi > phi0) & (phi < phi1)
        out[m] = umax / en * np.exp(1.0 / ((phi[m] - phi0) * (phi[m] - phi1)))
        if both_hemispheres:
            ms = (-phi > phi0) & (-phi < phi1)
            out[ms] = umax / en * np.exp(1.0 / ((-phi[ms] - phi0) * (-phi[ms] - phi1)))
        return out

    fine = np.linspace(-np.pi / 2, np.pi / 2, 4001); uf = uprof(fine); f = 2 * grid.Omega * np.sin(fine)
    dh = -grid.a * uf * (f + uf * np.tan(fine) / grid.a) / grid.g
    hf = np.concatenate([[0.0], np.cumsum(0.5 * (dh[1:] + dh[:-1]) * np.diff(fine))]); hf = h0 + hf - hf.mean()
    u = np.repeat(uprof(lats)[:, None], Nlon, axis=1); v = np.zeros_like(u)
    h = np.repeat(np.interp(lats, fine, hf)[:, None], Nlon, axis=1)
    if bump:
        LON, LAT = np.meshgrid(lons, lats); phi2 = np.pi / 4; alpha, beta = 1 / 3, 1 / 15
        h = h + 120.0 * np.cos(LAT) * np.exp(-((LON - np.pi) / alpha) ** 2) * np.exp(-((phi2 - LAT) / beta) ** 2)
        if both_hemispheres:
            h = h + 120.0 * np.cos(LAT) * np.exp(-((LON - 0.4 * np.pi) / alpha) ** 2) * np.exp(-((-phi2 - LAT) / beta) ** 2)
    return u, v, h
