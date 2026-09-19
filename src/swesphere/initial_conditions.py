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


def williamson_tc2(grid: SphereGrid, U0: float = 38.0,
                   H0: float = 2800.0
                   ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Williamson 1992 Test Case 2: solid-body rotation."""
    u = U0 * grid.cosL
    v = np.zeros_like(u)
    h = H0 - (grid.a * grid.Omega * U0 + 0.5 * U0 ** 2) / grid.g * grid.sinL ** 2
    return u, v, h


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
