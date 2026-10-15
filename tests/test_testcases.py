# -*- coding: utf-8 -*-
"""The standard test cases added in 0.5.0, and the Coriolis field they need.

The point of the rotated test case 2 is in `test_rotated_tc2_needs_its_own_coriolis`:
with the grid's own `f` the tilted flow is not a steady solution of anything, so a
model that cannot change `f` cannot run the test at all. Everything else here pins
the pieces that make that possible.
"""
from __future__ import annotations

import numpy as np
import pytest

from swesphere import SWEModel, presets
from swesphere.dynamics import rhs
from swesphere.initial_conditions import (rh_angular_phase_speed,
                                          rossby_haurwitz, tc2_coriolis,
                                          williamson_tc2)

DAY = 86400.0
LMAX = 16


def tendency_norm(model, x):
    """Size of dx/dt, scaled so the three variables are comparable."""
    du, dv, dh = rhs(*model.unpack(x), model.grid)
    k = model.n_pole_rows
    sl = slice(k, model.Nlat - k) if k else slice(None)
    return float(np.sqrt(np.mean(du[sl] ** 2 + dv[sl] ** 2)
                         + np.mean(dh[sl] ** 2) * (model.grid.g / 1e4)))


# ----------------------------------------------------------------------
# with_coriolis
# ----------------------------------------------------------------------
class TestWithCoriolis:
    def test_returns_a_copy_and_leaves_the_original_alone(self):
        g = SWEModel(LMAX=8).grid
        f0 = g.f_cor.copy()
        g2 = g.with_coriolis(np.zeros_like(g.f_cor))
        assert g2 is not g
        assert np.array_equal(g.f_cor, f0)
        assert np.all(g2.f_cor == 0.0)

    def test_geometry_is_untouched(self):
        g = SWEModel(LMAX=8).grid
        g2 = g.with_coriolis(np.zeros_like(g.f_cor))
        for name in ("LAT_ns", "LON_ns", "cosL", "sinL", "cos_safe", "kx"):
            assert np.array_equal(getattr(g, name), getattr(g2, name))
        assert (g2.LMAX, g2.Nlat, g2.Nlon, g2.a, g2.Omega, g2.g) == \
               (g.LMAX, g.Nlat, g.Nlon, g.a, g.Omega, g.g)

    def test_wrong_shape_is_rejected(self):
        g = SWEModel(LMAX=8).grid
        with pytest.raises(ValueError, match="shape"):
            g.with_coriolis(np.zeros((3, 3)))

    def test_model_set_coriolis_round_trips(self):
        m = SWEModel(LMAX=8)
        f0 = m.grid.f_cor.copy()
        assert m.set_coriolis(np.zeros_like(f0)) is m
        assert np.all(m.grid.f_cor == 0.0)
        m.set_coriolis(None)
        assert np.allclose(m.grid.f_cor, f0)


# ----------------------------------------------------------------------
# Test case 2
# ----------------------------------------------------------------------
class TestWilliamsonTC2:
    def test_alpha_zero_is_unchanged_from_0_4_0(self):
        """The default path must be bit identical, or every earlier result moves."""
        g = SWEModel(LMAX=LMAX).grid
        u, v, h = williamson_tc2(g, alpha_deg=0.0)
        U0, H0 = 38.0, 2800.0
        assert np.array_equal(u, U0 * g.cosL)
        assert np.array_equal(v, np.zeros_like(v))
        assert np.array_equal(
            h, H0 - (g.a * g.Omega * U0 + 0.5 * U0 ** 2) / g.g * g.sinL ** 2)
        assert np.array_equal(tc2_coriolis(g, 0.0), g.f_cor)

    def test_alpha_zero_never_touches_the_meridional_equation(self):
        """Why the angle exists: at alpha = 0, v is identically zero."""
        g = SWEModel(LMAX=LMAX).grid
        assert np.all(williamson_tc2(g, alpha_deg=0.0)[1] == 0.0)
        assert np.abs(williamson_tc2(g, alpha_deg=45.0)[1]).max() > 1.0

    @pytest.mark.parametrize("alpha", [0.0, 30.0, 45.0, 90.0])
    def test_wind_speed_is_solid_body_at_every_angle(self, alpha):
        g = SWEModel(LMAX=LMAX).grid
        u, v, _ = williamson_tc2(g, alpha_deg=alpha)
        # |U| = U0 cos(phi'), the latitude in the tilted frame.
        s = tc2_coriolis(g, alpha) / (2.0 * g.Omega)
        assert np.allclose(np.hypot(u, v), 38.0 * np.sqrt(1.0 - s ** 2),
                           atol=1e-9)

    @pytest.mark.parametrize("alpha", [30.0, 90.0])
    def test_rotated_tc2_needs_its_own_coriolis(self, alpha):
        """The whole reason `with_coriolis` exists.

        The test tilts the planet along with the flow, so it prescribes
        `f = 2 Omega sin(phi')`. With that field the state is steady. With the
        grid's own `f` the balance is not integrable, no height field holds the
        flow, and the tendency is large: what you would measure is the drift of
        a problem that was never in balance.
        """
        right = SWEModel(LMAX=LMAX, n_pole_rows=4)
        right.set_coriolis(tc2_coriolis(right.grid, alpha))
        x = right.pack(*williamson_tc2(right.grid, alpha_deg=alpha))

        wrong = SWEModel(LMAX=LMAX, n_pole_rows=4)
        assert tendency_norm(right, x) < 1e-2 * tendency_norm(wrong, x)

    def test_preset_installs_the_coriolis_field(self):
        m, _ = presets.williamson_tc2(alpha_deg=60.0, LMAX=8)
        assert np.allclose(m.grid.f_cor, tc2_coriolis(m.grid, 60.0))
        assert not np.allclose(m.grid.f_cor, 2.0 * m.grid.Omega * m.grid.sinL)

    def test_preset_at_alpha_zero_is_the_plain_earth(self):
        m, x0 = presets.williamson_tc2(LMAX=8)
        assert np.allclose(m.grid.f_cor, 2.0 * m.grid.Omega * m.grid.sinL)
        assert np.array_equal(x0, SWEModel(LMAX=8).initial_condition("tc2"))

    @pytest.mark.parametrize("alpha", [0.0, 90.0])
    def test_the_solution_stays_put(self, alpha):
        """Steady means steady: the discretization should barely move it.

        With no sponge, so this measures the scheme. How much the sponge adds
        on top is the interesting question and not one for a unit test: at
        LMAX 16 with 4 polar rows the same six hours moves the rotated case
        by 8 percent instead of 1.
        """
        m, x0 = presets.williamson_tc2(alpha_deg=alpha, LMAX=LMAX,
                                       n_pole_rows=0)
        h0 = m.unpack(x0)[2]
        h1 = m.unpack(m.propagate(x0, [0.0, 0.25 * DAY]))[2]
        rel = float(np.sqrt(np.mean((h1 - h0) ** 2))
                    / np.sqrt(np.mean(h0 ** 2)))
        assert rel < 2e-2


# ----------------------------------------------------------------------
# Test case 6
# ----------------------------------------------------------------------
class TestRossbyHaurwitz:
    def test_fields_are_finite_including_at_the_poles(self):
        g = SWEModel(LMAX=LMAX).grid
        for F in rossby_haurwitz(g):
            assert np.isfinite(F).all()

    def test_it_is_a_wavenumber_R_pattern(self):
        g = SWEModel(LMAX=LMAX).grid
        for R in (4, 5):
            h = rossby_haurwitz(g, R=R)[2]
            row = h[g.Nlat // 4] - h[g.Nlat // 4].mean()
            amp = np.abs(np.fft.rfft(row))
            assert amp.argmax() == R

    def test_it_is_not_zonal_and_not_a_tc2(self):
        g = SWEModel(LMAX=LMAX).grid
        u, v, h = rossby_haurwitz(g)
        assert np.abs(v).max() > 1.0
        assert (h - h.mean(axis=1, keepdims=True)).std() > 1.0

    def test_phase_speed_matches_the_published_value(self):
        g = SWEModel(LMAX=8).grid
        nu = rh_angular_phase_speed(g, R=4)
        assert nu > 0.0                                   # eastward
        assert 2.0 * np.pi / nu / DAY == pytest.approx(29.5, abs=0.5)

    def test_reachable_by_name(self):
        m = SWEModel(LMAX=8)
        g = m.grid
        for name in ("rossby_haurwitz", "tc6"):
            assert np.array_equal(m.initial_condition(name),
                                  m.pack(*rossby_haurwitz(g)))
        m2, x2 = presets.rossby_haurwitz(LMAX=8)
        assert np.array_equal(x2, m.initial_condition("rossby_haurwitz"))

    def test_it_holds_its_shape_for_a_day(self):
        m, x0 = presets.rossby_haurwitz(LMAX=LMAX, n_pole_rows=4)
        h0 = m.unpack(x0)[2]
        h1 = m.unpack(m.propagate(x0, [0.0, DAY]))[2]
        k = m.n_pole_rows
        a0 = np.abs(np.fft.rfft(h0[k:-k], axis=1))[:, 4].mean()
        a1 = np.abs(np.fft.rfft(h1[k:-k], axis=1))[:, 4].mean()
        assert 0.8 < a1 / a0 < 1.2
