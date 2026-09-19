import numpy as np, pytest
from swesphere import presets, SWEModel
DAY = 86400.0


def test_layout_and_masks():
    m = presets.waves(); x = m.initial_condition("rossby")
    assert x.size == m.dim and m.interior_mask().sum() == (m.Nlat - 2 * m.n_pole_rows) * m.Nlon
    u, v, h = m.unpack(x); assert np.allclose(m.pack(u, v, h), x)


def test_rk4_keeps_energy_over_two_days():
    """Two days of the waves preset: the h std changes by less than 5%."""
    m = presets.waves(); x = m.initial_condition("rossby")
    s0 = m.unpack(x)[2].std(); s2 = m.unpack(m.propagate(x, [0.0, 2 * DAY]))[2].std()
    assert abs(s2 - s0) / s0 < 0.05


def test_unknown_scheme_is_rejected():
    with pytest.raises(KeyError):
        SWEModel(scheme="not-a-scheme")


def test_two_jets_sustains_eddies():
    """Ten days are enough to see the instability grow from the 4 m bump toward tens of metres."""
    m, x0 = presets.two_jets()
    x = m.propagate(x0, [0.0, 10 * DAY]); h = m.unpack(x)[2][m.n_pole_rows:-m.n_pole_rows]
    eddy = (h - h.mean(axis=1, keepdims=True)).std()
    assert 20.0 < eddy < 400.0 and np.isfinite(x).all()


def test_galewsky_two_jets():
    """Both hemispheres carry a westerly jet of the same strength (the grid rows are not
    exactly mirror-symmetric about the equator, so compare the maxima, not the rows)."""
    m = presets.waves(); u, v, h = m.unpack(m.initial_condition("galewsky", bump=False))
    half = m.Nlat // 2
    assert abs(u[:half].max() - u[half:].max()) < 2.0 and u.max() > 70.0 and np.all(u >= -1e-9)


def test_mass_is_conserved():
    """Flux-form continuity: two unforced days change the global mean depth by less than 1e-4
    (what remains is the filter and the polar rows; dh = -h*div without the advection of h loses ~1e-3)."""
    from swesphere import diagnostics
    m = presets.waves(); x = m.initial_condition("rossby")
    m0 = diagnostics.invariants(m, x)["mass"]; m2 = diagnostics.invariants(m, m.propagate(x, [0.0, 2 * DAY]))["mass"]
    assert abs(m2 - m0) / m0 < 1e-4


def test_continuity_is_the_divergence_of_the_mass_flux():
    """dh/dt = -h*div - u.grad(h), checked pointwise away from the poles on a smooth state."""
    from swesphere.dynamics import rhs, divergence
    from swesphere.grid import d_lon, d_lat
    m = presets.waves(); G = m.grid; u, v, h = m.unpack(m.initial_condition("rossby")); v = v + 5.0 * np.cos(2 * G.LON_ns) * G.cosL
    dh = rhs(u, v, h, G)[2]; adv = u * d_lon(h, G) / (G.a * G.cos_safe) + v * d_lat(h, G) / G.a
    k = m.n_pole_rows; ref = -(h * divergence(u, v, G) + adv)
    assert np.abs(dh - ref)[k:-k].max() < 0.02 * np.abs(ref)[k:-k].max() and np.abs(adv)[k:-k].max() > 0.1 * np.abs(ref)[k:-k].max()


def test_tc2_is_steady_without_the_sponge():
    """Williamson test case 2: one day without the sponge keeps h within 1e-3 (l2) of the exact steady state."""
    from swesphere import diagnostics
    m = SWEModel(n_pole_rows=0); G = m.grid; x0 = m.initial_condition("tc2", U0=2 * np.pi * G.a / (12 * DAY), H0=2.94e4 / G.g)
    x = m.propagate(x0, [0.0, DAY]); m.n_pole_rows = 10
    assert diagnostics.error_norms(m, x, x0, "h")["l2"] < 1e-3


def test_perturbed_initial_condition():
    m = presets.waves(); a = m.initial_condition("rossby", seed=1, perturbed=True); b = m.initial_condition("rossby", seed=2, perturbed=True)
    assert np.isfinite(a).all() and np.abs(a - b).max() > 0 and np.abs(a - m.initial_condition("rossby")).max() > 0
