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
