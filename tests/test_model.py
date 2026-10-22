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
    x = m.propagate(x0, [0.0, 10 * DAY]); h = m.interior(m.unpack(x)[2])
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
    # A fixed band, not m.n_pole_rows. What this test must exclude is the two
    # rows at each pole where d_lat falls back to a one-sided first-order
    # difference, which is a property of the discretization and has nothing
    # to do with whether a sponge is switched on; tying it to n_pole_rows
    # made it select an empty array once the default became 0.
    k = 10; ref = -(h * divergence(u, v, G) + adv)
    assert np.abs(dh - ref)[k:-k].max() < 0.02 * np.abs(ref)[k:-k].max() and np.abs(adv)[k:-k].max() > 0.1 * np.abs(ref)[k:-k].max()


def test_tc2_is_steady_without_the_sponge():
    """Williamson test case 2: one day without the sponge keeps h within 1e-3 (l2) of the exact steady state."""
    from swesphere import diagnostics
    m = SWEModel(n_pole_rows=0); G = m.grid; x0 = m.initial_condition("tc2", U0=2 * np.pi * G.a / (12 * DAY), H0=2.94e4 / G.g)
    # Propagated with no sponge; n_pole_rows is then set only so that
    # error_norms(interior=True) scores away from the poles.
    x = m.propagate(x0, [0.0, DAY]); m.n_pole_rows = 10
    assert diagnostics.error_norms(m, x, x0, "h")["l2"] < 1e-3


def test_perturbed_initial_condition():
    m = presets.waves(); a = m.initial_condition("rossby", seed=1, perturbed=True); b = m.initial_condition("rossby", seed=2, perturbed=True)
    assert np.isfinite(a).all() and np.abs(a - b).max() > 0 and np.abs(a - m.initial_condition("rossby")).max() > 0


# ----------------------------------------------------------------------
# The contract of propagate.
#
# These exist because version 0.3.0 broke all four: the filter counter
# restarted on every call, so how a driver chopped up a forecast changed the
# model; a zero-length interval advanced a step; an interval that was not a
# whole number of steps landed somewhere else without saying so; and
# just_final_state=False raised NotImplementedError while the documentation
# said it returned the trajectory.
# ----------------------------------------------------------------------
def test_splitting_a_forecast_does_not_change_it():
    """The length of a driver's forecast window must not change the model."""
    m, x0 = presets.two_jets(LMAX=8)
    total = 12 * m.dt
    one = m.propagate(x0, [0.0, total])
    for per_call in (1, 2, 3, 4, 6):
        x = x0.copy()
        for k in range(12 // per_call):
            x = m.propagate(x, [k * per_call * m.dt, (k + 1) * per_call * m.dt])
        assert np.array_equal(x, one), (
            f"{per_call} step(s) per call disagrees with a single call"
        )


def test_the_filter_schedule_follows_absolute_time():
    """A call that starts mid-schedule keeps the schedule of the whole run."""
    m, x0 = presets.two_jets(LMAX=8)
    assert m.filter_every == 3
    # The step ending at t = 3*dt is a filter step; reaching it in a second
    # call must still filter there.
    one = m.propagate(x0, [0.0, 3 * m.dt])
    half = m.propagate(x0, [0.0, 2 * m.dt])
    two = m.propagate(half, [2 * m.dt, 3 * m.dt])
    assert np.array_equal(one, two)


def test_a_zero_length_interval_is_a_no_op():
    m, x0 = presets.two_jets(LMAX=8)
    assert np.array_equal(m.propagate(x0, [0.0, 0.0]), x0)
    assert np.array_equal(m.propagate(x0, [500.0, 500.0]), x0)


def test_a_partial_step_is_refused():
    m, x0 = presets.two_jets(LMAX=8)
    with pytest.raises(ValueError, match="whole number of time steps"):
        m.propagate(x0, [0.0, 0.5 * m.dt])
    with pytest.raises(ValueError, match="backwards"):
        m.propagate(x0, [m.dt, 0.0])


def test_just_final_state_false_returns_the_trajectory():
    m, x0 = presets.two_jets(LMAX=8)
    n = 5
    traj = m.propagate(x0, [0.0, n * m.dt], just_final_state=False)
    assert traj.shape == (n + 1, m.dim)
    assert np.array_equal(traj[0], x0)
    assert np.array_equal(traj[-1], m.propagate(x0, [0.0, n * m.dt]))


def test_members_are_independent_through_one_model_object():
    """An ensemble is propagated member by member on one model object.

    Nothing in propagate may carry over from one member to the next.
    """
    m, _ = presets.two_jets(LMAX=8)
    members = [m.initial_condition("rossby", seed=s, perturbed=True)
               for s in (1, 2, 3)]
    forward = [m.propagate(x, [0.0, 4 * m.dt]) for x in members]
    backward = [m.propagate(x, [0.0, 4 * m.dt]) for x in reversed(members)]
    for a, b in zip(forward, reversed(backward)):
        assert np.array_equal(a, b)


# ----------------------------------------------------------------------
# 0.6.0: the sponge is off by default
# ----------------------------------------------------------------------
def test_the_sponge_is_off_by_default():
    """Pinned, because this is the behaviour change of 0.6.0.

    EXP-10 of the model description paper ran 45 configurations for 120 days
    each, three presets by three resolutions by five sponge widths, and every
    one completed, including every sponge-free run. The sponge is not needed
    for stability; it costs nothing in time; it destroys the convergence of
    the scheme; and at LMAX 32 with the old default of ten rows it removed
    more than half of the midlatitude eddy amplitude.
    """
    m = SWEModel()
    assert m.n_pole_rows == 0
    assert np.all(m._sponge == 1.0)
    assert m.interior_mask().all()


def test_interior_returns_the_whole_field_with_no_sponge():
    """The trap the helper exists for: field[k:-k] is empty when k is 0."""
    m = SWEModel(LMAX=8)
    h = m.unpack(m.initial_condition("rossby"))[2]
    assert m.interior(h).shape == h.shape
    assert h[m.n_pole_rows:-m.n_pole_rows].size == 0      # the old idiom
    k = 4
    assert SWEModel(LMAX=8, n_pole_rows=k).interior(h).shape[0] == h.shape[0] - 2 * k


def test_the_sponge_still_works_when_asked_for():
    m = SWEModel(LMAX=8, n_pole_rows=3)
    assert m.n_pole_rows == 3
    assert m._sponge[0] < 1.0 and m._sponge[-1] < 1.0
    assert m._sponge[m.Nlat // 2] == 1.0
    assert m.interior_mask().sum() == (m.Nlat - 6) * m.Nlon
