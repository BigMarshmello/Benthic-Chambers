"""Stage 0 tests: the closed-form relations later stages must reproduce."""

import math

import pytest

from benthic import scaling as sc

PHI = 0.40
D0 = 1.8e-9
C0 = 0.25          # mol m-3, air saturation
RATE = 1.0e-4      # mol m-3(porewater) s-1
MU = 1.19e-3
RHO = 1025.0
G = 9.81


# ----------------------------------------------------------------- tortuosity


def test_boudreau_tortuosity_known_values():
    # theta^2 = 1 - ln(phi^2); at phi = 1 the sediment is open water.
    assert math.isclose(sc.tortuosity_squared(1.0), 1.0)
    assert math.isclose(sc.tortuosity_squared(0.40), 1.0 - math.log(0.16))
    assert sc.tortuosity_squared(0.40) == pytest.approx(2.8326, abs=1e-4)


def test_tortuosity_increases_as_porosity_falls():
    assert sc.tortuosity_squared(0.35) > sc.tortuosity_squared(0.45)


def test_sediment_diffusivity_reduces_d0():
    d_s = sc.sediment_diffusivity(D0, PHI)
    assert d_s < D0
    assert d_s == pytest.approx(D0 / 2.8326, rel=1e-4)


def test_tortuosity_rejects_nonphysical_porosity():
    with pytest.raises(ValueError):
        sc.tortuosity_squared(0.0)
    with pytest.raises(ValueError):
        sc.tortuosity_squared(1.5)


# -------------------------------------------- 1D diffusion + zero-order sink


def test_penetration_depth_matches_hand_calculation():
    d_s = sc.sediment_diffusivity(D0, PHI)
    L = sc.penetration_depth(d_s, C0, RATE)
    assert L == pytest.approx(math.sqrt(2 * d_s * C0 / RATE))
    assert L * 1e3 == pytest.approx(1.783, abs=1e-3)   # mm


def test_the_two_basis_conventions_agree_when_converted():
    """R_bulk = phi*R_porewater must give the same penetration depth.

    This is the discrepancy in the project brief: the brief's PDE uses
    phi*R (porewater basis) but its analytic L carries a phi inside the
    root (bulk basis). Only one can be right for a given R.
    """
    d_s = sc.sediment_diffusivity(D0, PHI)
    L_pw = sc.penetration_depth(d_s, C0, RATE)
    L_bulk = sc.penetration_depth_bulk(PHI, d_s, C0, PHI * RATE)
    assert L_pw == pytest.approx(L_bulk, rel=1e-12)


def test_naive_use_of_the_briefs_formula_is_wrong_by_sqrt_phi():
    """Guard against silently mixing the conventions."""
    d_s = sc.sediment_diffusivity(D0, PHI)
    correct = sc.penetration_depth(d_s, C0, RATE)
    naive = sc.penetration_depth_bulk(PHI, d_s, C0, RATE)  # same R, wrong basis
    assert naive == pytest.approx(correct * math.sqrt(PHI), rel=1e-12)
    assert naive < correct


def test_flux_is_consistent_with_the_profile_gradient():
    """J must equal phi*D_s*dC/dz at z=0 for C = C_0*(1 - z/L)^2."""
    d_s = sc.sediment_diffusivity(D0, PHI)
    L = sc.penetration_depth(d_s, C0, RATE)
    J = sc.diffusive_flux(PHI, d_s, C0, RATE)
    assert J == pytest.approx(2.0 * PHI * d_s * C0 / L, rel=1e-12)


def test_flux_equals_depth_integrated_consumption():
    """Steady state: what goes in must be consumed. Closes the mass balance."""
    d_s = sc.sediment_diffusivity(D0, PHI)
    L = sc.penetration_depth(d_s, C0, RATE)
    J = sc.diffusive_flux(PHI, d_s, C0, RATE)
    assert J == pytest.approx(PHI * RATE * L, rel=1e-12)


def test_damkohler_is_two_at_the_penetration_depth():
    d_s = sc.sediment_diffusivity(D0, PHI)
    L = sc.penetration_depth(d_s, C0, RATE)
    assert sc.damkohler(RATE, L, d_s, C0) == pytest.approx(2.0, rel=1e-12)


# ------------------------------------------------ Elliott & Brooks bed head


def test_elliott_brooks_exponent_switches_at_ratio_0p34():
    """3/8 below H/d = 0.34, 3/2 above; continuous at the switch."""
    d = 1.0
    below = sc.elliott_brooks_head_amplitude(0.2, 0.3399, d)
    above = sc.elliott_brooks_head_amplitude(0.2, 0.3401, d)
    assert below == pytest.approx(above, rel=1e-3)
    # at exactly H/d = 0.34 the bracket is 1, so both branches give 0.28*U^2/2g
    at = sc.elliott_brooks_head_amplitude(0.2, 0.34, d)
    assert at == pytest.approx(0.28 * 0.2**2 / (2 * 9.81), rel=1e-12)


def test_elliott_brooks_scales_with_velocity_squared():
    a = sc.elliott_brooks_head_amplitude(0.10, 0.03, 0.5)
    b = sc.elliott_brooks_head_amplitude(0.20, 0.03, 0.5)
    assert b == pytest.approx(4.0 * a, rel=1e-12)


def test_elliott_brooks_fails_at_zero_depth():
    # Intertidal site: the relation is undefined as the flat drains.
    with pytest.raises(ValueError):
        sc.elliott_brooks_head_amplitude(0.10, 0.03, 0.0)


def test_elliott_brooks_reproduces_huettel_gust_ripple_gradient():
    """Independent cross-check of two separate published relations.

    Huettel & Gust (1992) state that a 10 cm/s bottom flow over a 7 cm
    ripple produces ~0.2 Pa/cm. Elliott & Brooks, given the same flow and
    bedform and a plausible ripple wavelength, must land in the same place.
    """
    h_m = sc.elliott_brooks_head_amplitude(0.10, 0.07, 0.5, G)
    p_m = sc.head_to_pressure(h_m, RHO, G)
    for wavelength in (0.35, 0.50):
        grad_pa_per_cm = sc.sinusoid_max_gradient(p_m, wavelength) / 100.0
        assert 0.10 < grad_pa_per_cm < 0.30, grad_pa_per_cm


# --------------------------------------------------------- stirrer forcing


def test_stirrer_pressure_round_trips_with_rim_gradient():
    radius = 0.15
    dps = sc.stirrer_pressure_from_rim_gradient(20.0, radius)  # 0.2 Pa/cm
    assert dps == pytest.approx(1.5)
    assert sc.stirrer_rim_gradient(dps, radius) == pytest.approx(20.0)


def test_stirrer_forcing_is_comparable_to_field_forcing():
    """The chamber default must not be wildly stronger than the bed it mimics."""
    h_m = sc.elliott_brooks_head_amplitude(0.10, 0.03, 0.5, G)
    field_grad = sc.sinusoid_max_gradient(sc.head_to_pressure(h_m, RHO, G), 0.30)
    chamber_grad = sc.stirrer_rim_gradient(1.5, 0.15)
    assert 0.3 < chamber_grad / field_grad < 3.0


# ------------------------------------------- advection/diffusion crossover


def test_darcy_velocity_is_linear_in_permeability():
    q1 = sc.darcy_velocity(1e-11, 15.7, MU)
    q2 = sc.darcy_velocity(2e-11, 15.7, MU)
    assert q2 == pytest.approx(2.0 * q1)


def test_crossover_permeability_gives_unit_peclet():
    d_s = sc.sediment_diffusivity(D0, PHI)
    L = sc.penetration_depth(d_s, C0, RATE)
    grad = 15.68
    k = sc.crossover_permeability(MU, PHI, d_s, grad, L)
    q = sc.darcy_velocity(k, grad, MU)
    assert sc.peclet(q, L, PHI, d_s) == pytest.approx(1.0, rel=1e-12)


def test_users_quoted_thresholds_bracket_the_crossover():
    """The brief's 3e-12 and 2.6e-11 m2 thresholds are independently plausible.

    At the default bed forcing they sit either side of Pe = 1, i.e. they
    mark 'advection becomes non-negligible' and 'advection dominates'.
    """
    d_s = sc.sediment_diffusivity(D0, PHI)
    L = sc.penetration_depth(d_s, C0, RATE)
    grad = 15.68
    pe_lo = sc.peclet(sc.darcy_velocity(3.0e-12, grad, MU), L, PHI, d_s)
    pe_hi = sc.peclet(sc.darcy_velocity(2.6e-11, grad, MU), L, PHI, d_s)
    assert pe_lo < 1.0 < pe_hi
    assert 0.1 < pe_lo < 0.5
    assert 1.5 < pe_hi < 4.0


def test_irrigation_number_is_small_across_the_oxic_layer():
    """Irrigation is a parallel supply path, not a reshaping of the microprofile."""
    d_s = sc.sediment_diffusivity(D0, PHI)
    L = sc.penetration_depth(d_s, C0, RATE)
    assert sc.irrigation_number(5e-6, L, d_s) < 0.1
    # but over the burrow depth scale it is enormous
    assert sc.irrigation_number(5e-6, 0.30, d_s) > 100.0
