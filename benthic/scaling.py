"""Closed-form scaling relations used to sanity-check every later stage.

Nothing here solves a PDE.  These are the back-of-envelope results that
the finite-volume model must reproduce; if stage 2 or stage 3 disagrees
with these by more than a few percent in the limits where they apply,
the solver is wrong, not the algebra.

Conventions (see config/config.yaml header):
  z positive downward, z = 0 at the sediment-water interface
  fluxes positive INTO the sediment
  R and alpha are per unit volume of POREWATER
"""

from __future__ import annotations

import math

# --------------------------------------------------------------------
# Sediment diffusion
# --------------------------------------------------------------------


def tortuosity_squared(phi: float) -> float:
    """Boudreau (1996) modified Weissberg relation, theta^2 = 1 - ln(phi^2)."""
    # phi = 1 is the open-water limit and must return theta^2 = 1.
    if not 0.0 < phi <= 1.0:
        raise ValueError("porosity must be in (0, 1]")
    return 1.0 - math.log(phi * phi)


def sediment_diffusivity(d0: float, phi: float) -> float:
    """Effective porewater diffusivity D_s = D_0 / theta^2."""
    return d0 / tortuosity_squared(phi)


# --------------------------------------------------------------------
# 1D diffusion + zero-order consumption (the stage-1 analytic target)
# --------------------------------------------------------------------
#
# Steady state, no advection, no irrigation, R per unit porewater volume:
#
#     d/dz( phi * D_s * dC/dz ) = phi * R        =>   C'' = R / D_s
#
# with C(0) = C_0, C(L) = 0, C'(L) = 0 gives
#
#     C(z) = C_0 * (1 - z/L)^2,      L = sqrt(2 * D_s * C_0 / R)
#     J(0) = phi * sqrt(2 * D_s * C_0 * R)
#
# NOTE the absence of phi inside the square root for L.  The prompt's
# L = sqrt(2*phi*D_s*C_0/R) is the BULK-basis result (reaction term R,
# not phi*R).  Both are below; they differ by sqrt(phi) ~ 0.63.


def penetration_depth(d_s: float, c0: float, rate: float, basis: str = "porewater") -> float:
    """Oxygen penetration depth for zero-order consumption.

    basis='porewater': rate is per m3 porewater, term is phi*R  -> no phi
    basis='bulk':      rate is per m3 bulk sediment, term is R  -> phi appears
    """
    if basis == "porewater":
        return math.sqrt(2.0 * d_s * c0 / rate)
    if basis == "bulk":
        raise ValueError("bulk basis needs porosity; use penetration_depth_bulk")
    raise ValueError(f"unknown basis {basis!r}")


def penetration_depth_bulk(phi: float, d_s: float, c0: float, rate_bulk: float) -> float:
    """L = sqrt(2*phi*D_s*C_0/R_bulk) -- the form written in the project brief."""
    return math.sqrt(2.0 * phi * d_s * c0 / rate_bulk)


def diffusive_flux(phi: float, d_s: float, c0: float, rate: float) -> float:
    """Steady diffusive O2 uptake, mol m-2 s-1, positive into the bed."""
    return phi * math.sqrt(2.0 * d_s * c0 * rate)


# --------------------------------------------------------------------
# Bed pressure forcing
# --------------------------------------------------------------------


def elliott_brooks_head_amplitude(u: float, height: float, depth: float, g: float = 9.81) -> float:
    """Half-amplitude of the dynamic head over a bedform, in metres of water.

    h_m = 0.28 * (U^2 / 2g) * (H / (0.34*d))^n,  n = 3/8 if H/d <= 0.34 else 3/2

    Elliott & Brooks (1997), WRR 33:123-136, fitted to Fehlman (1985).
    U is the MEAN stream velocity of the flume calibration, H the bedform
    height, d the water depth.  Undefined as d -> 0.
    """
    if depth <= 0.0:
        raise ValueError("water depth must be positive (relation fails at d -> 0)")
    ratio = height / depth
    exponent = 3.0 / 8.0 if ratio <= 0.34 else 3.0 / 2.0
    return 0.28 * (u * u / (2.0 * g)) * (height / (0.34 * depth)) ** exponent


def head_to_pressure(h_m: float, rho: float = 1025.0, g: float = 9.81) -> float:
    """Convert a head amplitude (m) to a dynamic pressure amplitude (Pa)."""
    return rho * g * h_m


def sinusoid_max_gradient(p_amplitude: float, wavelength: float) -> float:
    """Peak |dp/dx| for p = p_m*sin(2*pi*x/lambda), in Pa m-1."""
    return p_amplitude * 2.0 * math.pi / wavelength


def stirrer_pressure_from_rim_gradient(gradient: float, radius: float) -> float:
    """delta_p_s giving a prescribed |dp/dr| at the chamber wall.

    For p(r) = delta_p_s * (r/R)^2 the rim gradient is 2*delta_p_s/R.
    Huettel & Gust (1992) set 0.2 Pa cm-1 = 20 Pa m-1.
    """
    return gradient * radius / 2.0


def stirrer_rim_gradient(delta_p_s: float, radius: float) -> float:
    """Inverse of :func:`stirrer_pressure_from_rim_gradient`."""
    return 2.0 * delta_p_s / radius


# --------------------------------------------------------------------
# Darcy velocity and the advection/diffusion competition
# --------------------------------------------------------------------


def darcy_velocity(k: float, grad_p: float, mu: float) -> float:
    """Darcy flux magnitude |q| = k*|grad p|/mu, m s-1 (volume per bulk area)."""
    return k * grad_p / mu


def peclet(q: float, length: float, phi: float, d_s: float) -> float:
    """Advective / diffusive flux ratio over ``length``.

    Advective O2 flux ~ q*C; diffusive ~ phi*D_s*C/length.  So
    Pe = q*length/(phi*D_s).  Using the oxygen penetration depth as the
    length scale answers the question that actually matters here: does
    porewater flow move more O2 across the oxic layer than diffusion does?
    """
    return q * length / (phi * d_s)


def damkohler(rate: float, length: float, d_s: float, c0: float) -> float:
    """Reaction vs diffusive transport over ``length``: Da = R*length^2/(D_s*C_0).

    Da = 2 exactly when length is the analytic penetration depth.
    """
    return rate * length * length / (d_s * c0)


def irrigation_number(alpha: float, length: float, d_s: float) -> float:
    """Irrigative vs diffusive supply: alpha*length^2/D_s."""
    return alpha * length * length / d_s


def crossover_permeability(
    mu: float, phi: float, d_s: float, grad_p: float, length: float
) -> float:
    """Permeability at which Pe = 1 over ``length``.

    Since q is linear in k, k_cross = mu*phi*D_s/(|grad p|*length).
    Below it the oxic layer is diffusion-controlled, above it
    advection-controlled.
    """
    return mu * phi * d_s / (grad_p * length)


# --------------------------------------------------------------------
# Chamber design: signal size and replicate variance
# --------------------------------------------------------------------
#
# Added at revision 1, once the real production geometry was known
# (168 mm bore, ~600 mm tube).  Both results below are properties of the
# chamber, not of the sediment, so they can be settled before any PDE is
# solved -- and they turned out to be the binding design constraints.


def chamber_area(bore: float) -> float:
    """Enclosed sediment area, m2. ``bore`` is the INTERNAL diameter."""
    return math.pi * (bore / 2.0) ** 2


def chamber_volume(bore: float, height: float) -> float:
    """Enclosed water volume, m3, for a plain cylinder."""
    return chamber_area(bore) * height


def drawdown_rate(flux: float, height: float) -> float:
    """dC_w/dt in mol m-3 s-1 for a uniform flux over the enclosed area.

    V dC_w/dt = -J*A and V = A*h, so the area cancels exactly:
    dC_w/dt = J/h.  The signal depends on the WATER COLUMN HEIGHT alone,
    not on the chamber diameter.  Diameter still matters, but for
    replicate variance (see prawn_count_cv), not for signal size.
    """
    if height <= 0.0:
        raise ValueError("water column height must be positive")
    return flux / height


def time_to_drawdown(flux: float, height: float, c0: float, fraction: float) -> float:
    """Seconds for C_w to fall by ``fraction`` of c0, at constant flux.

    Optimistic: it ignores the decline in flux as C_w falls, so the real
    time is longer.
    """
    return fraction * c0 / drawdown_rate(flux, height)


def aspect_ratio(height: float, bore: float) -> float:
    """Water column height divided by bore.

    Published stirred chambers run near 1.5-1.7; far above that, the
    well-mixed assumption and the swirl reaching the bed both weaken.
    """
    return height / bore


def signal_to_noise(
    flux: float, height: float, duration: float, resolution: float, drift: float = 0.0
) -> float:
    """Concentration change over a deployment, in units of sensor noise.

    ``resolution`` and ``drift`` are in mol m-3 and mol m-3 s-1.
    Below ~3 the flux is not measurable against the instrument.
    """
    signal = drawdown_rate(flux, height) * duration
    noise = resolution + drift * duration
    if noise <= 0.0:
        return float("inf")
    return signal / noise


def prawn_count_mean(density: float, area: float) -> float:
    """Expected number of prawns enclosed by one chamber."""
    return density * area


def prawn_zero_probability(density: float, area: float) -> float:
    """P(a chamber encloses no prawns), assuming Poisson-distributed burrows.

    Poisson is the right null model for randomly placed burrows. Real
    callianassid beds are patchy, which makes this an OPTIMISTIC bound:
    clustering raises P(zero) above the Poisson value.
    """
    return math.exp(-prawn_count_mean(density, area))


def prawn_count_cv(density: float, area: float) -> float:
    """Coefficient of variation of the enclosed prawn count, = 1/sqrt(mean)."""
    mean = prawn_count_mean(density, area)
    if mean <= 0.0:
        return float("inf")
    return 1.0 / math.sqrt(mean)


def replicates_for_cv(density: float, area: float, target_cv: float) -> int:
    """Chambers needed so the MEAN prawn count has CV <= target_cv.

    Counting noise only -- it is a floor on replication, not an estimate
    of total variance. Real scatter also carries sediment heterogeneity
    and per-prawn activity differences, so treat this as a lower bound.
    """
    mean = prawn_count_mean(density, area)
    if mean <= 0.0:
        raise ValueError("no prawns: replication cannot fix a zero mean")
    return math.ceil((1.0 / target_cv) ** 2 / mean)
