"""Print every number quoted in docs/stage0_physics_and_assumptions.md.

Run:  python3 scripts/stage0_report.py
Nothing here is a model result -- it is closed-form algebra on the
default (largely PLACEHOLDER) parameters, used to set expectations and
to give stages 1-3 something falsifiable to hit.
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from benthic import config, scaling as sc  # noqa: E402

SEC_PER_DAY = 86400.0


def main() -> None:
    cfg = config.load()

    phi = config.get(cfg, "sediment.porosity")
    d0 = config.get(cfg, "sediment.d0_o2")
    c0 = config.get(cfg, "site.o2_saturation")
    rate = config.get(cfg, "respiration.rate_zero_order")
    mu = config.get(cfg, "site.dynamic_viscosity")
    rho = config.get(cfg, "site.water_density")
    g = config.get(cfg, "site.gravity")

    u = config.get(cfg, "site.free_stream_velocity")
    depth = config.get(cfg, "site.water_depth")
    height = config.get(cfg, "field.mound_height")
    lam = config.get(cfg, "field.mound_wavelength")

    bore = config.get(cfg, "chamber.bore")
    radius = bore / 2.0
    dps = config.get(cfg, "chamber.stirrer_pressure")
    tube = config.get(cfg, "chamber.tube_length")
    d_ins = config.get(cfg, "chamber.insertion_depth")
    h_ch = tube - d_ins

    theta2 = sc.tortuosity_squared(phi)
    d_s = sc.sediment_diffusivity(d0, phi)

    print("=" * 66)
    print("STAGE 0 SCALING REPORT  (closed form, default parameters)")
    print("=" * 66)

    print("\n-- sediment diffusion -------------------------------------")
    print(f"  porosity phi                 {phi:.3f}")
    print(f"  tortuosity theta^2           {theta2:.4f}   [1 - ln(phi^2)]")
    print(f"  D_0 (O2, seawater)           {d0:.3e} m2/s")
    print(f"  D_s = D_0/theta^2            {d_s:.3e} m2/s  ({d_s/d0:.3f} x D_0)")

    print("\n-- 1D diffusion + zero-order consumption (stage 1 target) --")
    L = sc.penetration_depth(d_s, c0, rate)
    J = sc.diffusive_flux(phi, d_s, c0, rate)
    L_bulk = sc.penetration_depth_bulk(phi, d_s, c0, rate * phi)
    print(f"  C_0 (air saturation)         {c0*1e3:.0f} umol/L")
    print(f"  R (porewater basis)          {rate:.2e} mol m-3 s-1"
          f"  = {rate*SEC_PER_DAY:.2f} mol m-3 d-1")
    print(f"  penetration depth L          {L*1e3:.3f} mm")
    print(f"  diffusive flux J             {J:.3e} mol m-2 s-1"
          f"  = {J*SEC_PER_DAY*1e3:.2f} mmol m-2 d-1")
    print(f"  same L via bulk basis        {L_bulk*1e3:.3f} mm   (identical, as it must be)")
    print(f"  Damkohler at z = L           {sc.damkohler(rate, L, d_s, c0):.3f}  (must be 2)")

    print("\n-- FIELD: Elliott & Brooks bed forcing ---------------------")
    h_m = sc.elliott_brooks_head_amplitude(u, height, depth, g)
    p_m = sc.head_to_pressure(h_m, rho, g)
    grad = sc.sinusoid_max_gradient(p_m, lam)
    print(f"  U={u} m/s, H={height} m, d={depth} m, lambda={lam} m")
    print(f"  H/d                          {height/depth:.3f}"
          f"   -> exponent {'3/8' if height/depth <= 0.34 else '3/2'}")
    print(f"  head half-amplitude h_m      {h_m*1e3:.4f} mm")
    print(f"  pressure amplitude p_m       {p_m:.3f} Pa")
    print(f"  peak |dp/dx|                 {grad:.2f} Pa/m = {grad/100:.3f} Pa/cm")
    print(f"  advective cell depth ~lam/2pi {lam/(2*math.pi)*1e2:.1f} cm")

    print("\n-- cross-check vs Huettel & Gust (1992) --------------------")
    print("  They state 0.2 Pa/cm is the gradient from a 10 cm/s flow over")
    print("  a 7 cm ripple. Elliott & Brooks, same forcing, d = 0.5 m:")
    for lam_hg in (0.35, 0.50):
        h = sc.elliott_brooks_head_amplitude(0.10, 0.07, 0.5, g)
        p = sc.head_to_pressure(h, rho, g)
        gr = sc.sinusoid_max_gradient(p, lam_hg)
        print(f"    lambda={lam_hg:.2f} m -> p_m={p:.3f} Pa, |dp/dx|={gr/100:.3f} Pa/cm")
    print("  => independent sources agree to within ~30%. Both formulae stand.")

    print("\n-- CHAMBER: stirrer forcing --------------------------------")
    rim = sc.stirrer_rim_gradient(dps, radius)
    dps_hg = sc.stirrer_pressure_from_rim_gradient(20.0, radius)
    area = sc.chamber_area(bore)
    print(f"  bore {bore*1e3:.0f} mm, R = {radius:.4f} m   (area {area:.5f} m2)")
    print(f"  delta_p_s (default)          {dps:.2f} Pa")
    print(f"  rim gradient 2*dp_s/R        {rim:.2f} Pa/m = {rim/100:.3f} Pa/cm")
    print(f"  delta_p_s matching H&G 0.2 Pa/cm at the rim: {dps_hg:.2f} Pa")
    print("  sign: low p at centre, high at wall -> downwelling at the wall,")
    print("        upwelling at the centre. Falsifiable in stage 5.")

    print("\n-- advection vs diffusion: where does k matter? ------------")
    k_cross = sc.crossover_permeability(mu, phi, d_s, grad, L)
    print(f"  crossover permeability (Pe=1 over L): {k_cross:.2e} m2")
    print(f"  your quoted thresholds:               3.0e-12 and 2.6e-11 m2")
    print("  Pe over the oxic layer at the default bed forcing:")
    for k in (1e-12, 3e-12, 1e-11, 2.6e-11, 1e-10):
        q = sc.darcy_velocity(k, grad, mu)
        pe = sc.peclet(q, L, phi, d_s)
        print(f"    k={k:8.1e} m2 -> q={q*3.6e6:8.3f} mm/h,  Pe={pe:7.3f}"
              f"   {'advection-controlled' if pe > 1 else 'diffusion-controlled'}")

    print("\n-- irrigation number --------------------------------------")
    alpha = config.get(cfg, "bioirrigation.alpha0")
    print(f"  alpha0 {alpha:.1e} s-1 = {alpha*3.156e7:.0f} yr-1")
    print(f"  I = alpha*L^2/D_s at z=0     {sc.irrigation_number(alpha, L, d_s):.2e}")
    print("  (tiny at the oxic-layer scale: irrigation does not reshape the")
    print("   microprofile, it adds a parallel supply over the burrow depth)")

    print("\n-- chamber signal size (REVISION 1: real geometry) ---------")
    vol = sc.chamber_volume(bore, h_ch)
    dur = config.get(cfg, "chamber.deployment_duration")
    dcdt = sc.drawdown_rate(J, h_ch)
    ratio = sc.aspect_ratio(h_ch, bore)
    limit = config.get(cfg, "chamber.aspect_ratio_max")
    print(f"  tube {tube:.2f} m - insertion {d_ins:.2f} m -> water column {h_ch:.2f} m")
    print(f"  V = {vol*1e3:.1f} L, A = {area:.5f} m2, V/A = h = {h_ch:.3f} m")
    print(f"  aspect ratio h/bore          {ratio:.2f} : 1", end="")
    print(f"   *** exceeds calibrated {limit:.1f}:1 ***" if ratio > limit else "")
    print(f"  diffusion-only dC_w/dt       {-dcdt*3600*1e3:.3f} umol L-1 h-1"
          f"  ({dcdt*3600/c0*100:.3f} %/h)")
    print(f"  time to 20% drawdown         "
          f"{sc.time_to_drawdown(J, h_ch, c0, 0.20)/3600:.0f} h")
    print(f"  drawdown over {dur/3600:.1f} h            "
          f"{dcdt*dur/c0*100:.2f} % of initial")
    res = config.get(cfg, "chamber.sensor_resolution") * 1e-3
    drift = config.get(cfg, "chamber.sensor_drift") * 1e-3 / 3600.0
    print(f"  signal/noise over {dur/3600:.1f} h        "
          f"{sc.signal_to_noise(J, h_ch, dur, res, drift):.2f}"
          f"   (need >~3; diffusion-only case)")
    print("  NOTE: DIFFUSION-ONLY. Advection at k >= 1e-11 plus prawns")
    print("        gives several times this. But the tall water column is")
    print("        the binding constraint -- see 'if you shorten the tube'.")

    print("\n-- if you shorten the tube ---------------------------------")
    for h_alt in (0.50, 0.30, 0.20, 0.10):
        r_alt = sc.drawdown_rate(J, h_alt)
        print(f"    water column {h_alt:.2f} m -> {r_alt*3600*1e3:6.3f} umol L-1 h-1"
              f",  {sc.time_to_drawdown(J, h_alt, c0, 0.20)/3600:5.0f} h to 20%"
              f",  aspect {sc.aspect_ratio(h_alt, bore):.1f}:1")

    print("\n-- prawn-count lottery (enclosed area is small) ------------")
    print(f"  A = {area:.5f} m2")
    for dens in (10.0, 50.0, 100.0, 200.0):
        print(f"    {dens:5.0f} ind/m2 -> mean {sc.prawn_count_mean(dens, area):5.2f}"
              f" prawns, P(zero)={sc.prawn_zero_probability(dens, area)*100:5.1f}%"
              f", CV={sc.prawn_count_cv(dens, area)*100:5.1f}%"
              f", need {sc.replicates_for_cv(dens, area, 0.30):2d} chambers for CV<=30%")
    print("  Poisson is optimistic: real prawn beds are patchy, which")
    print("  raises P(zero) and CV above these values.")
    print()


if __name__ == "__main__":
    main()
