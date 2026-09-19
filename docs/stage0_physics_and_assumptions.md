# Stage 0 — Physics, assumptions, and what I think is wrong with the brief

**Status:** specification only. No transport solver exists yet. Every number below
is either closed-form algebra or a tagged placeholder.

---

## 1. Headline

1. **Your two contested formulae are both correct.** The Elliott & Brooks bed-head
   relation is right as you wrote it, including the 3/8 → 3/2 exponent switch at
   `H/d = 0.34`. The forced-vortex `p(r) = Δp_s (r/R)²` stirrer profile is the right
   functional form. Boudreau's `θ² = 1 − ln(φ²)` is right.
2. **One thing in the brief is genuinely inconsistent:** your governing PDE and your
   analytic penetration depth use *different* definitions of the respiration rate.
   They differ by `√φ` ≈ 0.63 — a 37% error in penetration depth if mixed. Fixed in
   §3.1; the fix is now a test.
3. **Two independent literature sources agree on the pressure forcing to ~30%.**
   Elliott & Brooks, driven with Huettel & Gust's own stated conditions (10 cm s⁻¹
   over a 7 cm ripple), gives 0.13–0.19 Pa cm⁻¹ against their stated 0.2 Pa cm⁻¹.
   That is a real cross-validation, not a coincidence, and it sets the chamber
   stirrer default at **Δp_s = 1.5 Pa**.
4. **Your permeability thresholds are independently reproduced.** At the default
   bed forcing, `k = 3×10⁻¹² m²` gives Péclet 0.28 and `k = 2.6×10⁻¹¹ m²` gives
   Péclet 2.4 across the oxic layer. They bracket Pe = 1 (at `k ≈ 1.1×10⁻¹¹ m²`),
   i.e. they mark "advection stops being negligible" and "advection takes over".
   Since you have no measured `k`, permeability is the model's **primary swept
   axis, not an input**.
5. **The brief is missing seven things**, two of them first-order for your actual
   question. The biggest: **no diffusive boundary layer** (§4.1) and **no
   photosynthesis** (§4.2) — and you have a light chamber, so without the second the
   light deployments cannot be modelled at all.
6. **Engineering result you can use now — and it is the binding constraint.**
   With the real production geometry (168 mm bore, ~600 mm tube, 100 mm inserted
   ⇒ a 0.50 m water column), diffusion-dominated sediment with no prawns draws
   down only **0.205% of initial DO per hour** (0.51 µmol L⁻¹ h⁻¹), and needs
   **97 h to reach 20% drawdown**. **Shorten the water column.** Full analysis in
   §12.

---

## 2. Physics as I will implement it

### 2.1 Domain and conventions

- `z` positive **downward**, `z = 0` at the sediment–water interface (SWI).
- Fluxes positive **into the sediment** (O₂ uptake is a positive flux).
- `p` is the **dynamic** (non-hydrostatic) pressure. Gravity is absorbed into the
  reference state, so it does not appear in Darcy's law. Constant density.

### 2.2 Porous-medium flow

```
q = −(k/μ) ∇p        ∇·q = 0        ⇒   ∇·(k ∇p) = 0
```

`q` is the Darcy flux (volume per unit **bulk** area per time); the interstitial
velocity is `v = q/φ`. Rigid, non-deforming bed. `k` isotropic and, for now,
depth-uniform.

### 2.3 Oxygen transport — corrected

```
φ ∂C/∂t = ∇·(φ D_s ∇C) − q·∇C + φ α(z)(C_w − C) − φ R(C)
```

This is your equation, unchanged. The `φ` in front of `α` and `R` means **both are
defined per unit volume of porewater** (the Boudreau/Meile convention, which is what
published `α` values use). That is now an explicit `respiration.basis` switch in the
config, because it changes the analytic check:

| Basis | Reaction term | Penetration depth `L` | Interfacial flux `J` |
|---|---|---|---|
| **porewater** (default) | `φR` | `√(2 D_s C₀ / R)` | `φ √(2 D_s C₀ R)` |
| bulk | `R_bulk` | `√(2 φ D_s C₀ / R_bulk)` | `√(2 φ D_s C₀ R_bulk)` |

**Your brief wrote the porewater PDE but the bulk `L`.** Both rows are correct
*individually*; mixing them is a `√φ` = 0.63 error. With `R_bulk = φR` the two rows
give identical answers — `tests/test_scaling.py::test_the_two_basis_conventions_agree_when_converted`
asserts exactly that, and a companion test pins the size of the error if you mix them.

Closure identity the solver must satisfy at steady state: `J = φ R L`. Also tested.

Transport closures:
- `D_s = D_0/θ²`, `θ² = 1 − ln(φ²)` — **Boudreau (1996)**, confirmed. At φ = 0.40,
  `θ² = 2.833`, so `D_s = 0.353 D_0`.
- `D_0(T,S)` from the Ramsing & Gundersen tables, Stokes–Einstein scaled.
- `C_sat(T,S)` from Garcia & Gordon (1992). **The implementation gets a stage-1
  acceptance test against a published table value** — I am not going to trust
  recalled coefficients.

### 2.4 Consumption

Default is **Monod**, `R = R₀ C/(K_m + C)`, not pure zero-order. Pure zero-order
drives `C` negative below the oxic layer and formally needs a free-boundary
treatment. With `K_m ≈ 1 µM` against `C_sat ≈ 250 µM`, respiration is
near-zero-order through ~99% of the oxic layer, so the stage-1 analytic check still
applies. Zero-order and first-order remain selectable.

### 2.5 Field bed forcing — verified

```
h_m = 0.28 · (U²/2g) · (H/(0.34 d))^n ,    n = 3/8  if H/d ≤ 0.34
                                           n = 3/2  if H/d >  0.34
```

**Verified against the published form** (Elliott & Brooks 1997 WRR 33:123–136,
fitted to Fehlman 1985). Your version is correct. Three things to be careful about:

- `h_m` is a **head half-amplitude in metres of water**. The pressure amplitude is
  `p_m = ρ_w g h_m`. Your brief jumped straight to `p_m` without the conversion.
- `U` in the calibration is the **mean (depth-averaged) stream velocity**, not the
  free-stream velocity you wrote. Record the measurement height with your velocity.
- Bed BC `p(x) = p_m sin(2πx/λ)`, periodic in `x`. Analytic check for stage 2:
  `p(x,z) = p_m sin(2πx/λ)·exp(−2πz/λ)` for a deep bed, and
  `p = p_m sin(2πx/λ)·cosh(2π(d_b−z)/λ)/cosh(2π d_b/λ)` for finite thickness `d_b`
  on an impermeable base. Both as you specified — correct.

### 2.6 Chamber forcing

`p(r) = Δp_s (r/R)²` is the forced-vortex (solid-body rotation) form — correct.
Sign matters and is a testable prediction: **pressure is lowest at the centre and
highest at the wall**, so the model must produce porewater **downwelling near the
chamber wall and upwelling at the centre**. That is the pattern Huettel & Gust
observed. If stage 5 gives the opposite, stage 5 is wrong.

Realistic range (Huettel & Gust 1992, MEPS 82:187–197): they ran a radial gradient
of **0.2 Pa cm⁻¹**, stated as the equivalent of a 10 cm s⁻¹ bottom flow over a 7 cm
ripple. For `p = Δp_s (r/R)²` the rim gradient is `2Δp_s/R`, so at R = 0.15 m that
setting is **Δp_s = 1.5 Pa**. Useful range **0–5 Pa**.

### 2.7 Chamber water balance

```
V dC_w/dt = −∮ (J_diff + J_adv + J_irr) dA − V R_water
```

Incompressibility plus impermeable walls forces `∮ q·n dA = 0` over the enclosed
SWI, so advection redistributes O₂ without net volume exchange — but it is *not*
flux-neutral, because inflowing water carries `C_w` and outflowing water carries
depleted porewater.

---

## 3. What is wrong in the brief

| # | Item | Verdict |
|---|---|---|
| 3.1 | `L = √(2φD_s C₀/R)` vs the `φR` PDE | **Inconsistent.** `√φ` = 37% error. Fixed, tested (§2.3). |
| 3.2 | "free-stream velocity U" in Elliott & Brooks | **Wrong definition.** Calibration uses mean stream velocity. |
| 3.3 | `p(x) = p_m sin(...)` with `h_m` from the relation | **Missing `ρg`.** `h_m` is head, not pressure. |
| 3.4 | Exponent 3/2 above `H/d = 0.34` | **Correct as written.** Verified. |
| 3.5 | `θ² = 1 − ln(φ²)` | **Correct.** Boudreau (1996). |
| 3.6 | `p(r) = Δp_s(r/R)²` | **Correct form.** Range now anchored (§2.6). |
| 3.7 | `exp(−2πz/λ)` + cosh finite-depth form | **Correct.** |
| 3.8 | Zero-order consumption as the default | **Numerically unsafe.** Monod default instead (§2.4). |

**Applicability caveat on 3.4 (not an error, a limitation):** Elliott & Brooks was
calibrated on 2D triangular **dunes** in a **unidirectional freshwater flume**.
Sandprawn mounds are roughly conical, irregular and three-dimensional, in a tidal
marine lagoon. Using it here is an extrapolation of the *geometry*, and the flume
calibration has no counterpart for biogenic roughness. It is still the best
available a-priori estimator — it is what the marine permeable-sediment literature
uses — but treat `p_m` as uncertain to a factor of ~2, and sweep it.

---

## 4. What is missing from the brief

### G1 — No diffusive boundary layer *(first-order; add in stage 1)*

There is a stagnant film above the SWI across which O₂ moves only by molecular
diffusion. When the bed is diffusion-dominated it can carry a large share of the
**total transport resistance**. Its thickness is set by near-bed shear — **which the
chamber changes by design**. This is a direct chamber-vs-field bias mechanism and it
is absent from your model. For a low-`k` site it may be the *largest* bias term.
Added as `sediment.diffusive_boundary_layer` and, separately,
`chamber.stirrer_dbl_thickness`, so the mismatch can be isolated.

### G2 — No photosynthesis *(first-order; you have a light chamber)*

A light chamber measures **net community metabolism (P − R)**, not respiration. With
`production.enabled = false` the model simply cannot represent it. Added: a
`P(z) = P_max·tanh(I/I_k)·exp(−z/L_light)` source with `L_light ≈ 0.6 mm`, so
production sits in a layer *thinner* than the oxic layer. Two consequences worth
anticipating:

- In permeable sand under light, advection drags O₂-supersaturated water downward,
  and the flux can **reverse sign**.
- Supersaturation in a closed chamber can **bubble out**, which silently destroys the
  mass balance of a real deployment. Flag it in your protocol: log for bubbles.

Also practical: **perspex is transparent, so your "dark" chamber must be actively
shrouded**, and a clear chamber in sun **heats**. Q10 ≈ 2 means a 5 °C rise is a
~40% respiration change. Log temperature inside both variants.

### G3 — No tidal forcing / air exposure

Langebaan is **intertidal**. Water depth `d` in Elliott & Brooks changes through the
tide and the relation is undefined as `d → 0`. Drainage and re-infiltration at low
tide drive porewater exchange that no term in this model represents. **The field
model is a snapshot at one tidal stage**, and the honest comparison is chamber vs
field *at the same tidal stage*.

### G4 — Bioirrigation driving concentration is an upper bound

`α(z)(C_w − C)` drives toward the **overlying water** concentration at every depth.
For a 1 m burrow ventilated in intermittent bouts, the water down there is depleted
and periodically anoxic, so the real driving concentration is well below `C_w`.
Consequence: **the classical form overestimates deep irrigative O₂ supply.**
Both a `reduced_at_depth` option and a `duty_cycle` factor are in the config.

**My assumed α profile, as you asked:**
`α(z) = α₀ exp(−z/z_irr)` with **`α₀ = 5×10⁻⁶ s⁻¹` (158 yr⁻¹)** and
**`z_irr = 0.05 m`**. The exponential mimics the decline in burrow cross-sections
with depth. `z_irr` is *deliberately far shallower than the 1 m burrow depth*, for
two reasons: (i) O₂ is consumed within millimetres of any supplied surface, and
(ii) burrow water is depleted and intermittently ventilated (G4). **I could not find
a species-specific `α` for _K. kraussi_ anywhere I can reach — this is the single
least constrained parameter in the model.** Treat both `α₀` and `z_irr` as fitted,
and let stage 6 rank them.

### G5 — The chamber severs the burrow system

_K. kraussi_ burrows reach **~1 m**. Your chamber inserts 5–20 cm. **It is physically
impossible for this chamber to enclose an intact burrow system** — it cuts every
burrow it encloses. Severed burrows may vent to the chamber, may be plugged by the
prawn, or may short-circuit under the wall edge. This is a *structural* limitation of
the method, not a parameter, and it should appear in your report. The model can
bound it (insertion depth sweep) but cannot resolve it without explicit burrow
geometry.

### G6 — Chamber area drives the replicate variance

At the real 168 mm bore, `A = 0.0222 m²`. At 100 ind m⁻² that is **2.2 prawns per
chamber**, and **11% of chambers catch none at all**; at 10 ind m⁻², **80% are
empty**. Poisson counting noise on prawn number dominates between-chamber scatter and
scales as `1/√(ρA)`, giving a **67% CV at 100 ind m⁻²** — so you need **6 chambers**
just to pull the counting noise alone down to 30%. Poisson is the optimistic case;
real callianassid beds are patchy, which makes it worse. Full table in §12.

### G7 — Smaller omissions, all now in the config

- Water-column respiration inside the chamber, indistinguishable from sediment uptake
  in `dC_w/dt`. A sealed water-only blank measures it.
- **The measured flux is not the in-situ flux even in a perfect chamber**: flux
  depends on `C_w`, so `dC_w/dt` decays and a linear fit under-reads. Config now
  carries `flux_estimator` and `max_drawdown_fraction` so the model can report *what
  you would infer from sensor data* alongside the true instantaneous flux. **That
  difference is one of the answers you are after.**
- Stirring can **resuspend** surface sand, changing `k` and releasing porewater.
- No waves (Langebaan is sheltered — stated as an assumption, not a finding).

---

## 5. Scaling results

Reproduce with `python3 scripts/stage0_report.py`. All from default (largely
**placeholder**) parameters — these set expectations and give stages 1–3 something
falsifiable to hit. They are **not model results**.

**Diffusion baseline** (φ = 0.40, D₀ = 1.8×10⁻⁹ m² s⁻¹, C₀ = 250 µM, R = 10⁻⁴ mol m⁻³ s⁻¹):

| Quantity | Value |
|---|---|
| `θ²` | 2.833 |
| `D_s` | 6.35×10⁻¹⁰ m² s⁻¹ |
| Penetration depth `L` | **1.78 mm** |
| Diffusive flux `J` | **6.2 mmol m⁻² d⁻¹** |
| Damköhler at `z = L` | 2.000 (exact, by construction) |

**Bed forcing** (U = 0.10 m s⁻¹, H = 0.03 m, d = 0.50 m, λ = 0.30 m):
`h_m` = 0.075 mm → `p_m` = **0.75 Pa** → peak `|∇p|` = **0.157 Pa cm⁻¹**;
advective cell depth ≈ λ/2π = 4.8 cm.

**Cross-check** — Elliott & Brooks driven with Huettel & Gust's stated conditions
(0.10 m s⁻¹ over a 0.07 m ripple): `p_m` = 1.03 Pa, `|∇p|` = **0.13–0.19 Pa cm⁻¹**
for λ = 0.35–0.50 m, against their stated **0.2 Pa cm⁻¹**. Two independent published
relations agree to ~30%.

**Where permeability matters** (this is the core of your design question):

| `k` (m²) | Darcy `q` | Péclet over `L` | Regime |
|---|---|---|---|
| 1×10⁻¹² | 0.05 mm h⁻¹ | 0.09 | diffusion |
| **3×10⁻¹²** *(your lower threshold)* | 0.14 mm h⁻¹ | **0.28** | diffusion, advection ~25% |
| 1×10⁻¹¹ | 0.47 mm h⁻¹ | 0.92 | crossover |
| **2.6×10⁻¹¹** *(your upper threshold)* | 1.23 mm h⁻¹ | **2.40** | advection |
| 1×10⁻¹⁰ | 4.74 mm h⁻¹ | 9.24 | advection dominant |

Pe = 1 at `k = 1.1×10⁻¹¹ m²`. **Your thresholds bracket it.** Note this crossover
moves with `U`, `H`, `λ` and `R` — it is not a property of the sediment alone, which
is itself a result worth stating in your report.

**Chamber sensitivity** — superseded by §12, which uses the real geometry.

---

## 6. Parameter inventory

60 parameters, all in `config/config.yaml` with SI units, default, range and source.
`python3 benthic/config.py` prints the audit.

| Status | Count | Meaning |
|---|---|---|
| PLACEHOLDER | 18 | **You must measure or supply this** |
| ASSUMPTION | 15 | My choice, defensible, unverified |
| NUMERICAL | 12 | Solver settings |
| LITERATURE | 7 | Cited source |
| DERIVED | 5 | Computed at runtime |
| MEASURED | 3 | Already known |

**Still completely unset (`value: null`) — the model will refuse to run:**
`sediment.permeability`, `sediment.median_grain_size`, `sediment.sorting_sigma`,
`production.surface_irradiance`, `field.prescribed_head_amplitude`.

`config.get()` raises on an unset placeholder rather than substituting a default.
That is deliberate: it makes "invent a precise value silently" a test failure.

**Ranked by what I expect to matter most, and what it costs you to get:**

| Priority | Parameter | How | Why it ranks here |
|---|---|---|---|
| 1 | `permeability` | constant-head permeameter on intact cores | Sets the entire regime; currently a sweep axis, not an input |
| 2 | `chamber.tube_length` | **a tape measure, today** | Bore is now known (168 mm). Tube length sets the water column, and the water column sets whether the deployment works at all (§12). |
| 3 | `respiration.rate_zero_order` | core incubation or microprofile fit | Sets `L` and the baseline flux |
| 4 | `mound_height`, `mound_wavelength` | ruler + scaled photo | Sets the field advective forcing |
| 5 | `prawn_density` | count burrow openings in a quadrat | Sets `α₀` *and* the replicate variance |
| 6 | `temperature` | logger **inside** the chamber | Q10 ≈ 2; perspex heats |
| 7 | `d50` + sorting | sieve / laser | Only needed if `k` is estimated not measured |
| 8 | `free_stream_velocity`, `water_depth` | ADV + tide record | Note: mean velocity, and record the tidal stage |

---

## 7. Limitations that could change a conclusion

| ID | Simplification | Could flip which conclusion? |
|---|---|---|
| L1 | 2D / axisymmetric | Yes for magnitude. 3D conical mounds focus flow more than 2D ridges; 2D likely **under**estimates advective exchange. |
| L2 | No explicit burrow geometry | Yes. Nonlocal `α` cannot represent a severed burrow venting into the chamber (G5). |
| L3 | Depth-uniform `k`, `φ` | Yes if a compacted or shelly layer exists. Sandprawn mounds are sorted material — near-surface `k` may differ from bulk. |
| L4 | Rigid bed, no resuspension | Matters at high stirring; would show up as apparent `k` change. |
| L5 | Steady bed forcing, no waves/tides | Yes for the field case (G3). Chamber case is much less affected — that asymmetry is itself a bias. |
| L6 | Continuous, time-averaged `α` | Fine for hour-long deployments, wrong for short ones. |
| L7 | `k` isotropic | Probably minor in well-sorted sand. |
| L8 | No temperature transient in the chamber | Could be significant for a clear chamber in sun; measurable, so measure it. |
| L9 | 16 mm wall treated as infinitely thin | The wall displaces and compacts a 16 mm annulus of sand exactly where stirrer-driven downwelling enters. Would bias chamber advection; direction unclear. |
| L10 | Stirrer swirl assumed to reach the bed undiminished | **Likely optimistic at 3:1 aspect ratio.** Sets chamber advective flux; see `stirrer_transfer_efficiency`. |

---

## 8. What I need from you

**Blocking nothing** — I can proceed to stage 1 on placeholders, since stage 1 is an
analytic verification that does not depend on your site values.

**Decisions I'd like before stage 5:**

1. ~~**Chamber geometry.**~~ **Answered (revision 1):** 168 mm bore, 200 mm OD,
   ~600 mm tube. Still need the **tube length measured**, and a decision on whether
   you can shorten the water column (§12).
2. ~~**Confirm the α profile.**~~ **Answered (revision 1):** you have no estimate, so
   the assumed profile in §G4 stands as the default and `α₀` and `z_irr` are swept.
   The model's job is now to report how much they matter; if stage 6 ranks them high,
   §G4 gives the tracer method for measuring `α`.
3. ~~**G1 / G2 / G6 in scope?**~~ **Answered (revision 1):** in. See §12 for exactly
   what "and other stuff" was taken to mean — say so if I have over- or under-read it.
4. **New, and the one that matters:** can the water column be shortened, or must the
   tube stay ~600 mm? This decides whether a deployment is hours or days (§12).

---

## 9. Verified this stage

- `config/config.yaml` parses; all 60 parameters have a valid status; every numeric
  parameter has units; every PLACEHOLDER/ASSUMPTION/LITERATURE entry cites a source;
  every default lies inside its stated range.
- `config.get()` raises on an unset placeholder (permeability specifically).
- Boudreau tortuosity: exact values, monotonicity, open-water limit `θ²(φ=1) = 1`.
- Penetration depth, interfacial flux, and depth-integrated consumption are mutually
  consistent (`J = 2φD_s C₀/L = φRL`); Damköhler = 2 at `z = L`.
- The two respiration bases agree when converted, and the mixed-convention error is
  exactly `√φ`.
- Elliott & Brooks: exponent switch is continuous at `H/d = 0.34`, scales as `U²`,
  raises on `d → 0`, and reproduces the Huettel & Gust ripple gradient.
- Stirrer `Δp_s` ↔ rim-gradient round-trip; chamber and field forcing within a factor
  of 3 at defaults.
- Crossover permeability gives Pe = 1 exactly; your two thresholds bracket it.

**31 tests, all passing** (`python3 -m pytest tests/ -q`).

## 10. Known limitations of this stage

- No PDE is solved. Everything above is algebra.
- `C_sat(T,S)` and `D_0(T,S)` are declared but **not yet implemented**; their
  fallback values are approximate and carry a stage-1 acceptance test against
  published tables.
- The Elliott & Brooks verification is from secondary sources quoting the relation
  verbatim; I could not reach the 1997 paper directly (publisher domains are blocked
  from this environment). The form, coefficient, both exponents and the threshold all
  match your version, and the Huettel & Gust cross-check is independent — but if you
  have library access, **confirm it against the original before you cite it**.
- No species-specific `α` for _K. kraussi_ was found. §G4 is my reasoning, not a
  citation.

## 11. Sources

- Elliott, A.H. & Brooks, N.H. (1997) *Transfer of nonsorbing solutes to a streambed
  with bed forms: Theory.* Water Resources Research 33(1):123–136. (after Fehlman,
  H.M. 1985, MSc thesis, Caltech)
- Huettel, M. & Gust, G. (1992) *Solute release mechanisms from confined sediment
  cores in stirred benthic chambers and flume flows.* Mar. Ecol. Prog. Ser. 82:187–197.
- Huettel, M. & Gust, G. (1992) *Impact of bioroughness on interfacial solute
  exchange in permeable sediments.* Mar. Ecol. Prog. Ser. 89:253–267.
- Huettel, M. & Webster, I.T. (2001) *Porewater flow in permeable sediments.* In: The
  Benthic Boundary Layer. (advection threshold `k` > ~10⁻¹² m²)
- Janssen, F., Faerber, P., Huettel, M., Meyer, V. & Witte, U. (2005) *Pore-water
  advection and solute fluxes in permeable marine sediments (I): Calibration and
  performance of the novel benthic chamber system Sandy.* Limnol. Oceanogr.
  50:768–778.
- Boudreau, B.P. (1996) *The diffusive tortuosity of fine-grained unlithified
  sediments.* Geochim. Cosmochim. Acta 60:3139–3142.
- Garcia, H.E. & Gordon, L.I. (1992) *Oxygen solubility in seawater: better fitting
  equations.* Limnol. Oceanogr. 37:1307–1312.
- Ramsing, N.B. & Gundersen, J.K. *Seawater and Gases: tabulated physical
  parameters.* Unisense technical report.
- Meile, C., Koretsky, C.M. & Van Cappellen, P. (2001) *Quantifying bioirrigation in
  aquatic sediments: an inverse modeling approach.* Limnol. Oceanogr. 46:164–177.
- Branch, G.M. & Pringle, A. (1987) *The impact of the sand prawn Callianassa kraussi
  Stebbing on sediment turnover and on bacteria, meiofauna, and benthic microflora.*
  J. Exp. Mar. Biol. Ecol. 107:219–235.
- Pillay, D. & Branch, G.M. (2011) *Bioengineering effects of burrowing thalassinidean
  shrimps on marine soft-bottom ecosystems.* Oceanogr. Mar. Biol. Annu. Rev. 49:137–192.
- Flemming, B.W. *Depositional processes in Saldanha Bay and Langebaan Lagoon,
  Western Cape, South Africa.* (site sedimentology; bimodal quartz/carbonate sands)

---

## 12. Revision 1 — real chamber geometry, and what it costs you

**Inputs that changed:** production chamber is **168 mm bore / 200 mm OD**, tube
**~600 mm** ("maybe slightly longer, definitely less than 1 m"), and you have **no
estimate for α**. Scope confirmed: diffusive boundary layer **in**, plus the rest of
the recommended set (§12.4).

### 12.1 Headline: the chamber is 3.3× less sensitive than my stage-0 placeholder

The enclosed area is **0.0222 m²**, not the 0.0707 m² I assumed — 31% of it. But
**area is not the problem**, because it cancels:

```
V dC_w/dt = −J·A ,   V = A·h    ⇒    dC_w/dt = J/h
```

**The signal depends only on the water column height `h`, not on the bore.** With a
600 mm tube inserted 100 mm, `h = 0.50 m` against the 0.15 m I had assumed.

| | Stage 0 placeholder | **Production chamber** |
|---|---|---|
| Bore | 300 mm | **168 mm** |
| Water column `h` | 0.15 m | **0.50 m** |
| Enclosed area | 0.0707 m² | **0.0222 m²** |
| Volume | 10.6 L | **11.1 L** |
| Aspect ratio `h`/bore | 0.5 : 1 | **3.0 : 1** |
| Diffusion-only drawdown | 1.71 µmol L⁻¹ h⁻¹ | **0.51 µmol L⁻¹ h⁻¹** |
| … as % of saturation | 0.68 %/h | **0.205 %/h** |
| Time to 20% drawdown | 29 h | **97 h** |

At a 900 mm tube it is worse still: `h = 0.80 m`, 0.32 µmol L⁻¹ h⁻¹, **156 h**.

### 12.2 Signal-to-noise: a short deployment will not work in the diffusive case

Over a 4 h deployment the diffusion-only signal is **2.1 µmol L⁻¹**. Against a
0.4 µmol L⁻¹ sensor resolution and 0.5 µmol L⁻¹ h⁻¹ of drift, that is a
**signal-to-noise of 0.86** — you would be fitting a slope to drift.

This is the diffusion-only floor, and it is the honest worst case. With advection at
`k ≥ 10⁻¹¹ m²` plus prawns, total flux could be 3–10× higher and a 3–4 h deployment
becomes workable. **But you do not yet know `k`, so you cannot currently rule out the
floor case** — and a chamber designed only for the favourable case will return
unusable data on the patches that matter most.

**What actually fixes it — shorten the water column:**

| Water column `h` | Drawdown | Time to 20% | Aspect ratio |
|---|---|---|---|
| 0.50 m *(as built)* | 0.51 µmol L⁻¹ h⁻¹ | 97 h | 3.0 : 1 |
| 0.30 m | 0.86 µmol L⁻¹ h⁻¹ | 58 h | 1.8 : 1 |
| **0.20 m** | **1.28 µmol L⁻¹ h⁻¹** | **39 h** | **1.2 : 1** |
| 0.10 m | 2.57 µmol L⁻¹ h⁻¹ | 19 h | 0.6 : 1 |

This needs no new fabrication if the lid can seat lower in the tube, or if the tube
is inserted deeper — though deeper insertion trades against G5 (burrow severance) and
changes the flow path under the wall. **A ~0.20 m water column is the sweet spot**:
2.5× the signal, and it brings the aspect ratio back inside the envelope the stirred-
chamber literature actually calibrated.

### 12.3 Two new problems the tall tube creates

**Aspect ratio is outside the calibrated envelope.** Huettel & Gust (1992) and Janssen
et al. (2005) both ran ~1.6 : 1 (30 cm tall, 19 cm ID). The production chamber is
**3.0 : 1**, up to 4.8 : 1 at a 900 mm tube. Two assumptions weaken:

- **Swirl may not reach the bed.** The forced vortex is generated at the stirrer and
  decays with distance. Over 3–5 bore diameters, the radial pressure gradient at the
  sediment surface is **lower than the stirrer setting — possibly far lower**. I have
  added `stirrer_transfer_efficiency` (default 1.0, the *optimistic* bound, so the
  model overstates rather than understates chamber advection until calibrated). **This
  is now a top-3 uncertainty and it is specific to your geometry.** Calibrate with a
  dye or brine release over an impermeable plate, or a centre-to-wall differential
  pressure measurement at the bed.
- **The well-mixed assumption is no longer a formality.** A 0.5 m tall, 0.168 m bore
  tube, clear, in sun, can stratify thermally. If it does, your DO sensor reads a
  local value, not `C_w`, and the flux is wrong regardless of everything else in this
  model. **Cheap check: a second DO sensor at a different height, or a dye test.**

**`Δp_s` rescaled.** The Huettel & Gust anchor is a pressure *gradient* (0.2 Pa cm⁻¹),
so the smaller bore needs a proportionally smaller centre-to-wall difference:
`Δp_s = 0.2 Pa cm⁻¹ × R/2` = **0.84 Pa**, down from 1.5 Pa.

### 12.4 Prawn-count lottery, at the real bore

`A = 0.0222 m²`. Assuming Poisson-distributed burrows:

| Prawn density | Mean per chamber | P(chamber is empty) | CV | Chambers for CV ≤ 30% |
|---|---|---|---|---|
| 10 ind m⁻² | 0.22 | **80%** | 212% | 51 |
| 50 ind m⁻² | 1.11 | 33% | 95% | 11 |
| 100 ind m⁻² | 2.22 | 11% | 67% | **6** |
| 200 ind m⁻² | 4.43 | 1% | 48% | 3 |

Poisson is the **optimistic** null model — real callianassid beds are patchy, which
raises both P(empty) and CV. **Count burrow openings at your site before fixing the
replicate number**; it is the cheapest measurement on the list and it determines
whether 6 chambers or 20 is the right answer.

Note this is counting noise *alone* — a floor on replication, not an estimate of total
variance. Sediment heterogeneity and per-prawn activity add to it.

### 12.5 What "and other stuff" was taken to mean

Confirmed in scope. **Say so if I have over- or under-read this:**

| Gap | In? | Where it lands |
|---|---|---|
| G1 diffusive boundary layer | **yes** | Stage 1 — changes the 1D model directly |
| G2 photosynthesis / light chamber | **yes** | Stages 1 & 5 |
| G6 prawn-count variance | **yes** | Stage 6 (scaffolding already in `scaling.py`) |
| G7 flux-estimator bias | **yes** | Stage 5 |
| G4 α uncertainty | **swept, not assumed** | Stages 4 & 6 |
| G3 tidal forcing | **no** — out of scope | Stated as a limitation (L5) |
| G5 burrow severance | **bounded, not resolved** | Insertion-depth sweep, stage 5 |

### 12.6 What changed in the repo

- `config/config.yaml`: chamber section rewritten — `bore`, `outer_diameter`,
  `wall_thickness`, `tube_length`, derived `height_above_sediment` and `water_volume`,
  `aspect_ratio_max`, `stirrer_transfer_efficiency`, `well_mixed_assumption`,
  `sensor_resolution`, `sensor_drift`, `replicates`. `stirrer_pressure` 1.5 → 0.84 Pa.
  `deployment_duration` 1 h → 4 h. 69 parameters now.
- `benthic/scaling.py`: chamber area/volume, drawdown rate, time-to-drawdown, aspect
  ratio, signal-to-noise, and the Poisson prawn-count functions.
- `tests/`: 47 passing (was 31). The headline numbers above are pinned as tests, so
  they cannot drift silently as parameters change.

### 12.7 Three things to do before Stage 5

1. **Measure the tube length**, and tell me whether the water column can be shortened.
2. **Count burrow openings** in a few quadrats — sets both `α₀` and the replicate number.
3. **Get your DO sensor's short-term stability and drift** (not its absolute accuracy).
   With this water column, drift is the same size as the signal.
