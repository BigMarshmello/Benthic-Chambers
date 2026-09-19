# Benthic Chambers

Process model for dissolved-oxygen flux across the sediment–water interface in
permeable intertidal sand, built to answer one engineering question:

> **Which transport processes dominate O₂ exchange, and how much does the chamber
> itself — enclosure plus stirring — bias the measured flux relative to the
> undisturbed bed?**

Site: sandprawn (*Kraussillichirus kraussi*) flats, Langebaan Lagoon, South Africa.

## Processes represented

1. Molecular diffusion in porewater, with Boudreau (1996) tortuosity
2. Topography-driven advection — bed pressure over mounds drives Darcy porewater flow
3. Stirrer-driven advection in a closed chamber (radial pressure gradient)
4. Bioirrigation by burrowing sandprawns (nonlocal, depth-dependent)
5. O₂ consumption (zero-order / first-order / Monod)
6. Photosynthesis, for light-chamber deployments *(added at stage 0)*

## Configurations

- **FIELD** — 2D Cartesian (x, z), periodic, repeating mounds, free-stream flow
- **CHAMBER** — 2D axisymmetric (r, z), stirred, closed volume, partial wall insertion

## Status

| Stage | Description | State |
|---|---|---|
| 0 | Physics, assumptions, parameter inventory | **done** — [note](docs/stage0_physics_and_assumptions.md) |
| 0.1 | Revision 1: real chamber geometry | **done** — [§12](docs/stage0_physics_and_assumptions.md#12-revision-1--real-chamber-geometry-and-what-it-costs-you) |
| 1 | 1D diffusion + zero-order consumption | not started |
| 2 | 2D Darcy pressure and flow field | not started |
| 3 | Darcy ⇄ O₂ transport coupling | not started |
| 4 | Bioirrigation | not started |
| 5 | Chamber mode, closed-chamber `C_w(t)` | not started |
| 6 | Sensitivity analysis, Morris screening | not started |
| 7 | Figures and report section | not started |

## Layout

```
config/config.yaml   every parameter: SI units, default, range, status, source
benthic/config.py    loader + validator; refuses unset placeholders
benthic/scaling.py   closed-form relations every later stage must reproduce
scripts/             runnable reports
tests/               pytest suite, one group per stage
docs/                one markdown note per stage
```

## Use

```bash
pip install -r requirements.txt
python3 benthic/config.py          # parameter audit: what is measured vs assumed
python3 scripts/stage0_report.py   # scaling report
python3 -m pytest tests/ -q        # test suite
```

## Conventions

- `z` positive **downward**, `z = 0` at the sediment–water interface
- Fluxes positive **into the sediment** (O₂ uptake is positive)
- `p` is **dynamic** pressure; gravity is absorbed into the reference state
- `R(C)` and `α(z)` are per unit volume of **porewater** (they enter multiplied by `φ`)

## Data status

Permeability, grain size, irradiance, tube length and sensor noise are **not yet
measured**. `config.get()` raises rather than substituting a value, so nothing
downstream can quietly invent them. See §6 of the stage 0 note for the measurement
priority list, and §12.7 for the three that block Stage 5.

Chamber geometry as built: **168 mm bore, 200 mm OD, ~600 mm tube**. The resulting
0.50 m water column makes the diffusion-only signal marginal — see §12.
