# Inputs and units

The runner reads [current.toml](../configs/current.toml). Copy that file for a
different case, then pass its path with `--input`. Unknown input names and
invalid dimensions are rejected. Values use SI units except ppm for CO₂,
µmol/s for net exchange and CFM for the nominal fan rating.

| Group | What to change |
|---|---|
| `chamber` | Inside x/y/z dimensions, floor elevation, maximum gap and assumed physical skin thickness |
| `domain` | Exterior box extents; the Cartesian grid rounds these upward to an even number of cells |
| `air` | Wind speed and direction, temperature/pressure, ambient CO₂ and initial closed-chamber CO₂ |
| `cycle` | Sealed and unsealed durations, travel times and cycle count |
| `palm` | Height, crown radii and trunk scale for the prescribed porous drag and schematic illustration |
| `exchange` | Whole-tree net CO₂ rate applied continuously to leaflets; negative means uptake |
| `fans` | Count, support position, heights, discharge angles, nominal flow, disc size and force width |
| `numerics` | Cell size, maximum time step, save cadence, viscosity, Smagorinsky coefficient and plant drag |
| `waterlily` | Immersion width, minimum numerical wall thickness and native pressure-solver controls |
| `display` | Frame rate and fixed concentration-scale preferences |

`chamber.width_m` is the x dimension along sliding travel. The design's
6 m width is the solver's y dimension (`depth_m`). Heights are measured above
the chamber floor, which is distinct from the external ground.

Wind and fan angles point **toward** the flow direction: 0° is +x and 90° is
+y, counterclockwise viewed from above. These are not meteorological wind-from
bearings. All three fans use the same support x/y position.

Omitting `fans.heights_m` places fans at the centres of equal height bands.
Omitting `fans.azimuths_deg` spaces their directions uniformly around 360°.
To place them explicitly, add these entries inside `[fans]`:

```toml
heights_m = [0.7, 2.0, 3.3]
azimuths_deg = [30.0, 150.0, 270.0]
```

Both lists must match `fans.count`. Remove both entries when changing count
for automatic placement. `count=0` disables fan forcing. `momentum_factor`
scales thrust, not RPM; changing `model` alone changes only the descriptive name.

The net-exchange input comes from February 2025 C2 daytime cycles with
Rg ≥ 10 W/m²: Q95 uptake strength, equivalent to Q05 signed net exchange.
It is not a current-year observation or dynamic leaf photosynthesis.
These inputs affect the CO₂ extension. Molecular diffusivity is in m²/s;
the dimensionless turbulent Schmidt number gives `D_t = nu_t / Sc_t`.
CO₂ is a passive scalar: it does not change air density, buoyancy or stomatal
behaviour. The airflow-only runner deliberately leaves scalar inputs unapplied
and labels its results accordingly.

Changing the chamber dimensions does not rescale plant size or gas exchange.
The plant envelope must fit inside the configured chamber. The numerical walls
extend outward by at least `numerical_wall_cells × cell_m`; the exterior box
must contain them throughout travel. That thickness is reported separately
from physical skin thickness and is a numerical approximation.
