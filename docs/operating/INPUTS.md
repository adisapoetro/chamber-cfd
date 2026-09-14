# Operating-chamber inputs

Edit `configs/operating/current.yaml`, or copy it and give the copy to `--input`.
Every field is optional with a validated default; the checked-in file spells out
all defaults. Unknown names, wrong types, infinities and invalid geometry fail.
Use SI units except CO₂ in ppm, exchange in µmol/s and the stated fan CFM rating.

| Group / variable | Default | Meaning |
|---|---:|---|
| `schema_version` | 1 | Input contract version |
| `title` | Operating chamber… | Figure title; use a short description |
| `chamber.width_m` | 4 | Closed interior x dimension; shells slide along x |
| `chamber.depth_m` | 6 | Closed interior y dimension |
| `chamber.height_m` | 4 | Nominal airspace height, excluding the elevated floor |
| `chamber.floor_m` | 0.6 | Floor elevation above exterior ground |
| `chamber.maximum_gap_m` | 2 | Gap between halves; each half travels gap/2 |
| `domain.x_m`, `y_m`, `z_m` | null | Null gives width + gap + 8, depth + 4, floor + height + 4.7 m; these are exterior extents, not chamber size |
| `air.wind_speed_m_s` | 1.5 | Imposed external horizontal wind speed; zero = calm |
| `air.wind_to_deg` | 0 | Flow **toward** this direction: 0° = +x, 90° = +y; counterclockwise, no compass bearing |
| `air.ambient_co2_ppm` | 400 | Incoming ambient concentration |
| `air.initial_co2_ppm` | 400 | Initial concentration inside the closed chamber |
| `air.temperature_k` | 298.15 | Constant air temperature for gas conversion/density |
| `air.pressure_pa` | 101325 | Constant reference pressure |
| `cycle.closed_s` | 300 | Sealed interval at the start of each cycle |
| `cycle.open_phase_s` | 600 | Entire unsealed phase, **including both travel intervals** |
| `cycle.opening_travel_s` | 20 | Opening travel time |
| `cycle.closing_travel_s` | 20 | Closing travel time, ending at the period boundary |
| `cycle.cycles` | 2 | Number of complete cycles; default period = 900 s |
| `palm.height_m` | 2.39 | Schematic plant envelope above the chamber floor |
| `palm.crown_radius_x_m`, `crown_radius_y_m` | 1.52 | Horizontal envelope radii, not measured leaf area |
| `palm.trunk_radius_m` | 0.144831 | Trunk illustration/porous drag scale |
| `exchange.net_co2_umol_s` | −44.40047791098563 | **Whole-tree net** exchange: negative removes CO₂; positive releases it; zero disables the scalar source |
| `exchange.label` | Prescribed February 2025… | Provenance/assumption label shown on all videos; update when changing the rate |
| `exchange.model` | constant_leaflet_net_exchange | Only currently supported forcing adapter |
| `fans.count` | 3 | Nonnegative integer; zero means no fans |
| `fans.pole_x_m`, `pole_y_m` | 0.3, 0 | Fixed support position relative to the trunk at (0,0) |
| `fans.heights_m` | null | Null = centres of equal height bands; explicit list length must match count; measured above floor |
| `fans.azimuths_deg` | null | Null = evenly spaced around 360°; explicit list length must match count |
| `fans.first_azimuth_deg` | 0 | First automatic direction, in the same flow-to convention as wind |
| `fans.model` | Delta AFC1212DE | Descriptive model name; changing this alone does not change physics |
| `fans.diameter_m` | 0.12 | Disc diameter used to calculate nominal momentum |
| `fans.free_air_cfm_per_fan` | 148.3 | Assumed free-air volume flow **per fan**, not verified installed flow |
| `fans.momentum_factor` | 1 | Multiplier on nominal momentum forcing; 0 turns forcing off; not an RPM fraction |
| `fans.kernel_sigma_m` | 0.3 | Gaussian force regularization width; ≥ fan radius |
| `fans.schedule` | closed_only | Only supported schedule; no forcing while opening/open/closing |
| `numerics.cell_m` | 0.5 | Target maximum block-cell size; fitted cells need not all have this width |
| `numerics.dt_s` | 0.5 | Maximum time step; steps split at saved times and wall-motion events |
| `numerics.save_every_s` | 10 | Saved-field cadence; final state also saved |
| `numerics.duration_s` | null | Optional shorter diagnostic duration; null runs all configured cycles |
| `numerics.viscosity_m2_s` | 1.5e−5 | Molecular kinematic viscosity |
| `numerics.molecular_diffusivity_m2_s` | 1.6e−5 | CO₂ molecular diffusion |
| `numerics.smagorinsky` | 0.12 | Smagorinsky coefficient |
| `numerics.turbulent_schmidt` | 0.7 | Eddy viscosity / scalar eddy diffusivity ratio |
| `numerics.canopy_drag_per_m` | 0.3 | Prescribed porous canopy drag coefficient |
| `numerics.trunk_drag_per_m` | 20 | Prescribed trunk drag coefficient |
| `screens.enabled` | true | Run finer-grid and smaller-time-step diagnostics |
| `screens.duration_s` | null | Null = first closure, limited to main duration |
| `screens.grid_refinement` | 1.5 | Divide cell size by this factor |
| `screens.timestep_refinement` | 2 | Divide time step by this factor |
| `screens.relative_tolerance` | 0.05 | Relative RMS screen threshold; engineering criterion, not field validation |
| `display.co2_min_ppm`, `co2_max_ppm` | 380, 405 | Preferred fixed linear viridis scale; expanded outward if actual extrema exceed it |
| `display.fps` | 6 | Encoded frames per second; physical simulation time is printed separately |

## Geometry decision, 14 September 2026

The owner adopted the recent documentation's nominal **6 m design width ×
4 m design depth × 4 m clear height**. Solver x is the sliding direction, so
the configuration lists `width_m=4`, `depth_m=6`, `height_m=4`: **96 m³**.
The August 28 LIBZ budget request describes an existing 6 × 4 × 4 m chamber;
the September upsizing presentation labels the same size as approximate
width × depth × clear height. These are documentary dimensions, not an
inside-to-inside field survey. The previous 3.7 m height came from the CAD
glazing envelope and is preserved with its earlier result.

The floor remains a prescribed 0.6 m above ground, placing the new roof at
4.6 m. Automatic fan heights are **0.667, 2.000 and 3.333 m above the floor**;
the automatically sized external domain is 14 × 10 × 9.3 m. The palm and
prescribed net uptake retain their previous values. This geometry change
does not reprocess or rescale the historical observation-derived uptake.

The local documentary evidence and source hashes are retained in
`simulation/maintenance/2026-09-13-chamber-size-inspection/`. Those private
documents are not required by, or included in, the portable source repository.

Example overrides (quote lists and text at the shell):

```bash
python scripts/operating_chamber.py all --input configs/operating/current.yaml \
  --set chamber.width_m=5 --set chamber.height_m=4.2 \
  --set air.wind_speed_m_s=2 --set air.wind_to_deg=90 \
  --set fans.count=4 --set fans.free_air_cfm_per_fan=120 \
  --set exchange.net_co2_umol_s=-30 \
  --set 'exchange.label=Prescribed diagnostic net uptake' \
  --output runs/03_reference_chamber/my_changed_inputs
```

For three explicitly placed fans use, for example,
`--set 'fans.heights_m=[0.7,1.8,3.0]' --set 'fans.azimuths_deg=[30,150,270]'`.
Reset both lists to null when changing count for automatic placement.
Chamber changes do not silently rescale the palm. Actual source triangles must
fit inside the chamber. Small palms may be incompatible with the fixed
representative leaflet template; a clear validation failure requires a geometry
adapter, rather than silently relocating uptake to the trunk.

Changing fan flow by a factor r changes the nominal force by r². There is no
measured RPM-to-flow or installed pressure-flow curve in this model. Supplying a
fan RPM without that relation would not define a reproducible physical input.
