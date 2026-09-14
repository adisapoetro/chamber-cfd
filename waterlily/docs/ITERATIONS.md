# Audit and improvement record

The simulation remains a diagnostic. A numerical change is accepted only after
the affected tests and coupled inventory checks pass. The 10% airflow-to-scalar
correction screen remains unchanged, including failures at moving-wall endpoints.

## Iteration 1: baseline

Repository revision `109263776b209c2804c7cef177542071720a1127` contains the first
complete 1,800 s WaterLily/CO₂ bundle. Its local full result is preserved in
`results/archive/iteration_01/`, with inputs, frozen code, fields and checksums.
It passed 320 implementation assertions and conserved CO₂ inventory, but its
largest chamber-region flux correction was 163.81%. The largest aperture
velocity correction was 134.61 m/s. The coarse field was not grid-converged.

## Iteration 2: timing, stagnant diffusion and 3D presentation

The original mapping combined the new native velocity with midpoint geometry.
Moving faces now pair both endpoint velocities with geometry at their respective
times and integrate the reconstructed fluid velocity against the changing area.
A face supported at only one endpoint uses that velocity, with this extension
recorded in the audit. Native WaterLily momentum and pressure are unchanged.

Stagnant far-field faces now exchange CO₂ by diffusion. An independent analytical
diffusion problem with ambient reservoirs verifies this boundary treatment and
its refinement behaviour. Manufactured moving-face tests verify time pairing,
opening/closing support and the area–velocity integral. All 337 Julia assertions
pass, including the full accelerated coupled cycle.

The cycle illustration now includes three perpendicular CO₂ sections and sampled
three-component airflow arrows. A separate large 3D view uses the same saved
fields. White backgrounds, fixed viridis, the simple palm and leaflet uptake
are retained. Coloured sections are not a reconstructed continuous volume.
Fan-level panels show the fan assigned to that saved height; grey arrows mark
inactive fans. The review also removed a shared colour-scale callback problem
that produced plotting-library cleanup errors during long exports.

The shortened moving-wall probe improves the final closing-step relative
correction from about 173% to 92%, but still requires about 134 m/s correction
in a vanishing aperture. The first opening step worsens from about 58% to 102%:
the timing correction does not improve every interface state. The probe uses
real 20 s travel times with only one second
sealed and one second fully open, so it is a coupling diagnostic, not an
operating-cycle result. The full 1,800 s rerun lowers the worst relative
correction from 163.81% to 87.92%, with both versions still failing the screen.
Its native velocity and geometry fields remain byte-identical. Full-cycle
results and audit comparisons are recorded in [RESULTS.md](RESULTS.md).

The geometry-only check isolates an interface discontinuity. At a 1 µm gap,
the sampled diffuse-kernel integral exceeds the exactly closed value by
1.367 m³ on the 0.5 m grid and 1.061 m³ on the 0.25 m grid, with the same 1.5 m
numerical walls and exterior box. Exact sharp fluid volume remains constant.
At the roof seam the kernel mobility approaches 0.5, then becomes zero on
closure. These sampled integrals are diagnostic measures of kernel support,
not WaterLily fluid mass or measured chamber volume. They identify a remaining
geometry problem; halving the cell size has not removed it.

Reproduce the separate probes from this project folder, using fresh paths:

```bash
julia --project=. --threads=1 scripts/check_coupling.jl results/my_coupling.json
julia --project=. --threads=1 scripts/check_gap_geometry.jl results/my_gap_geometry.json
julia --project=. --threads=1 scripts/check_sensitivity.jl results/my_sensitivity.json
```

[Coupling evidence](../evidence/coupling_iteration_02.json) ·
[Gap-geometry evidence](../evidence/gap_geometry_iteration_02.json) ·
[Refinement evidence](../evidence/short_sensitivity.json)

To compare the two local complete bundles:

```bash
python scripts/compare_iterations.py results/archive/iteration_01 \
  results/current_chamber results/my_iteration_comparison
```

## Remaining qualification

The unresolved closing gap and thick numerical shells are still material
limitations. Conservation alone cannot establish realistic mixing. A qualified
interface treatment, fixed-geometry refinement, wall/domain sensitivity and
measured fan/airflow comparisons are needed before using this model for chamber
performance decisions. Test outcomes and example media do not establish those
scientific claims.
