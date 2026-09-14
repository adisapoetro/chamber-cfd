# Reference, verification and limits

The operating reference is `chamber-cfd` revision
`b31facf08685fcd6134e3fc349f22cbd46889363`: nominal 96 m³ airspace, three central
fans, two 900 s cycles and constant leaflet net exchange. Its portable input
copy is [reference_python_input.yaml](../configs/reference_python_input.yaml).
This Julia project is independent of the original Python solver.

| Component | Implementation |
|---|---|
| Momentum, pressure, immersed boundaries | Native WaterLily 1.8.0 |
| Shell motion, fan forcing, plant drag | Adapters using native hooks |
| Smagorinsky momentum diffusion | Native `sgs!` helper |
| CO₂ advection, diffusion, moving fluid inventory | Separate Julia cut-cell extension |
| Whole-tree exchange allocation | Representative leaflet surface intersections |
| Rendering | Python reads saved Julia fields; no Python flow solver |

The released WaterLily core does not supply the full moving-wall scalar
treatment. The [upstream scalar proposal, PR #290](https://github.com/WaterLily-jl/WaterLily.jl/pull/290),
checked on 14 September 2026, remains unmerged at
`de68dc235cb000a6dac217381213b3cbf2cf1165`. This project does not depend on that
experimental implementation.

## Verification

Run `julia --project=. --threads=1 -e 'using Pkg; Pkg.test()'`.

Tests cover independent sinusoidal advection and diffusion solutions with grid
refinement; native Taylor–Green vortex decay; bounds and sink exhaustion;
impermeable walls; exact sealed-volume depletion; uniform and nonuniform
inventory conservation during opening and closing; source area and integrated
uptake; and a complete accelerated cycle coupling native airflow to transport.
Accelerated cycles test implementation, not real chamber performance.

Output evidence distinguishes:

- **Airflow:** native pressure convergence, divergence and wall-core velocity.
- **Transport:** inventory balance, geometric law, bounds and coupling correction.
- **Reproduction:** input/source hashes, package pins and media decoding.

Read `summary.json` and `transport_audit.csv` beside each video. A completed
execution can still fail a coupling or numerical-sensitivity screen.

## Remaining scientific limits

A 0.5 m mesh cannot resolve a millimetre-scale skin. Its three-cell numerical
wall extends 1.5 m outward, substantially altering external blockage. Physical
wall thickness is an assumption. Mapping native BDIM velocity to sharp scalar
apertures requires an audited correction. A larger exterior domain, refinement
and an appropriate thin-wall treatment remain necessary for quantitative mixing
comparisons.

First-order upwinding adds numerical diffusion. Benchmark convergence establishes
code behaviour, not convergence of this chamber case. Cell size, time step,
external domain, wall representation and subgrid coefficients need sensitivity
studies before inferring refresh time or fan adequacy.

The source is February 2025 C2 daytime Q95 net uptake, based on 188 cycles with
Rg ≥ 10 W/m² and no upper radiation cutoff. It is not a 2026 observation, a leaf
gas-exchange law or a calibrated plant model. Leaf surfaces and porous drag are
schematic; installed fan performance and support position are unverified.

Passing numerical tests does not establish experimental validity. Previous
Python outputs, raw observations, plant-model checkouts and archives are unchanged.
