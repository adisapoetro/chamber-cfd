# CFD engine and scientific audit

Audit date: **15 September 2026**. Baseline repository revision:
[`191a1aa`](https://github.com/adisapoetro/chamber-cfd/tree/191a1aa542cc92283fa61f27090f37145d823f3b).
Exact source hashes, fresh test results and numerical measurements are in
[evidence.json](audit/2026-09-15/evidence.json).

**The WaterLily implementation uses real, unmodified WaterLily for 3D airflow.
The Python implementation uses project-owned CFD code, without PhiFlow.
Neither reference is presently qualified to predict physical chamber mixing,
refresh time or fan adequacy quantitatively.**

This audit covers source and dependency identity, units and assumptions,
numerical tests, analytical benchmarks, saved-result integrity, conservation,
resolution and moving-wall coupling. It does not supply experimental validation
or a separate external scientific review. Publication/privacy checks are a
different audit and provide no evidence of physical correctness.

## Which engine generated which result?

| Workflow and result | Airflow and pressure | CO₂ |
|---|---|---|
| Root `scripts/operating_chamber.py`; `examples/central_fans`; local `c2_central_fans` | Project-owned Python finite-volume solver on an arbitrary Lagrangian–Eulerian (ALE) moving mesh; NumPy arrays and SciPy linear algebra | Project-owned conservative Python transport and leaflet source |
| `waterlily/scripts/run_co2.jl`; `waterlily/examples/current_chamber`; full `waterlily/results/current_chamber` | Official WaterLily 1.8.0 native 3D momentum, pressure and immersed-boundary operations | Project-owned Julia finite-volume extension, including the airflow-to-scalar interface |

The published Python workflow imports
[`operating/workflow.py`](../src/fluid_dynamic/operating/workflow.py), which calls
[`recreated/fan_solver.py`](../src/fluid_dynamic/recreated/fan_solver.py) and
[`operators.py`](../src/fluid_dynamic/recreated/operators.py). Its tests execute
in an environment without PhiFlow installed. NumPy and SciPy are numerical
libraries; they do not make this a PhiFlow calculation. Older archived projects
may have different dependencies and are outside this reference audit.

The Julia call chain is
[`build_flow` → `advance_flow!`](../waterlily/src/Flow.jl) →
`WaterLily.Simulation` → `WaterLily.sim_step!` → native momentum and multigrid
pressure. `CheckedPoisson` delegates to WaterLily's `MultiLevelPoisson`, sets
its tolerance/iteration limit and checks its residual. It does not implement a
replacement pressure algorithm. Chamber geometry, fan forcing, drag and the
choice of boundary conditions are supplied by this project.

The loaded package's complete Git tree hash is
`8e1d973f428df4bae450e1a45f19d0dba3ac6857`. It matches both the pinned
[`Manifest.toml`](../waterlily/Manifest.toml) and the official
[`v1.8.0` release](https://github.com/WaterLily-jl/WaterLily.jl/tree/v1.8.0),
whose commit is `2a2b2ba2eeacc885d84e346f743b72a34e0b713f`.
The audit used Julia 1.12.5 and the SIMD CPU backend.

WaterLily is a published incompressible-flow solver with its own benchmark and
validation evidence. That evidence concerns the package and the problems
studied in its paper; it does not validate this chamber adapter or CO₂ extension.
See Weymouth and Font (2025),
[WaterLily.jl: A differentiable and backend-agnostic Julia solver for incompressible viscous flow around dynamic bodies](https://www.sciencedirect.com/science/article/pii/S0010465525002504).

## What is standard, and what needs qualification?

Both implementations use recognizable incompressible-flow and conservative
advection–diffusion formulations. Air motion transports CO₂; molecular and
modelled unresolved diffusion spread gradients. Fans supply momentum, and the
prescribed leaflet source removes CO₂. Opening does not reset the chamber to
ambient concentration. The [Python method](operating/METHOD.md) and
[Julia method](../waterlily/docs/METHOD.md) describe the equations and units.

The WaterLily CO₂ extension follows
[`advance_scalar!`](../waterlily/src/Coupled.jl) →
[`waterlily_face_flux`](../waterlily/src/Transport.jl) →
[`compatible_flux`](../waterlily/src/CutCells.jl) → `transport_step`.
It maps WaterLily's diffuse immersed velocity onto sharp fluid apertures, then
corrects those fluxes to conserve volume as the walls move. The deblending
omits the BDIM first-moment contribution and is therefore approximate.
CO₂ uses the corrected flux; the displayed arrows show native WaterLily velocity.
Their agreement cannot be assumed where the correction is large.

The correction metric is

$$r_Q=\frac{\lVert Q_{corrected}-Q_{mapped}\rVert_2}
{\max(\lVert Q_{mapped}\rVert_2,10^{-30})},$$

evaluated on faces touching the fixed chamber region. The existing 10% screen
is an implementation safeguard, not a universal physical-accuracy tolerance.
A small correction would still require convergence and experimental checks.
The current large correction fails even this preliminary screen.

The Julia scalar scheme is first-order upwind with backward Euler and adds
numerical diffusion. Numerical smoothing can therefore resemble physical
mixing. Python and Julia also use different exterior boundary and wall
treatments; their matching input values do not make them an otherwise identical
solver comparison.

## Checks repeated in this audit

| Check | Fresh result | What it establishes |
|---|---|---|
| Python component suite | **79 tests pass** | Tested operators, configuration, conservation and workflow behaviour |
| Julia component suite | **337 assertions pass** in 13 test sets | Native flow integration, scalar operators, geometry and a short accelerated coupled cycle |
| Native 3D periodic ABC flow, 8³ → 16³ → 32³ | Relative velocity L₂ error **1.358% → 0.378% → 0.0959%**; divergence below 2 × 10⁻⁸ s⁻¹ | Refinement toward a known analytical incompressible solution with all three velocity components |
| Julia 3D periodic scalar wave, 8³ → 12³ → 16³ | RMS error **0.0457 → 0.0345 → 0.0273 ppm**; bounded and conservative | Analytical advection and diffusion in all three directions without a moving body |
| Python saved fields | **243 snapshots and 243 source maps**, across main/grid/time-step cases | File integrity, recorded statistics, source allocation and sealed-volume balance |
| WaterLily full bundle | **809 manifested files**, 181 snapshots, 25,686 transport steps | Source integrity, saved statistics and CO₂ inventory replay against logged source/boundary flux |

The new native benchmark uses the decaying periodic Beltrami field
`u = (sin z + cos y, sin x + cos z, sin y + cos x) exp(−νt)` on `[0, 2π]³`,
with `ν = 0.05` and `t = 0.1 s`. The nonlinear acceleration is a pressure
gradient and the Laplacian is `−u`, so this is an analytical solution of the
full incompressible equations. The separate scalar benchmark uses a periodic
sine wave with known translation and diffusive decay, evaluated as cell means.
Neither test contains a chamber, fan, canopy or moving interface.

The Python suite includes a manufactured viscous channel profile, but that test
isolates the viscous operator. This audit does not establish a full coupled 3D
canonical-flow benchmark for the custom Python momentum algorithm. Likewise,
the Julia suite's accelerated cycle is not a test of realistic refresh time.

Independent WaterLily inventory replay differs from the logged integrated
source and external flux by at most **6.85 × 10⁻⁹ ppm m³**. This checks saved
concentration/volume against logged fluxes; it does not independently reconstruct
every native face flux. For Python, the external scalar flux is not saved
separately, so the open-phase global budget cannot be independently reconstructed
from the bundle. Its solver-reported residual and saved-field statistics were
checked, and the first sealed mean agrees with its analytical value within
**1.47 × 10⁻⁷ ppm**.

All nine current Julia core files match the full reference's frozen source
byte for byte. The Python numerical operators and workflow also match. Three
Python `src/` files differ: one configuration comment and presentation text in
two rendering modules. Configuration syntax trees are identical; the actual
byte mismatches remain recorded in the evidence instead of being hidden.

## Findings that prevent quantitative use

| Finding | Evidence and consequence |
|---|---|
| **WaterLily-to-CO₂ coupling fails** | The full reference reaches **87.92%** chamber-region flux correction at 900 s, against the 10% screen. Its advecting scalar flux is materially altered from the native-flow mapping during closure. |
| **WaterLily spatial mixing is not grid-converged** | At 10 s, halving cells from 0.5 to 0.25 m changes CO₂ spatial SD from **0.2612 to 0.3507 ppm (+34.30%)**, with the same domain and numerical wall thickness. The smaller time-step case changes SD by only +0.45%; that does not resolve the grid dependence. No refined full 1,800 s cycle was executed in this audit. |
| **The closing interface is discontinuous** | A geometry-only probe from a 1 µm gap to exact closure changes the sampled diffuse-kernel integral by **1.367 m³** at 0.5 m cells and **1.061 m³** at 0.25 m cells. Sharp fluid volume remains constant. These are kernel-support diagnostics, not physical air-volume or mass losses. |
| **Numerical walls differ greatly from the input skin** | At 0.5 m cells the native numerical wall is **1.5 m**, compared with a **0.003 m assumed skin**. Nominal inside volume is retained, but outside blockage and the narrow closing-gap flow are altered. The actual exterior box is 14 × 10 × 10 m after grid rounding. |
| **Python also fails its grid screen** | The first 300 s spatial-SD comparison changes by **5.834% relative RMS**, exceeding the existing 5% threshold. Its time-step comparison is 1.124%. Opening-grid convergence remains unestablished. |
| **Installed airflow and geometry remain unverified** | No chamber tracer-decay or anemometry validation is supplied. Fan thrust, mounting, clear dimensions, canopy drag and leakage are not qualified by field measurements here. |

The largest auxiliary aperture-velocity correction is **133.95 m/s** in a
vanishing aperture. This is a correction diagnostic, not a predicted native
chamber airspeed. It must not be used as a physical velocity or dismissed merely
because the inventory balance still closes. A separate one-second crosswind,
zero-fan probe also fails the relative coupling screen; near-zero internal flow
makes that ratio sensitive, and the failure remains visible.

## Inputs and biological interpretation

The reference assumes **4 × 6 × 4 m inside airspace in solver x/y/z coordinates
(96 m³)**. Design width/depth are 6/4 m; x follows sliding travel. The floor is
0.6 m above ground. These are nominal inputs, not a newly surveyed installation.
There are two cycles of 300 s sealed and 600 s unsealed; each unsealed phase
includes 20 s opening and 20 s closing travel.

Three 0.12 m fans act at 0.667, 2.000 and 3.333 m above the floor, at a support
0.30 m beside the trunk, directed toward 0°, 120° and 240°. They run only while
sealed. Their force uses a nominal free-air rating of 148.3 CFM per fan and
`F = ρQ²/A`, about 0.513 N each. This does not establish installed thrust or
flow. A 0.5 m cell cannot resolve the 0.12 m fan; the Gaussian momentum-source
width is 0.3 m. Fan housings and support blockage are omitted.

The specified whole-tree net exchange is **−44.40047791098563 µmol/s** in all
phases. Its recorded provenance is February 2025 C2 daytime Q95 uptake strength,
using Rg ≥ 10 W/m² among 188 qualifying cycles. The audit verifies the value
used in the computation, not a fresh reanalysis of the private observations.
It is not a 2026 measurement, gross photosynthesis or a responsive XPalm/VPalm
physiology calculation.

The palm is a 2.39 m schematic envelope with 240 representative leaflet
triangles. Their approximately 0.733 m² proxy area distributes the total source;
it is not a measured whole-palm leaf area. Trunk and bare frond axes have no
assigned uptake. Temperature, humidity, buoyancy and responsive stomata are
not solved.

For a sealed volume at fixed pressure and temperature, the mean must obey

$$\bar C_{ppm}(t)=\bar C_{ppm}(0)+q_{tree}t\frac{RT}{PV}.$$

With the stated rate this gives a **3.39461 ppm decrease in 300 s**. The result
follows from the imposed rate and volume, irrespective of fan mixing. Matching
that curve verifies a balance; it does not verify concentration gradients,
fan performance or biological response.

## Reproduce the audit

Install the pinned Python environment using the [running guide](operating/REPRODUCING.md).
From the repository root, run the component tests:

```bash
python -m pytest tests/operating tests/recreated -q
julia --project=waterlily --threads=1 -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

# Native package identity, 3D analytical airflow and 3D scalar transport.
julia --project=waterlily --threads=1 waterlily/scripts/audit_backend.jl \
  waterlily/results/backend_audit.json

# Short matched-domain sensitivity and geometry-only closing-gap probes.
mkdir -p waterlily/results
julia --project=waterlily --threads=1 waterlily/scripts/check_sensitivity.jl \
  waterlily/results/short_sensitivity_audit.json
julia --project=waterlily --threads=1 waterlily/scripts/check_gap_geometry.jl \
  waterlily/results/gap_geometry_audit.json
```

These JSON paths must be fresh. The analytical benchmark script exits with an
error if its implementation checks fail. The sensitivity and geometry scripts
record diagnostic results; a successful process does not mean convergence.

For a saved-field replay, first generate full bundles with the documented
[Python](operating/REPRODUCING.md) and [Julia](../waterlily/README.md) workflows.
The compact committed examples omit the large binary fields and cannot serve
as substitutes for full run folders. For example, after generating
`runs/central_fans` and `waterlily/results/my_current_chamber`:

```bash
python scripts/audit_saved_fields.py \
  --python-run runs/central_fans \
  --waterlily-run waterlily/results/my_current_chamber \
  --output runs/saved_field_audit.json
```

This independent reader imports NumPy, not either solver. It writes a separate
report without modifying the bundles. `implementation_replay: PASS` means its
integrity and balance checks passed; resolution/coupling failures and
`physical_qualification: NOT_ESTABLISHED` are reported separately. Do not run
it with Python's `-O`, which disables assertions. Fresh runs may have additional
frozen documentation/helper files and different timings or roundoff while
retaining the same numerical core.

## Work required before physical claims

1. Qualify a consistent moving-wall representation and its airflow-to-scalar
   mapping. Check closing gaps, swept volume and flux corrections through real
   20 s travel, with native momentum and pressure preserved.
2. Establish full-cycle grid and time-step convergence, then assess numerical
   wall thickness, exterior-domain size and transport diffusion separately.
   Resolve or calibrate the fan momentum model against installed performance.
3. Verify installed dimensions and operating conditions, then compare with
   independent airspeed and tracer/CO₂ recovery measurements using declared
   acceptance criteria. Use separately documented exchange inputs for any new
   period or biological claim.

The next implementation task is the closing-interface qualification in step 1.
This audit adds reproducible checks and documentation; it changes no existing
solver operators, reference inputs, fields, videos or prior evidence records.
