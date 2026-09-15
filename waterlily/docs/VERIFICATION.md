# Numerical evidence for this implementation

The Julia package test run passes **337 assertions** under Julia 1.12.5 and
WaterLily 1.8.0. See [the test record](../evidence/tests.json) and
[the executable tests](../test/runtests.jl). These are implementation checks,
not experimental validation.

The [15 September engine audit](../../docs/ENGINE_AUDIT.md) repeats this suite,
adds independent analytical 3D benchmarks, verifies the official package tree,
and replays the full reference fields. Its dated evidence supplements the
earlier records below; the chamber qualification failures remain open.

New tests cover zero-wind diffusion to ambient reservoirs and moving-face
velocity/geometry timing. The [iteration record](ITERATIONS.md) includes the
retained baseline and reproducible endpoint probes.

Short sealed runs use the same actual 14 × 10 × 10 m exterior box and hold
numerical wall thickness at 1.5 m. They apply identical plant, wind, fans and
net exchange. [Machine-readable results](../evidence/short_sensitivity.json)
include runtime, time-step count, inventory errors and coupling corrections.

| Ten-second case | CO₂ spatial SD | Change from reference | Largest chamber-region flux correction |
|---|---:|---:|---:|
| 0.5 m cells, maximum Δt 0.5 s | 0.261153 ppm | — | 63.51% |
| 0.5 m cells, maximum Δt 0.025 s | 0.262327 ppm | +0.45% | 85.41% |
| 0.25 m cells, maximum Δt 0.5 s | 0.350737 ppm | +34.30% | 0.590% |

All three means agree with analytical sealed-volume depletion to within
1.1 × 10⁻¹² ppm. The coarse reference's correction is largest during startup;
at ten seconds it is 1.607%. The refined case passes the 10% coupling screen
over this short interval. The large change in spatial variation means the
coarse chamber result is **not grid-converged**.

The maximum time step is a cap; WaterLily's native CFL and the explicit
momentum forcing can require smaller steps. The time-step comparison used
129 versus 400 actual steps over ten seconds.

A separate one-second crosswind case with zero fans executes conservatively
but fails the relative correction screen. With very little physical internal
flow, that relative measure is sensitive to the initial wall mapping. Failure
remains recorded rather than being removed from the evidence.

The full-cycle result carries its own `verification.json`, `summary.json` and
`transport_audit.csv`. Its fields are independently replayed to check inventory
against the integrated source and external boundary flux. Media decoding and
white backgrounds are checked separately.

The initial full coarse run failed the coupling screen at motion endpoints:
the chamber-region correction was approximately 49% when opening finished
and 164% when closing finished. The endpoint-paired mapping reduces those
first-cycle values to about 22% and 88%; they still fail. This is not just a
startup issue. Exact sharp fluid
volumes and the native diffuse wall field behave differently as a narrow gap
closes. The numerically thick shells also produce an exaggerated squeeze-flow
region. Conservation is retained, but the corrected scalar flux cannot be
treated as a small perturbation of native airflow there.

These endpoint failures prevent claiming that the coarse example quantitatively
recreates the physical chamber. Refinement and wall/interface qualification must
address them before using its spatial gradients or refresh behaviour as evidence.

The next resolution check is a complete cycle at 0.25 m with fixed numerical
wall thickness, using [refined_fixed_wall.toml](../configs/refined_fixed_wall.toml):

```bash
julia --project=. --threads=1 scripts/run_co2.jl \
  --input configs/refined_fixed_wall.toml --output results/refined_fixed_wall
```

This is substantially more expensive than the coarse case. It has not been
executed over the full 1,800 s record. Even a converged fixed-wall result would
still require a physical wall/domain sensitivity study and field validation.
