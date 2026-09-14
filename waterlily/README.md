# Whole-tree chamber in WaterLily

Reproducible Julia simulation of a sliding whole-tree chamber with three central
fans and CO₂ uptake on representative palm leaflets. **WaterLily 1.8.0 computes
airflow; a separately tested Julia extension computes CO₂ transport.** Python
only renders saved results.

The reference has 4 × 6 × 4 m inside airspace (96 m³), a 2.39 m palm, 1.5 m/s wind
and two cycles of 5 minutes sealed followed by 10 minutes unsealed. Fans operate
while sealed. Prescribed net uptake continues in every phase.

This is an **engineering diagnostic, not an experimentally validated chamber
model**. Numerical wall thickening and first-order transport affect mixing.
Conservation and coupling checks are separate from those physical limits.
The coarse example **fails the airflow-to-scalar coupling screen during wall
motion** and is not grid-converged. Use it to inspect the implementation and
operating sequence, not to judge real fan adequacy or chamber refresh time.

## Run and view the result

Use Julia 1.12 (tested with 1.12.5), Python 3.12 and FFmpeg. Start in this folder:

```bash
julia --project=. --threads=1 -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

# Two complete cycles: 1,800 physical seconds.
julia --project=. --threads=1 scripts/run_co2.jl \
  --input configs/current.toml --output results/my_current_chamber

python3.12 -m venv .venv-render
source .venv-render/bin/activate
python -m pip install -r requirements-render.lock
python scripts/render_co2.py results/my_current_chamber
python scripts/audit_bundle.py results/my_current_chamber
```

For a quick execution check, append `--duration 10` and choose another output
folder. That check stays sealed; tests exercise movement separately. Every run
requires a fresh output folder and preserves its inputs, source and environment.

`LocalPreferences.toml` selects WaterLily's SIMD CPU backend. Julia and BLAS use
one thread for this small grid. Registered dependency pins in the manifest do
not depend on another local WaterLily checkout.

## Results

The compact repository example includes [the animated cycle](examples/current_chamber/videos/c2_202502_cycle.gif),
[the MP4](examples/current_chamber/videos/c2_202502_cycle.mp4), all four media
views, 24 stills and the numerical audits. See [its contents](examples/current_chamber/README.md).

The full example is written to `results/current_chamber/`:

| Location within a result | Contents |
|---|---|
| `videos/c2_202502_cycle.mp4` and `.gif` | Palm, concentration slice, mean CO₂ and spatial variation |
| `videos/fan_sections.*` | Concentration and native airflow at the three fan heights |
| `videos/operating_timeseries.*` | Concentration, cycle, fan and exchange histories |
| `videos/leaf_source_location.*` | Representative leaflets and their source allocation |
| `figures/` | Six physical-time stills for each of the four views |
| `timeseries.csv` | Mean, spatial standard deviation, extrema and operation |
| `transport_audit.csv` | Every-step mass balance, geometric conservation and airflow corrections |
| `fields/` | Concentration, fluid volume, native velocity, distance and leaflet area arrays |
| `summary.json`, `inputs.json` | Run parameters, outputs and numerical evidence |
| `configs/source/`, `manifest.json` | Frozen generating source and SHA-256 checksums |

All figures have white backgrounds. CO₂ uses fixed viridis, normally 380–405 ppm;
the renderer expands the limits if saved concentrations exceed them. It never
changes the concentration data. `progress.json` shows the last saved physical
time during execution.

The native airflow check is separate in `results/diagnostics/airflow_only/`.
Reproduce it with `scripts/run.jl` and `scripts/render_airflow.py`; it contains
no CO₂. Existing Python results remain in their own folders.

## Change the experiment

Copy [configs/current.toml](configs/current.toml), edit it and pass the new path
with `--input`. Chamber dimensions, plant envelope, wind, fan count and placement,
thrust, cycle timing and exchange are independent inputs.
The [input guide](docs/INPUTS.md) explains units and conventions.

The reference rate is **−44.40047791098563 µmol CO₂/s**, from February 2025 C2
daytime Q95 uptake strength with Rg ≥ 10 W/m². It is prescribed net exchange,
not a 2026 measurement or a dynamic photosynthesis model. The 240 source
triangles represent leaflets; trunk and bare frond axes have no source.

[Equations and numerical method](docs/METHOD.md) ·
[Inputs](docs/INPUTS.md) · [Verification and porting limits](docs/PORTING.md)
[Measured numerical checks](docs/VERIFICATION.md)

Short sensitivity checks can be reproduced with
`julia --project=. --threads=1 scripts/check_sensitivity.jl results/my_sensitivity.json`.
They hold the numerical wall thickness and actual exterior box fixed when
refining the grid. They do not establish convergence over the complete cycle.
