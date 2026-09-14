# Reproduce or change the operating-chamber study

Run these commands from the `fluid-dynamic` repository root. A fresh checkout
needs only repository files: no PalmTwin database, prior output, notebook,
absolute user path, Julia model, Blender file or private source table.

## Install the tested environment

Use Python **3.12** (reference 3.12.13) and FFmpeg/ffprobe on PATH (reference 8.1,
with libx264 and GIF support). The runtime dependency lock includes tests.

```bash
python3.12 -m venv .venv-operating
.venv-operating/bin/python -m pip install -r requirements-operating.lock
.venv-operating/bin/python scripts/operating_chamber.py validate
.venv-operating/bin/python -m pytest tests/operating tests/recreated/test_conservation.py
ffmpeg -version
```

If using uv: `uv venv --python 3.12 .venv-operating`, then
`uv pip install --python .venv-operating/bin/python -r requirements-operating.lock`.
The legacy `uv.lock`/project dependencies support earlier PhiFlow workflows;
`requirements-operating.lock` is the explicit environment for this workflow.
The launcher sets BLAS/OMP/MKL to one thread. Each output records library and
platform versions, configuration hashes and every executed source file.

## Run the complete bundle

```bash
.venv-operating/bin/python scripts/operating_chamber.py all \
  --input configs/operating/current.yaml \
  --output runs/03_reference_chamber/reproduction --workers 2
```

A fresh directory is mandatory. The active result stays at
`runs/03_reference_chamber/c2_central_fans`; its previous 3.7 m version is
preserved under `archive/2026-09-14-before-size-revision/c2_central_fans`.
The default geometry is now nominal 4 × 6 × 4 m (96 m³), adopted by the owner
on 14 September 2026. It runs 1800 s with
0.5 m target cells and 0.5 s maximum steps, plus a 300 s 1/3 m grid screen and
a 300 s 0.25 s time-step screen. CFD takes a few minutes on the reference machine;
rendering the four videos takes longer. Do not confuse movie duration with
simulated physical time.

The new output contains:

```text
README.md, METHOD.md, INPUTS.md, REPRODUCING.md
configs/             resolved user inputs, solver cases, frozen source, environment
geometry/            original illustration and representative leaflet source
videos/              cycle, fan_sections, operating_timeseries, leaf_source_location
                     each as MP4 and GIF
figures/             stills, CO2 history, schedule, source and fan schematics
data/cfd_runs/       main, grid and dt: histories, saved fields, source-area maps
verification/        numerical/media checks, hashes and full bundle manifest
logs/                one simulation log per case
```

The four views match the chosen bundle's subjects. The revised active bundle
and new runs use generic filenames; the archived 3.7 m bundle retains its
historical `c2_202502_cycle` / `scenarios_timeseries` names and original bytes.

For a short engineering check:

```bash
.venv-operating/bin/python scripts/operating_chamber.py all \
  --set cycle.closed_s=2 --set cycle.open_phase_s=2 \
  --set cycle.opening_travel_s=0.5 --set cycle.closing_travel_s=0.5 \
  --set cycle.cycles=1 --set numerics.cell_m=1 \
  --set numerics.dt_s=0.25 --set numerics.save_every_s=0.5 \
  --set screens.enabled=false \
  --output runs/03_reference_chamber/smoke
```

This short check exercises opening and closure; it is not equivalent to the
30-minute operating study. For changed dimensions/wind/fans/rate see INPUTS.md.
`validate` prints resolved solver inputs and checks dimensions/fan placement;
`prepare` additionally checks the actual leaflet mesh before creating output.

## Resume without changing frozen inputs

```bash
python scripts/operating_chamber.py prepare --output runs/03_reference_chamber/my_run
python scripts/operating_chamber.py run --output runs/03_reference_chamber/my_run
python scripts/operating_chamber.py render --output runs/03_reference_chamber/my_run
python scripts/operating_chamber.py audit --output runs/03_reference_chamber/my_run
```

Complete numerical cases are hash-checked and skipped when resuming. An incomplete
case directory is never overwritten: preserve it and start a new output. Changed
inputs require `prepare`/`all` in another folder. A code change after preparation
fails the source check; rerun with its frozen source or prepare a new study.
For stills only, add `--preview` to render/all. A preview audit explicitly leaves
media unchecked; it does not claim a complete video delivery. An audit raises an
error on conservation/integrity failures. A failed sensitivity screen is recorded
as FAIL separately from a successful engineering execution; inspect checks.json.

For exact source recovery, `configs/source/` is a portable source tree containing
the launcher, numerical dependencies, mesh and dependency lock. Install its lock,
then run its launcher on the existing frozen bundle with `run` or `render`.
For a new rerun, copy `configs/inputs.json` to a YAML/JSON input and use `all`;
the methods and default input are also included in the frozen source tree.
JSON is accepted as YAML. Bit-identical video/ZIP bytes are not promised across
FFmpeg/NumPy versions or platforms. Compare numerical arrays and histories within
documented tolerance, keeping runtime versions alongside the comparison.

## Adapter boundary for another model

`operating.config.compile_study()` resolves units, geometry, timing, fans and net
exchange. The current numerical entry is `fan_solver.run_case(cfg, destination,
mesh_factory=..., fans=...)`. It returns a summary and writes time histories plus
snapshots (`time_s`, `gap`, cell centres/faces/volumes, CO2 ppm, velocity m/s).
The leaflet source adapter supplies cell-intersection areas. Renderers consume
saved fields, and audits independently recompute fixed-region statistics and
source totals. A different flow model must write that same coordinate/unit/time
contract or provide an explicit translator; it must not reuse this solver's
verification status. A dynamic photosynthesis model would need a new named
exchange adapter and conservation tests; `exchange.model` currently rejects
unsupported models. No implicit XPalm coupling is claimed.

## Git and recovery

Push source, documentation, dependency pins, tests and explicit diagnostic inputs.
Results, archives, raw observations, CAD files and virtual environments remain
local. Git reproduces results; it is not an off-machine backup of their videos
or the private observation provenance. The archive guide records physical moves
and the local maintenance record verifies their hashes. Never force-push or add
all workspace changes indiscriminately.
