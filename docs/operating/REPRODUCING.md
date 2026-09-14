# Running a simulation

All commands below run from the repository root. Inputs, the schematic plant
mesh and Python dependencies are included; the example does not require an
external plant database.

## Install

Use Python **3.12** and FFmpeg with `ffprobe`, libx264 and GIF support.
The reference environment used Python 3.12.13 and FFmpeg 8.1.

```bash
python3.12 -m venv .venv-operating
source .venv-operating/bin/activate
python -m pip install -r requirements-operating.lock
ffmpeg -version
ffprobe -version
python scripts/operating_chamber.py validate
```

On Windows, activate the environment with `.venv-operating\Scripts\activate`
and use the same `python` commands. The dependency lock includes the tests.
The launcher limits numerical libraries to one thread per worker and records
package versions with each run.

## Start with a short run

```bash
python scripts/operating_chamber.py all \
  --input examples/quick_start.yaml --output runs/quick_start
```

This runs one four-second cycle on a coarse grid, including opening and
closing, and generates the four video views. It checks that the calculation
and rendering work together; it is too short to study mixing performance.

## Reproduce the reference example

```bash
python scripts/operating_chamber.py all \
  --input examples/central_fans/input.yaml \
  --output runs/central_fans --workers 2
```

The main calculation covers 1800 s with 0.5 m target cells and 0.5 s maximum
time steps. Two additional runs cover the first 300 s: a 1/3 m grid and a
0.25 s time step. These check sensitivity to numerical resolution. The full
workflow takes several minutes on the reference machine, including rendering.

Open `runs/central_fans/README.md` for links to the generated results.
The [committed example](../../examples/central_fans/) contains the media,
histories and checks for comparison. The full CO₂/velocity fields and frozen
source tree are generated in the run folder.

Choose a new output directory for each run. Input validation rejects an existing
folder, invalid units or values, and a plant mesh that does not fit the chamber.
Use [INPUTS.md](INPUTS.md) to change dimensions, wind, fans or net exchange.

## Run individual stages

```bash
python scripts/operating_chamber.py prepare \
  --input configs/operating/current.yaml --output runs/my_case
python scripts/operating_chamber.py run --output runs/my_case --workers 2
python scripts/operating_chamber.py render --output runs/my_case
python scripts/operating_chamber.py audit --output runs/my_case
```

`prepare` freezes inputs, source code, geometry and environment information.
The remaining commands read those frozen inputs. Do not pass `--input` or
`--set` to `run`, `render` or `audit`.

Completed numerical cases are checksum-verified and skipped when resuming.
An interrupted, incomplete case requires a new output folder; retain the old
folder for inspection. If source files have changed since preparation, use
the snapshot in the run's `configs/source/` directory:

```bash
python runs/my_case/configs/source/scripts/operating_chamber.py render \
  --output runs/my_case
```

The snapshot includes its own dependency lock. Installing that lock gives the
closest match to the original environment. For another complete run, use
`configs/inputs.json` from a result as `--input`; JSON is also accepted.

`render --preview` produces stills only. A preview audit leaves video checks
marked as not checked. A full audit raises an error on integrity or
conservation failures. Grid and time-step sensitivity failures are recorded
separately; inspect `verification/checks.json` even when the command succeeds.

## Check the installation and published files

```bash
python -m pytest tests/operating tests/recreated
python scripts/check_example.py
python scripts/check_example.py --decode
```

The example checker verifies the committed file list and SHA-256 hashes.
`--decode` also reads every image and movie with FFmpeg. It checks file
integrity, not the physical validity of the calculation.

Compare reproduced time histories and numerical checks before comparing video
bytes. Encoders and plotting libraries can change file bytes across platforms
without changing the simulated fields. The reference
[provenance record](../../examples/central_fans/provenance.json) identifies
inputs, source hashes, environment and the files supplied with the example.

## Connect another model

`operating.config.compile_study()` resolves dimensions, timing, wind, fans and
net exchange. The numerical entry point is
`fan_solver.run_case(cfg, destination, mesh_factory=..., fans=...)`.
It returns a summary and writes histories plus snapshots containing time,
wall gap, cell coordinates and volumes, CO₂ in ppm, and velocity in m/s.
Leaflet-source maps record intersected leaf area and exchange per cell.

Another flow solver needs an adapter for this coordinate, unit and output
contract. The plotting code can then read its saved fields. Its conservation
and convergence checks must be established separately. A dynamic photosynthesis
model likewise needs an exchange adapter; the current `exchange.model` accepts
only `constant_leaflet_net_exchange` and does not run XPalm or another plant model.
