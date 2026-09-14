# Operating whole-tree chamber CFD

**Current result:** [c2_central_fans](runs/03_reference_chamber/c2_central_fans/).
Open the [cycle video](runs/03_reference_chamber/c2_central_fans/videos/cycle.mp4)
or [operating timeseries](runs/03_reference_chamber/c2_central_fans/videos/operating_timeseries.mp4).
The owner selected this study on 12 September 2026 and revised its nominal size
on 14 September. The [previous 3.7 m version](archive/2026-09-14-before-size-revision/)
and [earlier studies](archive/2026-09-12-other-studies/) are preserved locally.

The documented nominal chamber is **6 × 4 × 4 m**, or **96 m³**, with an assumed
0.6 m elevated floor. The solver lists x × y × z as 4 × 6 × 4 m because x is
the sliding direction. Installed inside dimensions remain unverified. Three fans on a
fixed pole beside the trunk, at 0.667, 2.000 and 3.333 m above the floor, blow
in directions 120° apart. Cycles have **5 minutes
sealed and 10 minutes in the open phase**, including opening and closing travel.
CO₂ exchange acts continuously on representative leaflets along the fronds.

**Edit one file:** [configs/operating/current.yaml](configs/operating/current.yaml).
It exposes chamber and palm dimensions, wind speed/direction, fan count and
placement/directions/nominal flow, net exchange, cycle times and numerical settings.
[Every input and unit](docs/operating/INPUTS.md) ·
[Plain-language physics and equations](docs/operating/METHOD.md) ·
[Reproduction, outputs and model adapters](docs/operating/REPRODUCING.md).

```bash
python3.12 -m venv .venv-operating
.venv-operating/bin/python -m pip install -r requirements-operating.lock
# Install FFmpeg/ffprobe with libx264 on PATH; reference version 8.1.
.venv-operating/bin/python scripts/operating_chamber.py all \
  --input configs/operating/current.yaml \
  --output runs/03_reference_chamber/reproduction --workers 2
```

The output folder must be new. Each run freezes inputs, code, geometry and
environment, then writes fields, CSV histories, plots, four MP4/GIF views and an
audit. The familiar four-panel presentation uses white backgrounds and viridis.
The previous result is archived unchanged. The revised current result and new runs use generic filenames
(`cycle`, `fan_sections`, `operating_timeseries`, `leaf_source_location`).

```bash
.venv-operating/bin/python scripts/operating_chamber.py validate \
  --set chamber.width_m=5 --set air.wind_to_deg=90 --set fans.count=4 \
  --set exchange.net_co2_umol_s=-30
.venv-operating/bin/python -m pytest tests/operating tests/recreated/test_conservation.py
```

**Scientific status:** diagnostic, not validated. The default source is the
prescribed February 2025 Q95 net uptake, **−44.40047791098563 µmol/s**, not a 2026
measurement or dynamic leaf physiology. Fan performance and pole placement are
assumed; pole blockage is unresolved. The selected study's first-closure grid
screen fails at about 5.83% spatial variation against a 5% criterion. Conservation
passes; that does not establish field accuracy or opening convergence.

**Implementation:** the existing NumPy/SciPy conservative ALE solver in
`fluid_dynamic.recreated`, fitted geometry in `fluid_dynamic.scenarios`, and the
portable `fluid_dynamic.operating` workflow. The original local workspace also
retains legacy PhiFlow sources; they do not generate this selected study.

**Remote scope:** the dedicated operating branch contains source, a schematic
mesh, locked runtime dependencies, tests, documentation and explicit diagnostic
parameters. Results, archives, private observations, historical data-selection
files, CAD assets and virtual environments remain local. Result links above
therefore work in the original workspace or after a local reproduction, not in
an empty clone on GitHub. A fresh clone can generate the complete bundle without
any PalmTwin data release, old study or notebook.

[Current result and archive layout](docs/operating/RESULTS.md) ·
[Verification record](docs/operating/VERIFICATION.md).
