# Reading the results

The [central-fan example](../../examples/central_fans/) contains ready-to-view
results. A full run produces this structure in the requested output folder:

```text
README.md                   Links to the run's results and checks
METHOD.md, INPUTS.md,
REPRODUCING.md               Documentation captured with the run
configs/                    Resolved inputs, source snapshot and environment
geometry/                   Schematic palm and representative leaflet surfaces
videos/                     Four views, each as MP4 and GIF
figures/                    Still frames, histories, source and fan diagrams
data/cfd_runs/              Main, finer-grid and smaller-time-step cases
verification/               Integrity, conservation, sensitivity and media checks
logs/                       Calculation log for each case
```

## Videos and figures

| File stem | What it shows |
|---|---|
| `cycle` | 3D CO₂ slices, airflow, a horizontal section and concentration histories |
| `fan_sections` | Airflow and CO₂ in a vertical section through the fan support |
| `operating_timeseries` | Mean concentration, spatial variation, net exchange and opening/fan states |
| `leaf_source_location` | The representative leaf surfaces and where exchange enters the calculation |

The timestamp on each frame is physical simulation time. Playback slows around
wall movement, so video duration is not simulated duration. Figure names such
as `cycle_0300s.png` identify the corresponding physical time.

The CO₂ colour scale is fixed within a video and expands if needed to include
all calculated concentrations. Purple indicates lower concentrations and
yellow indicates higher concentrations. The source-location view has a
separate scale for depth-integrated removal; it is not a CO₂ concentration map.

## Histories and saved fields

Each case's `timeseries.csv` records time, opening gap, chamber-region mean and
standard deviation, airflow statistics, source rate, fan state and numerical
balance diagnostics. The reported region remains the original closed-chamber
volume even while the two shells move apart.

The mean is volume weighted. Spatial standard deviation measures differences
between air cells; it is not a confidence interval or sensor uncertainty.
The net-exchange rate is in µmol/s for the whole plant: negative means removal
from air. Source allocation follows representative leaflet area within each cell.

A full run also stores compressed NumPy snapshots with CO₂, velocity, cell
coordinates and volumes, alongside leaflet-area maps. `summary.json` lists
saved times and paths. The small committed example provides histories and
plots; rerun its input to generate all spatial arrays.

## Read the checks before interpreting the flow

`verification/checks.json` separates integrity/conservation checks from
numerical sensitivity. A run can complete successfully while a sensitivity
check fails. In the reference case, the first-closure grid check exceeds the
5% criterion; see [VERIFICATION.md](VERIFICATION.md).

`verification/bundle_manifest.json` contains the full run's checksums. The
committed example has its own `manifest.json`, covering only the files included
in the repository. Check those with `python scripts/check_example.py`.
