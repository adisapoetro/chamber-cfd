# Numerical checks

The [central-fan example](../../examples/central_fans/) completes two operating
cycles and passes conservation and file-integrity checks. Its first-closure
grid-sensitivity check exceeds the specified 5% criterion. These checks assess
the implementation and numerical behaviour; they do not validate the model
against measurements.

## Reference calculation

The chamber has a nominal rectangular air volume of **96 m³**: 4 × 6 × 4 m
in solver coordinates. The main run covers 1800 s with 0.5 m target cells and
0.5 s maximum time steps. Additional cases cover the first sealed 300 s using
a finer grid or a smaller time step.

| Check | Main-run result |
|---|---|
| Enclosed volume, independently calculated from saved cell intersections | 96 m³ |
| Mean CO₂ after the first 300 s closure | 396.605389 ppm |
| Maximum error against the first-closure analytical mean | 1.46 × 10⁻⁷ ppm |
| Maximum normalized CO₂ inventory residual | 2.41 × 10⁻⁸ |
| Maximum incompressibility residual | 7.09 × 10⁻¹⁰ s⁻¹ |
| Maximum moving-volume residual | 1.22 × 10⁻¹² m³ |
| Maximum difference between saved-field and recorded means | 1.14 × 10⁻¹³ ppm |
| Maximum leaflet-source allocation error | 7.11 × 10⁻¹⁵ µmol/s |

The analytical mean check follows directly from the prescribed net exchange,
temperature, pressure and sealed volume. It is independent of how fans
redistribute CO₂ within that volume. The finer-grid and smaller-time-step
cases also pass the conservation and source-allocation checks.

## Sensitivity to grid and time step

These comparisons use the relative RMS difference from the main run over
the **first 300 s only**. Each monitored quantity must stay below 5% to pass.

| Case | Target cell size | Maximum step | Spatial CO₂ standard deviation | Mean air speed | Result |
|---|---:|---:|---:|---:|---|
| Finer grid | 1/3 m | 0.5 s | 5.834% | 4.245% | **FAIL** |
| Smaller time step | 0.5 m | 0.25 s | 1.124% | 0.541% | PASS |

The grid result means spatial variation is still sensitive to resolution at
the chosen threshold. Passing the time-step comparison alone does not establish
convergence. Neither comparison tests the opening phase. A quantitative study
of mixing time or fan placement needs further resolution studies and comparison
with measurements.

## Reproduction and media

The examples were run on 14 September 2026 using Python 3.12.13, the pinned
dependencies and FFmpeg 8.1. Both documented commands completed: the four-second
installation example and the full reference calculation.

The standalone repository's **79 tests pass**. A repeat of the reference
calculation exactly reproduced the three case histories, all **243 saved
CO₂/velocity fields**, their **243 leaflet-source maps**, and the seven palm
geometry arrays from source revision `cf7728660855923457626e6deaaa8f43824bf151`.
This equality was checked in the same pinned environment; other platforms
may introduce floating-point or encoding differences.

The published example contains **4 MP4 files, 4 GIF files and 24 PNG figures**.
All decode successfully. The four views were visually inspected for readable
labels, white backgrounds, fan positions and the opening/closing sequence.

Exact numerical results are in
[checks.json](../../examples/central_fans/checks.json). The
[provenance record](../../examples/central_fans/provenance.json) records the
environment, source hashes and replay comparison. Verify the published file
list, checksums and media from the repository root:

```bash
python scripts/check_example.py --decode
```

## Physical assumptions

Installed inside dimensions, fan performance and support position still need
measurement. The model uses schematic leaflets, omits support and fan-housing
blockage, and prescribes whole-tree net exchange independently of the local
environment. The historical February 2025 rate is a diagnostic input, not a
current-year measurement. See [the model formulation](METHOD.md) for the
equations and remaining physical limits.
