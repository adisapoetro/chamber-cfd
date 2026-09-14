# Operating-chamber verification

## Current size revision: 14 September 2026

The owner adopted documented nominal **6 × 4 × 4 m** geometry. The solver's
x/y/z dimensions are **4 × 6 × 4 m**, with 96 m³ rectangular airspace. Three
fan levels follow the equal-band rule: 0.667, 2.000 and 3.333 m above the floor.
The fixed support, wind, palm, prescribed uptake, cycle and numerical operators
retain their previous settings. Automatic exterior height increases to 9.3 m.

| Check | Current result |
|---|---|
| Main run | Complete: 1800 s, 0.5 m target cells, 0.5 s maximum steps |
| Numerical screens | Complete: first 300 s, 1/3 m grid and 0.25 s step |
| Saved fields and leaf-source maps | 243 of each across the three cases |
| Simulated closed-region volume | 96 m³ in all cases, independently checked from saved cell intersections |
| First-closure analytical mean | 396.605389 ppm; maximum main-run error 1.46e-07 ppm |
| Main CO₂ inventory residual | Maximum fraction 2.41e-08 |
| Main incompressibility residual | Maximum 7.09e-10 s⁻¹ |
| Moving-volume residual | Main maximum 1.22e-12 m³ |
| Leaflet allocation and saved-field statistics | PASS; uptake is continuous and excluded from trunk/bare frond axes |
| Grid sensitivity | **FAIL**: 5.834% spatial-standard-deviation change against 5%; mean speed changes 4.245% |
| Time-step sensitivity | PASS: 1.124% spatial-standard-deviation change; mean speed changes 0.541% |
| Component tests | 137 passed, 1 skipped, 400.13 s; skipped historical-fixture parity test is separate from this complete run |
| Media | All four MP4/GIF pairs decode; white backgrounds and unclipped viridis scales pass |
| Visual review | Cycle closure/opening, fan section, operation plot, source location and updated fan-height schematic inspected |
| Palm/source preservation | All seven geometry arrays exactly equal the previous result |
| Previous result preservation | All 592 files retain their pre-move SHA-256 hashes |

Detailed current records are `verification/checks.json` and
`verification/bundle_manifest.json` inside the result. The local revision and
recovery record is `simulation/maintenance/2026-09-14-operating-size-revision/`.

## Historical 3.7 m result

The 13 September portable-workflow verification concerned the previous
**88.8 m³** geometry. Its numerical replay was exactly equal across 243 saved
fields; the grid screen failed at 6.867% spatial variation, and the time-step
screen passed at 0.559%. The earlier source-only checkout passed 79 tests.
These historical values do not describe the revised 96 m³ run. The previous
bundle is now in `archive/2026-09-14-before-size-revision/c2_central_fans`;
its inputs remain recoverable from source commit
`39ad6840cfe712cf0e986116701f9f81ba547615` and the earlier maintenance record.

## Scientific, privacy and recovery limits

These are implementation and numerical checks, not empirical validation.
Installed clear dimensions, fan performance and pole position remain assumed;
pole blockage is unresolved, leaflets are schematic and physiology is prescribed.
The −44.40047791098563 µmol/s source retains its February 2025 Q95 provenance;
it is not a new 2026 observation. No observation reprocessing or rate rescaling
was performed when changing the chamber height. The grid screen still fails,
and opening-grid convergence remains unestablished.

Only the existing portable source allowlist is eligible for the authorized
private Git update. Raw observations, holdouts, CAD, private data-selection
workflows, outputs, archives and environments remain local. Existing unrelated
working-branch changes remain preserved. A source commit is not a remote backup
of the result videos or archived evidence.
