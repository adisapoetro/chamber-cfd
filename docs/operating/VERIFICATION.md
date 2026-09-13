# Verification of the portable operating workflow

Completed 13 September 2026; work and recovery record began 12 September.
These are implementation and numerical checks, not empirical validation.

| Check | Result |
|---|---|
| Default full numerical replay | All histories and all 243 saved fields exactly equal the preserved central-fan result; maximum absolute difference 0 |
| Palm/source geometry replay | Triangle coordinates, masks, display mesh and frond axes exactly equal |
| Main calculation | 1800 s, 0.5 m target cells, 0.5 s maximum steps |
| Numerical screens | 300 s grid refinement and 300 s smaller-step runs reproduced |
| Source and inventory checks | PASS: leaf area, total exchange, first sealed analytic balance, moving volume, incompressibility and saved-field statistics |
| Grid sensitivity | **FAIL**: 6.867% spatial standard-deviation change against 5%; mean speed changes 4.053% |
| Time-step sensitivity | PASS: 0.559% spatial standard-deviation change; mean speed changes 0.550% |
| Original-workspace component tests | 138 passed, including independent transport/flow checks and changed-input tests |
| Clean source-only checkout | 79 included tests passed with only the locked public dependencies |
| Full default rendering | Four MP4/GIF pairs; full decoding, colour-range and white-corner checks pass |
| Changed-input complete bundles | Zero-fan/zero-source control and independent one-fan/positive-release run pass, including all four media views |
| Physical archive | 24,120 pre-change files preserved; two replaced navigation guides retained in recovery |
| Path lookup tests | 3 passed |

The independent fresh-environment run changed chamber width/height, wind speed
and direction, fan count and azimuth, exchange sign and operating timing. Separate
integration tests also exercise four fans, uptake, release and zero exchange.
No older result, private data table, notebook or installed PalmTwin package was
available in that source-only checkout. It uses Python 3.12.13 and the committed
runtime pins; the reference media encoder is FFmpeg 8.1.

The selected `c2_central_fans` result remains byte-identical. Temporary reruns,
logs, complete field comparisons, media hashes, before-state records and archive
move checks are local under
`simulation/maintenance/2026-09-12-operating-workflow/`.
The default reproduction's `verification/checks.json` and
`verification/bundle_manifest.json` contain its detailed checks and hashes.

Privacy: no observation rows, holdouts, private data-selection files, Blender/CAD
assets, local paths to raw observations, generated runs or archives are included
in the publication allowlist. The single explicit prescribed rate is an authorized
diagnostic input, not a new observational data release. The publication uses a
separate root branch containing only the reviewed workflow and dependencies; the
existing local branch/history and unrelated dirty work remain preserved.

Limits remain: assumed installed fan performance and pole position, unresolved
pole blockage, schematic leaflets, conditional source/geometry identity, no
dynamic leaf physiology, and no established opening-grid convergence. A code
replay or successful media audit cannot remove the failed grid screen.
