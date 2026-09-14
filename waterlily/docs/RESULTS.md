# Operating chamber: generated diagnostic

The complete 1,800 s WaterLily + Julia CO₂ run contains 181 saved times, covering
two cycles of 300 s sealed and 600 s unsealed. The fixed inside volume is 96 m³.
Prescribed net exchange is −44.40047791098563 µmol/s on representative leaflets,
including while open. Three central fans operate during closure.

The [repository example](../examples/current_chamber/README.md) contains all four
MP4/GIF views and 24 PNG stills. Full arrays and the frozen generating source
are in the local `results/current_chamber/` folder or can be regenerated using
the project README. Short checks are grouped under `results/diagnostics/`.

| Check | Result |
|---|---|
| Julia package tests | 320 assertions pass |
| Completed physical duration | 1,800 s |
| Saved concentration/velocity snapshots | 181 |
| Maximum geometric-law residual | 4.51 × 10⁻¹³ m³/s |
| Maximum cumulative CO₂ budget error | 6.82 × 10⁻⁹ ppm m³ |
| Lowest concentration across all solver steps | 390.159 ppm |
| Largest chamber-region scalar flux correction | 163.81% — **screen fails** |
| Short-run spatial SD change on grid refinement | +34.30% — **not grid-converged** |

At the end of the first closure, mean CO₂ is 396.6054 ppm: the expected 3.3946 ppm
decline from the initial 400 ppm. The second closure starts below ambient after
incomplete recovery and loses the same mean amount under the constant source.
These are numerical outputs, not observations.

The closing-gap correction is large. Its maximum velocity change over the
time-averaged sharp aperture reaches about 134.6 m/s in a vanishing aperture.
That is a scalar-coupling diagnostic, **not a credible predicted chamber air
speed**. It exposes disagreement between coarse native BDIM flow and sharp
transport geometry as the thick shells meet. The video arrows show native
WaterLily velocity; the corrected scalar fluxes are separately audited.

Mass conservation and positivity do not remove this limitation. The coarse
example is suitable for inspecting the workflow and operating sequence, but
its mixing rates and gradients should not be used as quantitative evidence.
[The verification notes](VERIFICATION.md) describe the fixed-wall refinement
check and remaining wall/interface qualification.
