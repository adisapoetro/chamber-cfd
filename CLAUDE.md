# Operating-chamber CFD contract

The owner's 12 September 2026 instruction selects
`runs/03_reference_chamber/c2_central_fans` as the current operating-design
reference and archives every other existing result study. Preserve that bundle.
Use the physical folder layout; do not create a results HTML dashboard.

Read README.md, docs/operating/INPUTS.md, METHOD.md and REPRODUCING.md.
The active input is configs/operating/current.yaml and the launcher is
scripts/operating_chamber.py. Execute in a new output folder. The dedicated
operating Git branch is a portable source snapshot; existing local legacy code,
branches, private data-selection workflows and generated outputs are preserved.
Do not sweep them into a push. No force push, raw-data upload or broad git add.

The runtime is the existing NumPy/SciPy conservative moving-mesh solver, not
PhiFlow or WaterLily. Chamber dimensions are 4 × 6 × 3.7 m at default settings.
The prescribed February 2025 net source is not a current-year observation or
gross leaf photosynthesis. Fans add momentum only; uptake continues all phases
on representative leaflets. No pole blockage or dynamic leaf physiology is solved.
Retain failed sensitivity screens and the diagnostic/not-validated status.

Use Python 3.12 and requirements-operating.lock plus FFmpeg/ffprobe. The portable
checks are tests/operating and the included numerical tests in tests/recreated.
The original mixed workspace contains additional legacy tests requiring their
own inputs/dependencies. Changes to operators require their component tests and
numerical replay; presentation changes require media decode and visual review.
