# One active operating-design result

Current result, revised 14 September 2026: **nominal 6 × 4 × 4 m**, 96 m³.
In the original PalmTwin workspace, open:

```text
fluid-dynamic/
  runs/
    README.md
    03_reference_chamber/
      README.md
      c2_central_fans/             CURRENT: nominal 4 m air height
        videos/                   Four MP4/GIF pairs
        figures/                  White-background figures
        data/                     Saved fields and histories
        configs/, geometry/, verification/
  archive/
    2026-09-14-before-size-revision/
      README.md
      c2_central_fans/             Previous 3.7 m height; original bytes
    2026-09-12-other-studies/
      README.md
      runs/
        01_chamber_designs/
        02_wind_directions/
        03_reference_chamber/     Earlier frond/fan/palm/reference variants
    …                             Older archives unchanged
```

The current bundle uses `videos/cycle.mp4`, `fan_sections.mp4`,
`operating_timeseries.mp4` and `leaf_source_location.mp4`, each with a GIF.
These are the four views produced by the portable workflow. The archived
3.7 m version keeps `c2_202502_cycle` and `scenarios_timeseries` filenames.

The revised geometry follows the owner's adoption of recent documentary nominal
dimensions. It is not an inside-to-inside survey or field validation. Read
INPUTS.md for the dimension/axis convention and VERIFICATION.md for the checks.

Original manifests and frozen READMEs retain their source-era paths and claims.
The current canonical folder now denotes the revised run; the old bundle's
historical path is recorded separately in `run_locations.json` under
`operating_size_revision`. Do not redirect the active path to the archive.

Local recovery records are in
`simulation/maintenance/2026-09-14-operating-size-revision/` and, for the prior
workflow publication and broad housekeeping,
`simulation/maintenance/2026-09-12-operating-workflow/`. They contain hashes,
before-state files, checks and move records. No result was deleted. Restore an
archived bundle only to an absent destination after preserving the current one.

Git publication contains rerunnable code and inputs, not videos, raw observations
or local archives. Use a fresh folder for each rerun. A newly executed run does
not automatically replace the current reference; the 14 September replacement
was explicitly requested by the owner. Keep temporary verification runs in the
maintenance directory and retain the physical folder layout without an HTML hub.
