# One active operating-design result

In the original PalmTwin workspace, open:

```text
physical-platform/whole-tree-chamber/simulation/fluid-dynamic/
  runs/
    README.md
    03_reference_chamber/
      README.md
      c2_central_fans/             CURRENT operating-design reference
        videos/                   Four MP4/GIF pairs
        figures/                  White-background figures
        data/                     Saved fields and histories
        configs/, geometry/, verification/
  archive/
    2026-09-12-other-studies/
      README.md
      runs/
        01_chamber_designs/        All earlier design comparisons
        02_wind_directions/       Earlier wind-direction study
        03_reference_chamber/
          c2_frond_uptake/         Earlier side-mounted fans
          fan_comparison/         Earlier fan/no-fan bundle
          c2_vpalm_leaves/         Earlier realistic-palm variant
          results/                Earlier reference bundle
          slides/                 Earlier explanatory slide/guide
    …                             Older archives unchanged
```

The original central-fan bundle retains its files and filenames. Start with
`videos/c2_202502_cycle.mp4`; `videos/scenarios_timeseries.mp4` is the operating
plot. Folder selection is an owner decision, not new field validation.

The dated archive preserves original relative subtrees. Original manifests and
frozen READMEs retain their source-era paths and statements. To locate moved
references, consult `run_locations.json` in the original workspace or this table;
the frozen result is not rewritten to pretend it was generated under a new path.
Do not replace archives or create compatibility symlinks in the active runs root.

The local maintenance record is
`simulation/maintenance/2026-09-12-operating-workflow/`: it records pre-change
hashes, the move mapping, file verification, numerical replay, tests and recovery.
Move a listed archive directory back to its recorded original path only when that
path is absent, then reverse the corresponding maintained path mapping. No result
is deleted by this housekeeping. Git publication contains code and rerunnable
inputs, not the local archive/video files or their private observation provenance.

Use a fresh output directory when running another diagnostic. After reviewing
it, explicitly choose whether it should replace the active reference; the
workflow does not silently promote a new run. Temporary verification reruns live
in the maintenance record, leaving the active runs folder uncluttered.
