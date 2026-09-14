# WaterLily chamber project

This folder is an independent Julia environment inside chamber-cfd. Its solver
is native WaterLily 1.8.0 with a separate conservative CO₂ extension. The parent
repository's Python operating workflow remains a separate reference.

- Run `julia --project=. --threads=1 -e 'using Pkg; Pkg.test()'` after code changes.
- Preserve native momentum and pressure. Audit scalar interface mapping,
  geometric conservation, closed-component compatibility and CO₂ inventory.
- Keep exact input units, chamber/plant dimensions, rate provenance and cycle
  timing in the TOML configuration. Net exchange acts on leaflets in all phases.
- Generate into a fresh result directory. Preserve frozen source, package pins,
  audits and SHA-256 manifests. Rendering reads saved fields.
- Keep failed numerical screens visible. A successful run or test is not
  experimental validation, grid convergence or calibrated photosynthesis.
- Keep the white backgrounds and fixed viridis concentration scale. No HTML
  results dashboard is needed.
- Version source and the compact `examples/current_chamber` bundle. Generated
  full fields and local depots stay outside Git. Do not include raw research data
  or unrelated parent-repository changes.
