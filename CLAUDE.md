# Repository conventions

Start with README.md and docs/operating/INPUTS.md, METHOD.md and REPRODUCING.md.
The entry point is scripts/operating_chamber.py; configs/operating/current.yaml
contains the reference inputs. Use Python 3.12, requirements-operating.lock
and FFmpeg/ffprobe.

Keep the numerical operators separate from configuration and presentation.
Changes to operators require the component tests and numerical comparison.
Changes to figures require media decoding and visual inspection. Run the
portable suite with `python -m pytest tests/operating tests/recreated`.

Generated runs belong in a fresh runs/ directory. The committed example under
examples/central_fans contains selected result files, inputs and provenance;
update its checksum manifest when replacing those files. Preserve existing
results and commit only the intended paths.

Keep units, coordinate conventions, source signs and numerical limitations
explicit. Nominal dimensions, schematic leaflets and assumed fan performance
must remain distinguishable from measurements. A passing test is not field
validation. Do not include private observations, field documents or credentials
in example outputs. Keep repository visibility and existing Git history intact.
