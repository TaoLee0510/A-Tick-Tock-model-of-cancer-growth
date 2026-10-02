# Portable output directories

All shipped YAML configurations use relative, distinct output directories.
They are relative to the process working directory. ABM configurations use
`output.directory`; continuum configurations use the same key. Nutrient
wrappers may override the base model's directory with
`nutrient.output.directory`. Structured-PDE wrappers may override their
continuum directory with `output.directory`. These overrides affect output
only and do not change checkpoint dynamics identity.

The four executables (`atcg3d`, `atcg3d_nutrient`, `atcg3d_continuum`, and
`atcg3d_structured_pde`) accept `--output-root PATH`. A relative YAML directory
is appended to that root. An empty root is rejected; omitting the option
preserves the configured directory. An explicitly absolute directory in a
local YAML file retains its meaning. Keep machine-specific paths out of
committed configurations.

Set `EXTERNAL_RUNS` locally to the desired directory on an external volume,
then run, for example:

```sh
build-codex/atcg3d --config configs/atcg3d_smoke_test_v3.yaml \
  --output-root "${EXTERNAL_RUNS:?set EXTERNAL_RUNS}" --dry-run
build-codex/atcg3d_structured_pde \
  --config ATCG3D_StructuredPDE/config/structured_smoke_r20_v1.yaml \
  --output-root "${EXTERNAL_RUNS:?set EXTERNAL_RUNS}"
```

`--dry-run` shows the effective output path without creating output directories.
Run metadata records the effective path; relocation does not affect biological
or numerical fingerprints. Existing run-directory overwrite guards still
apply.

Input paths (`base_config`, `continuum_config`, and checkpoint paths) retain
their existing resolution relative to the referencing YAML file. The output
root does not rewrite them. For continuation from an external volume, use a
local resume YAML with the correct checkpoint path. The shipped resume
examples point to their source run's relative directory and use a distinct
output directory for the resumed segment.

## Nutrient metrics continuation

New `nutrient/metrics.csv` files have 15 columns, including assembled r and K
consumption rates. Resuming a run with the legacy 13-column header preserves
that header and omits the two new columns from subsequent rows. A recognized
current header is appended with all 15 columns. CRLF headers are accepted.
Unknown or unreadable headers are rejected before the file is opened for
append, so existing data is preserved.

## Legacy PNG rendering

The optional historical 2D executable reads its font from `ATCG_FONT_FILE`.
If unset or empty, it looks for `Calisto MT.ttf` in the working directory.
Supply a locally licensed font file; no user-specific font path or font binary
is stored in the repository. This affects PNG annotation only.
