# ATCG3D population continuum model

`ATCG3D_Continuum` is a runnable population-level reduction of the ATCG3D
individual-cell model. It evolves four cell-number density fields

```text
r_small, r_large, K_small, K_large
```

together with the effective-nutrient field. It is intended for comparing
spatial population distributions under crowding and vascular nutrient supply;
it is not a cell-by-cell replacement for the PDMP model.

## Build and run

```sh
cmake -S . -B build-continuum -DCMAKE_BUILD_TYPE=Release \
  -DATCG_BUILD_LEGACY_2D=OFF \
  -DATCG_BUILD_3D=ON \
  -DATCG_BUILD_3D_CONTINUUM=ON \
  -DATCG3D_ENABLE_HDF5_CHECKPOINT=ON
cmake --build build-continuum --parallel

build-continuum/atcg3d_continuum \
  --config build-continuum/configs/continuum_smoke_v1.yaml \
  --dry-run

build-continuum/atcg3d_continuum \
  --config build-continuum/configs/continuum_smoke_v1.yaml
```

The smoke profile uses the exact 64-cell ABM initialization and an off-centre
synthetic vessel so nutrient-driven spatial heterogeneity is exercised in a
short run. `continuum_legacy_v1.yaml` uses the production mapped biology on a
coarser three-dimensional continuum grid.

`continuum_legacy_2d_2000_v1.yaml` is the native legacy-plane profile. It uses
a `2000 x 2000 x 1` unit-spaced lattice, runs to 2160 hours, and represents the
synthetic vessel as a y-directed line at `x = -40`. Its r population uses a
20-fold mobility multiplier on faces at or above 0.9 occupied fraction; the
exclusion mobility on those faces remains active. The corresponding base
configuration is `configs/atcg2d_legacy_native_v3.yaml`; its unscaled 2D density
limits and carrying capacities are used by the reaction term.

## ABM-derived initialization

`initialization.mode: base_model` constructs the exact deterministic initial
ATCG3D cell state from the referenced schema-v3 base configuration and then
coarse-grains it conservatively. Each large footprint is distributed across
its occupied sites before binning, so integrated cell number and occupied
volume are both preserved.

To begin from an evolved ABM state, use a matching HDF5 checkpoint:

```yaml
initialization:
  mode: abm_checkpoint
  abm_checkpoint: /absolute/path/checkpoint_....h5
vascular:
  source_mode: abm_perfusion
```

The checkpoint's cells and perfused vascular raster are imported directly.
HDF5 support must be enabled at build time. The configured continuum domain
must contain every imported cell footprint. `abm_plus_synthetic_line` can be
used for a controlled vascular perturbation on top of imported vessels.

Continuum schema v2 changes the nutrient sink from occupied-volume demand to
cell-number demand. A large and small cell of the same type therefore consume
the same integrated amount, while the supplied v2 structured profile sets r
demand to 1.2 times K. Occupancy, exclusion and carrying-capacity calculations
remain volume-weighted. The simple continuum model does not carry persistent
direction states; use `ATCG3D_StructuredPDE` for the aligned directed-migration
rule.

Continuum schema v3 is the shared-nutrient closure used by structured-PDE v5.
It evolves nutrient transiently, starts from an explicit configured value, and
holds selected planar edges and/or vessel voxels at maximum supply. Every cell
has the same per-cell demand, r and K use common density limits/carrying
capacity, and nutrient scales positive growth instead of changing hard spatial
capacity. Older schemas keep their previous equations.

## Output

The output directory contains:

```text
metrics.csv                         population totals and spatial summaries
fields/field_XXXXXXXX.csv           four densities, occupancy, nutrient, vessel
profiles/profile_XXXXXXXX.csv       radial population/nutrient profiles
checkpoints/*.continuum.bin         exact continuum restart state
config/requested.yaml               supplied wrapper
config/effective.json               resolved configuration identity
final.json                          final summary and checksum
```

`metrics.csv` reports r/K totals, size-resolved totals, occupied volume,
maximum occupied fraction, r/K mean nutrient exposure, r/K mean radius, vessel
volume, and a state checksum. Field snapshots contain only nonzero rows. The
optional `output.field_stride` integer samples every Nth grid coordinate in
each dimension when writing CSV fields; it reduces output size without changing
the simulated grid or dynamics.

To resume, keep the same dynamics and set:

```yaml
run:
  mode: resume
  resume_checkpoint: /absolute/path/checkpoint_....continuum.bin
```

The binary checkpoint stores all five evolving fields, the vascular source,
clock, nutrient refresh schedule, and a dynamics fingerprint. A different end
time is allowed; a dynamics mismatch is rejected.

## Interpretation boundary

This deterministic PDE model predicts coarse-grained mean density. It does not
retain cell UID, lineage, finite-number fluctuations, exact random clocks, or
individual hard-conflict outcomes. Those remain outputs of `atcg3d_nutrient`.
The two models can share initial cells, vascular state, density-growth law, and
nutrient parameters, which makes ensemble and continuum comparisons explicit.

See [`docs/model_spec.md`](docs/model_spec.md) for equations and numerical
details.
