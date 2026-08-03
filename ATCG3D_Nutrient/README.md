# ATCG3D Nutrient-Coupled Model

`ATCG3D_Nutrient` contains the runnable nutrient-coupled ATCG3D model. It
preserves individual cells, hard voxel occupancy, event clocks, the fixed
26-direction migration rules, and the discrete vascular graph from `ATCG3D`,
while adding a continuous effective-nutrient field.

The intended mathematical model is a hybrid system:

```text
cell and vessel PDMP + effective-nutrient reaction-diffusion PDE
```

The production `ATCG3D/` source tree remains the reference implementation. It
contains only the optional environment interface needed by the common event
engine; the baseline `atcg3d` executable creates no environment and retains its
existing behavior. Nutrient configuration, solver, output, and application code
remain isolated here.

## Build and run

```sh
cmake -S . -B build-nutrient -DCMAKE_BUILD_TYPE=Release \
  -DATCG_BUILD_LEGACY_2D=OFF \
  -DATCG_BUILD_3D=ON \
  -DATCG_BUILD_3D_NUTRIENT=ON \
  -DATCG3D_ENABLE_HDF5_CHECKPOINT=ON
cmake --build build-nutrient --parallel

build-nutrient/atcg3d_nutrient \
  --config build-nutrient/configs/nutrient_smoke_v1.yaml \
  --dry-run

build-nutrient/atcg3d_nutrient \
  --config build-nutrient/configs/nutrient_smoke_v1.yaml
```

Use `nutrient_legacy_v1.yaml` for the full production base profile. The smoke
profile deliberately uses a shorter halo, fewer fixed solver iterations, and a
four-hour refresh interval so end-to-end validation remains fast.

Wrapper paths are resolved relative to the wrapper file. For a copied or
mirrored bundle, if that relative path is unavailable, the loader also accepts
the referenced base configuration beside the wrapper; CMake copies both into
`build-nutrient/configs/`.

## Compatibility target

The implementation replaces the prescribed distance-only vascular
density relief with a supplied, diffused, and cell-consumed nutrient field. It
should keep the current biological scope:

- nutrient affects density-dependent growth and death scheduling;
- cell occupancy and migration feasibility remain discrete and unchanged;
- vessel nodes, segments, tips, collisions, and displacement remain discrete;
- the no-nutrient-extension compatibility limit reproduces the reference
  ATCG3D transition rules;
- a fixed seed, deterministic field solver, and restored checkpoint reproduce
  the same trajectory within the new model.

See [`docs/model_spec.md`](docs/model_spec.md) for the equations,
coupling contract, and implementation boundaries.

## Source layout

```text
ATCG3D_Nutrient/
  app/          `atcg3d_nutrient` entry point
  config/       strict schema-v1 wrapper and supplied profiles
  field/        sparse deterministic nutrient grid and PDE solver
  io/           metrics, field snapshots, and checkpoint sidecars
  docs/         model and numerical specification
```

## Checkpoint and resume

The standard HDF5 checkpoint remains authoritative for cells, event clocks,
lesions, and the vascular graph. Every nutrient run checkpoint has a required
same-stem sidecar:

```text
checkpoint_...h5
checkpoint_...nutrient.bin
```

The sidecar stores the exact sparse field, refresh schedule, configuration
fingerprint, base-state checksum, and clock identity. Resume refuses a missing,
mismatched, or truncated sidecar. This keeps existing HDF5 schemas readable by
the baseline model while giving the coupled model exact continuation.

## Implemented numerical contract

- quasi-steady reaction-diffusion equation;
- fixed-count relaxed Jacobi iterations for deterministic results;
- sparse blocks allocated around cells and perfused vessel voxels plus a halo;
- Michaelis-Menten r/K consumption assembled from individual footprints;
- vessel exchange assembled from the existing perfused capsule voxels;
- capacity multiplier in `[1, M_max]`, with `M_max=2` reproducing the old
  maximum relief of `0.5`;
- deterministic `environment_refresh` events that integrate existing division
  work before rescheduling growth and death;
- optional nutrient metrics and nonzero-voxel CSV snapshots.

The v1 model is a compatibility model: nutrition represents vascular support
above the unresolved avascular baseline. It does not yet make below-baseline
starvation a separate phenotype.
