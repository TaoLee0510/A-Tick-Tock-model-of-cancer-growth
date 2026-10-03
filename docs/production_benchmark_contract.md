# Production benchmark protocol

The three integrated inputs use a 2000 by 2000 thin grid, r200, transient
nutrient, shared VEGF lattice angiogenesis and a configured horizon of
2160 hours. They select the new prepared sector means and sparse population
storage. The complete nutrient and vascular fields remain dense. These are
production workload definitions; their existence does not establish long-time
model equivalence, acceptable memory at every tumor size, or independence
from the finite boundary.

`benchmark_production.py` records independent native runs of ABM, structured
PDE and adaptive hybrid. For each model it runs a reference prefix, writes a
checkpoint at half that prefix, and restores in a separate process with eight
threads. The configured biological horizon remains 2160 hours in all runs.
The stopping time must be a resource macro boundary; ABM stops after its
resource event at that boundary. Exact state and resource checksums, mass,
active mass and vascular diagnostics must match after restore. No numerical
allowance applies to this comparison. The ABM checkpoint includes HDF5 cell
state and its resource sidecar; hybrid includes all three native files.

Timing includes construction, resource-limited initialization or restore,
all operators, final diagnostics, optional checkpoint writing and teardown.
Peak resident memory includes derived sector caches and concurrent seeding
state. Checkpoint sizes include every sidecar. Reference and resumed timings
have different work scopes and are reported separately, without a speedup
claim. Record operating system, architecture, thread counts and raw values.

The initial opt-in measurement is a 1-hour prefix with a 0.5-hour checkpoint.
Run the same command with `--stop-hours 2160` for the full workload. Every
report includes both the actual stopping time and configured horizon. A short
prefix cannot validate the full duration. The unchanged interior screening
criterion is r99 below 0.8 of the planar half width at each observation. A
radius uses all r and K cell mass, voxel centers (agent anchor plus 0.5), and
unit-width radial bins with the lower bin edge as the reported quantile.
This common convention is used for all three production models. A
failure is reported as boundary contact; it must not be made a passing
equivalence result. This benchmark does not replace the multi-seed TOST
validation or the opt-in long ensembles.

The default CTest checks the restart path on a smaller 256-square fixture
with r20 and dynamic vasculature. It retains the 2160-hour configured horizon,
runs 0.5 hours, and restores a 0.25-hour checkpoint. Large measurements have
the opt-in `benchmark_production` label and require the HDF5 build for ABM.

```sh
cmake -S . -B build-production -DCMAKE_BUILD_TYPE=Release \
  -DATCG_BUILD_LEGACY_2D=OFF -DATCG3D_ENABLE_HDF5_CHECKPOINT=ON \
  -DATCG3D_ENABLE_LARGE_TESTS=ON
cmake --build build-production -j
ctest --test-dir build-production -L benchmark_production --output-on-failure
python3 scripts/benchmark_production.py \
  --exe build-production/atcg3d_production_benchmark \
  --output build-production/production-full.json --stop-hours 2160
```

The production duration grid retains a 0.5-hour bin width and expands its
maximum to 2048 hours. A 32-hour span rejected resource-limited initial active
clocks in the smaller fixture (maximum 344.588 hours). This is a clock-domain
choice, not a statistical allowance. Import still rejects clocks exceeding
the explicit span. Initial geometric work and resource-dependent division
times are stochastic, so other seeds or later states can require a larger
span or a differently resolved clock model. Long-time duration and spatial
closure accuracy remain subject to a separate refinement/equivalence study.
