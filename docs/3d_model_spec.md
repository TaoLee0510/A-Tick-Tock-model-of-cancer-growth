# ATCG three-dimensional model specification

## Scope and compatibility

The historical `ATCG/` implementation and `atcg` target remain the legacy 2D
model. The 3D implementation is isolated under `ATCG3D/` and builds as
`atcg3d`. PNGWriter, PNG, Freetype, Blitz, and the legacy dense `VisualRange`
are not dependencies of the 3D target.

The one shared biological primitive is
`ATCG/common/density_growth_rule.{hpp,cpp}`. The legacy header
`ATCG/density_growth_rate_calculation_1.hpp` continues to expose the same
`calculate_density_growth_rate()` API, while 3D calls the same pure function
after obtaining unique-anchor counts from its block index.

## Fixed direction set

Direction 0 is stay. Directions 1–8 are the original clockwise z=0 layer:

| ID | vector | ID | vector |
|---:|:---:|---:|:---:|
| 1 | (-1,-1,0) | 5 | (1,1,0) |
| 2 | (-1,0,0) | 6 | (1,0,0) |
| 3 | (-1,1,0) | 7 | (1,-1,0) |
| 4 | (0,1,0) | 8 | (0,-1,0) |

Directions 9–17 use z=-1 in the same planar order plus `(0,0,-1)`;
18–26 use z=+1 plus `(0,0,+1)`. `DirectionSet26` is a constexpr table.
Movement, angle, opposite-direction, cone, and footprint operations are generic
vector algorithms; there are no 26 copied branches.

## Geometry and stages

- Stage 0 (large): anchor plus `{0,1}^3`, exactly eight voxels.
- Stage 1 (small): one anchor voxel.
- Stage 2 (ultrasmall): one anchor plus an explicit co-location group.
- A large-cell move tests only `translated footprint - current footprint`.
  Axis, planar-diagonal, and spatial-diagonal moves enter 4, 6, and 7 voxels.
- Stage-1 recovery enumerates the eight 2×2×2 anchors containing its current
  site and chooses uniformly from feasible candidates.
- Removing from a co-location group updates the remaining group to stage 2 for
  two or more cells, stage 1 for one, and releases the voxel for zero.

`SparseChunkGrid3D` allocates 32³ chunks on demand by default. Each allocated
voxel stores one `uint32` primary slot. Only genuinely co-located voxels own an
overflow vector. Chunk lookup is hash-based; individual voxels are not map
entries. Floor division handles negative coordinates.

Thin-layer mode restricts cell anchors and migration/division directions to
z=0. The occupancy policy permits z=1 only for the upper half of a large
2×2×2 footprint. This makes the direction rule reduce exactly to the expected
2D 90/5/5, 90/10, and 50/50 cases without changing the 3D stage definition.

## Migration

K cells and random migration choose uniformly over feasible lattice directions
when `direction.distance_weight_exponent=0`. A future distance correction uses
weight `|d|^-exponent` without changing direction IDs.

For an r cell with no persistent direction, each feasible direction is filtered
by a forward cone satisfying Chebyshev distance 1 through `density_radius` and
angle at most `density_half_angle_deg`. Directions with estimated density at or
below the threshold remain eligible. The production estimator incrementally
counts unique biological anchors per block. Exact whole-store estimators exist
only for small tests.

For persistent direction `d0`, eligible turns are feasible directions other
than `d0` within `turn_half_angle_deg`. If `d0` and turns exist, continuation
has the configured probability and the remainder is uniform over turns. If
only `d0` exists it continues; if only turns exist they are uniform; otherwise
the cell stays and resets persistence. Density filtering is not repeated during
persistence unless `direction.persistence_uses_density=true`.

Defaults are radius 5, density half-angle 45°, threshold 0.60, turn half-angle
45°, continuation 0.90, and distance exponent 0.

## Division and death

Large division enumerates the radius-2 Chebyshev shell around the mother
anchor, exactly `5^3-3^3=98` daughter anchors, and retains complete available
footprints. If none exists and shape reduction is enabled, the mother footprint
is cleared and two distinct sites are selected from the 64 sites in
`[-1,2]^3`; mother and daughter become stage 1.

Stage-1 division selects one of 26 feasible neighbor anchors. When blocked, an
r mother follows the legacy non-survival rule. A blocked K mother creates an
explicit stage-2 co-location when ultrasmall behavior is enabled. Stage-2
separation moves one cell at a time and then derives the remaining stage from
group size.

Density-dependent growth uses a 6×6×6 anchor window by default. The initial 3D
capacity mapping multiplies the 2D limits/capacities by six. All values are
configuration fields and require scientific calibration; the factor six is a
starting interpretation, not a validated biological constant.

Initialization maps the disk and annulus to a sphere and spherical shell. It
generates actual sphere/shell candidate anchors by z/y slices rather than
scanning the global sparse domain or a dense X×Y×Z array. Footprint availability
is checked before every placement.

## Storage, identity, and memory

`CellStore3D` is a typed structure-of-arrays store with stable `uint32` slots
and a free list. UID and parent UID are immutable 64-bit biological identities;
UID never equals slot by contract and never changes after slot reuse. Position
anchors are three `int32` columns. Type, stage, viability, flags, and last
direction are byte columns. Rates are floats; event times remain doubles.
Event sequence and schedule generation are explicit mutable state.

The logical column width is 90 bytes per allocated slot. The measured synthetic
10-million-cell core, including CellStore, sparse occupancy, and density index,
was 99.74 bytes per cell on the tested macOS arm64 machine.

## Event scheduler and determinism

There is no fixed 0.005-hour global scan. Each cell schedules migration,
division, and death times. Density/growth updates are lazy and local occupancy
changes refresh only affected neighborhoods. Generation values invalidate stale
queue entries. Stable slots are reused; the store is not compacted or globally
sorted per event.

Simultaneous migrations run as proposal → stable conflict ordering → atomic
commit. A multi-voxel footprint is reserved as one proposal. Conflict priority
is keyed by seed, exact or configured time bucket, UID, event kind/sequence, and
finally UID. Per-cell random draws use `(seed, uid, event_type, event_sequence)`.
OpenMP parallelizes proposal creation only; deterministic ordering controls the
commit. Tests confirm equal final checksums for one and four threads.

The legacy 70×70 density-dependent migration activation maps to a configurable
70³ unique-anchor query. Full density blocks aggregate directly; boundary
blocks inspect compact per-anchor entries. Query centers are cached on
configurable 32³ activation blocks and invalidated by density generation, so
the implementation neither scans every cell nor traverses 70³ voxels per cell.

## Configuration and build

The strict versioned profile is `configs/atcg3d_legacy_like_v1.cfg`. Unknown
keys, unsupported strategies, NaN/infinity, negative times, invalid ranges, and
contradictory domain settings fail with nonzero exit status. CLI overrides use
repeatable `--set key=value`. Strategy strings are validated once at startup;
hot loops use concrete types and contain no string lookup or per-cell virtual
dispatch.

Typical builds:

```sh
cmake -S . -B build -DATCG_BUILD_LEGACY_2D=ON -DATCG_BUILD_3D=ON
cmake --build build -j
ctest --test-dir build --output-on-failure

cmake -S . -B build-3d \
  -DATCG_BUILD_LEGACY_2D=OFF \
  -DATCG3D_ENABLE_HDF5_CHECKPOINT=ON \
  -DATCG3D_ENABLE_VTKHDF=ON
```

The second command requires development installations of HDF5 C++ and VTK
with `IOHDF`. A build without those options retains the simulation core and
fails clearly if the corresponding output capability is requested.

Run and resume examples:

```sh
./build-3d/atcg3d --config configs/atcg3d_legacy_like_v1.cfg \
  --set output.enabled=true --set output.directory=my_run

./build-3d/atcg3d --config configs/atcg3d_legacy_like_v1.cfg \
  --resume my_run/checkpoints/checkpoint_0000000000010000.h5 \
  --set simulation.end_time_hours=48 --set output.directory=resumed_run
```

Resume may change end time, event cap, thread count, and output policy. It
rejects silent changes to biological or numerical dynamics.

## Parameters still requiring scientific calibration

- 3D carrying-capacity scale and r/K limits.
- Directional cone radius, angle, threshold, and block-estimator edge.
- Whether persistent r migration should repeat density filtering.
- Distance correction for axis/planar/spatial lattice steps.
- Large-cell display radius versus physical voxel size.
- 3D initial sphere/shell radius, thickness, and cell-stage mixture.
- Division timing and growth/death delays after mapping from the 2D model.
