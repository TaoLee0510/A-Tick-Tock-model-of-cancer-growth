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

The legacy r-to-K rule is applied only after daughter placement succeeds. For
an r mother, production evaluates a stage-aware 70³ unique-anchor density at
the pre-division anchor via the incremental block estimator. At density
`>=0.50`, the daughter converts to K with probability `0.05`. The decision is
stateless and keyed by `(seed, mother_uid, division_event_sequence)`. A
converted daughter receives a fresh deterministic K growth/migration-rate
sample keyed by its own UID. Inherited rates are capped at `1.3171805` for r
and `0.99505180` for K. All thresholds, block/window edges, probability, and
caps are strict YAML fields. Failed division attempts change neither cell
state nor event sequence.

Density-dependent growth uses a 6×6×6 anchor window by default. The initial 3D
capacity mapping multiplies the 2D limits/capacities by six. All values are
configuration fields and require scientific calibration; the factor six is a
starting interpretation, not a validated biological constant.

Initialization maps the disk and annulus to a sphere and spherical shell. It
generates actual sphere/shell candidate anchors by z/y slices rather than
scanning the global sparse domain or a dense X×Y×Z array. Footprint availability
is checked before every placement.

The angiogenesis surface index is updated only near occupancy changes. Its
centroid sums are incremental. Root sampling excludes sealed internal cavities
by retaining only faces with an unobstructed outward-normal axis ray, then uses
a deterministic stable-hash top-K sample. It is O(S + S log K), K≤4096, and
does not copy or sort the complete exposed-face set.

When a vessel becomes perfused, growth refresh does not query the global AABB
of the whole curved network. Capsule voxels are mapped to the union of density
blocks intersecting their configured relief cutoff; each block is visited once,
then anchors are filtered by exact Euclidean distance and nonzero relief. This
keeps dirty work local and prevents cells lying inside a large network bounding
box but outside every vessel gradient from changing schedule or RNG state.

## Storage, identity, and memory

`CellStore3D` is a typed structure-of-arrays store with stable `uint32` slots
and a free list. UID and parent UID are immutable 64-bit biological identities;
UID never equals slot by contract and never changes after slot reuse. Position
anchors are three `int32` columns. Type, stage, viability, flags, and last
direction are byte columns. Rates, remaining division work, the death deadline,
and the migration/activation/division/update time columns are floats internally.
The public snapshot API exposes times as doubles; HDF5 stores the four
activation/migration/division/update columns as Float64 and the death deadline
as Float32, preserving their already-quantized values.
Scheduled times round upward to the next representable float so quantization
cannot place a new event before the current clock. Last-update times use nearest
rounding, and lazy progress integrates between stored quantized values.
Event sequence, the three event-kind schedule generations, and the remaining
float division-work value are explicit mutable state.

Checkpoint/resume preserves the allocated slot count, every live cell's
original slot, and the free-list LIFO order exactly. These are deterministic
continuation state: compacting live cells during restore would change the slot
reused by a later birth and could change subsequent conflicts.

The logical column width is 98 bytes per allocated slot after adding persistent
cell-cycle work plus finite migration-state rate/timing. Older memory
measurements predate these fields and must not be
reported as the current ten-million-cell result; the scale benchmark reports
the actual current allocation when rerun.

## Event scheduler and determinism

There is no fixed 0.005-hour global scan. Each cell schedules migration,
division, and death times. Density/growth updates are lazy and local occupancy
changes refresh only affected neighborhoods. Elapsed work is deducted with the
previous density rate; density changes adjust the predicted completion time but
never redraw the cell cycle. Generation values invalidate stale
queue entries. When the heap grows by 25% (with a small fixed minimum slack),
it is deterministically rebuilt from canonical next times and generations;
this bounds inert entries without consuming RNG, changing generations, or
performing a full rebuild per event. Stable slots are reused; the cell store is
not compacted or globally sorted per event.

Simultaneous migrations, divisions, and stage recoveries run as proposal →
stable conflict ordering → atomic commit. Daughter and recovered large-cell
footprints are reserved as complete voxel sets; a co-location source group also
has an exclusive recovery lock. A losing division only schedules its retry and
cannot fall through to r-cell death or K-cell co-location. Same-time deaths are
judged against one density snapshot, committed together, and trigger one
stage-recovery batch. Conflict priority is keyed only by seed, exact or
configured time bucket, UID, and event kind, with UID as the final tie-break.
Per-cell biological random draws separately use
`(seed, uid, event_type, event_sequence)`. OpenMP parallelizes read-only death
decisions plus division, migration, and vessel-growth proposal creation;
deterministic ordering controls every commit. `simulation.threads` is the hard
maximum. In `parallel.mode: adaptive_cells_and_events_v1`, a configurable
population threshold selects a fraction of that maximum and
`min_events_per_thread` independently caps workers for small same-time batches.
The lower of those limits is used, bounded by the available OpenMP workers.
Thread selection consumes no simulation RNG. Tests cover threshold boundaries
and confirm equal final checksums for one and multiple threads.

Same-time divisions also build immutable proposals before any contender
mutates occupancy. Proposals reserve the complete daughter/shape-reduction
footprint and lock a stage-2 source co-location group; a conflict loser only
reschedules and cannot fall through to the r-cell death or K-cell co-location
fallback. Same-time deaths are judged against one density snapshot, committed
as a batch, and followed by one stage-recovery proposal batch. Recovery uses
the same whole-footprint reservations and one-move-per-co-location-group lock.

The legacy 70×70 density-dependent migration trigger maps to a configurable
70³ unique-anchor query (70² in exact thin-layer mode). Full density blocks
aggregate directly; boundary
blocks inspect compact per-anchor entries. Query centers are cached on
configurable 32³ activation blocks and invalidated locally by anchor changes.
Simulation keeps a derived two-bit activation-class cache per occupied query
block (small/ultrasmall and large). Local occupancy changes recompute only the
bounded set of affected block classes, normally at most 27 for the 70/32
mapping. Residents of a block are visited in bulk only when the corresponding
class crosses the threshold; moved and newly born anchors are always refreshed
directly so cross-block moves cannot retain the source flag. Stage-recovery
before/after anchors are included in the dirty set. The class cache is rebuilt
deterministically after initialization/resume and is not checkpoint state.
Migration is not an on/off gate. In normal state every cell has a queued
migration event: r cells use `Beta(5,5) * 0.5`, while K cells use their own
normal/base rate, and both choose uniformly among feasible directions. A high
density threshold crossing starts a finite active interval. During it the
inherent/base rate is used and r cells use the directional density/persistence
rule. Its duration is
`Beta(alpha=0.005,beta=0.011666666666666667) * remaining_cycle_hours`, whose
configured mean fraction is 0.30. An exact activation-end event returns the
cell to normal state and resamples the r normal rate; low density never ends an
active interval early, and later high density may activate the cell again.
Migration and activation-end events share a generation, so stale events are
invalid without a scan. The steady path therefore neither scans every cell,
all residents of unchanged query blocks, nor 70³ voxels per cell.

## Configuration and build

The complete production profile is
`configs/atcg3d_legacy_2d_mapped_v2.yaml`; the separate
`configs/atcg3d_smoke_test_v2.yaml` profile is only for small CI/developer
runs. Schema v2 rejects missing, unknown, and duplicate keys, loose boolean
spellings, integer overflow, NaN/infinity, negative times, invalid ranges, and
contradictory settings. There are no command-line value overrides: all model,
run, output, and resume parameters come from the selected YAML file. Strategy
strings are validated once at startup; hot loops use concrete typed values and
contain no string lookup or per-cell virtual dispatch. The 3D target requires
the `yaml-cpp` development package.

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

Run and validation examples:

```sh
./build-3d/atcg3d --config configs/atcg3d_legacy_2d_mapped_v2.yaml --dry-run
./build-3d/atcg3d --config configs/atcg3d_legacy_2d_mapped_v2.yaml
```

To resume, create a YAML run file with `run.mode: resume` and
`run.resume_checkpoint` set to the checkpoint path. End time, event cap,
thread count, and output policy are also edited in that YAML file. Resume
rejects silent changes to biological or numerical dynamics.

## Mapped defaults and scientific calibration

No counterpart parameter is left unset. The current production YAML maps the
2D values proportionally into 3D: density capacity is multiplied by six,
6×6 becomes 6×6×6, 70×70 becomes 70×70×70 with block aggregation, disk/ring
geometry becomes sphere/shell geometry, the 2D directions remain IDs 1–8 of
the 26-direction set, and legacy growth, migration, division, and death
distributions retain their 2D numeric values where their units do not change.
Display radii are 1.0, 0.5, and 0.25 voxels for stages 0, 1, and 2.

"Requires scientific calibration" therefore means that these implemented
mapping assumptions have not yet been fitted to independent 3D experimental
data; it does not mean the simulator is missing values. In particular, future
calibration may revise the factor-six carrying capacities, cone parameters,
distance weighting, activation threshold/block edge, sphere/shell proportions,
division/death timing, and physical voxel/display scale. Angiogenesis has no
2D predecessor, so its trigger volume, Poisson intensity, vessel speeds and
diameter, connection distance, direction bias, and density-relief gradient are
explicit starting hypotheses and should be fitted separately.
