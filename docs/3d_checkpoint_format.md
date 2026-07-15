# ATCG3D checkpoint format

Checkpoint files are versioned HDF5 state files and are not visualization
frames. Schema v2 is the first schema that can resume angiogenesis. Numeric
datasets and attributes use explicit little-endian fixed-width HDF5 types;
native types are used only for the in-memory transfer.

The writer creates `NAME.h5.tmp` in the destination directory, globally
flushes and closes it, reads it back through the strict v2 reader, and then
atomically renames it to `NAME.h5`. A failed write or verification removes the
temporary file. Existing final checkpoints are not overwritten.

## Metadata and statistics

`/meta` contains scalar attributes:

```text
schema_version = 2             UInt32
dimension = 3                  UInt32
completed_events               UInt64
time_hours                     Float64
next_uid state_checksum        UInt64
effective_config_json          UTF-8 string
dynamics_config_json           UTF-8 string
```

`/stats` contains UInt64 scalar attributes for every `SimulationStats3D`
counter:

```text
migration_attempts migration_commits divisions deaths conflict_rejections
angiogenesis_seed_attempts angiogenesis_roots angiogenesis_seed_rejections
vessel_growth_attempts vessel_growth_commits vessel_anastomoses
vascular_displacements
```

End time, maximum events, thread count, run/output paths, intervals, sampling
limits, and display radii may change on resume. `dynamics_config_json` excludes
those run-policy fields. A change to geometry, biology, density, migration,
division, death, RNG seed, angiogenesis, or other numerical dynamics is
rejected.

## Cell table

`/cells` contains equal-length, one-dimensional columns:

```text
slot                                            UInt32
x y z                                           Int32
uid parent_uid event_sequence                   UInt64
clone_id                                        UInt32
type stage viability flags last_direction       UInt8
inherent_growth_rate density_growth_rate
migration_rate normal_migration_rate division_work_remaining
death_deadline                                   Float32
next_migration_time migration_activation_end_time
next_division_time
last_update_time                                Float64
migration_schedule_generation
division_schedule_generation
death_schedule_generation                       UInt32
```

The group also contains `slot_count` as a UInt64 scalar attribute and a
`free_slots` UInt32 dataset. `slot` maps each UID-sorted live row back to its
stable store slot. `free_slots` records the exact LIFO order, and the live and
free slot sets must form a disjoint, complete partition of `[0, slot_count)`.

`migration_rate` is the persistent inherent/base rate;
`normal_migration_rate` is the ordinary low-density rate. Active state requires
both the `kMigrationActive` flag and an activation-end time strictly after the
checkpoint clock; normal state requires a zero activation-end time.

Migration, division, and death generations are intentionally separate. The v1
aggregate `schedule_generation` is not written. Queue entries are reconstructed
from the three next-event times and their matching generations, so a stale event
of one kind cannot invalidate or recreate another kind. Restore rejects
activation/rate/time combinations that omit a required normal/active migration
or pair an active flag with an invalid end time. The migration generation also
guards the activation-end event; heap maintenance state is derived and is not
checkpointed.

`division_work_remaining` is the unfinished unit-rate work for the already
sampled cell cycle. Together with `density_growth_rate` and `last_update_time`,
it permits exact lazy progress integration across density changes without
drawing a new cycle after migration or a neighborhood refresh.

Activation-end, next-migration, next-division, and last-update remain Float64
in HDF5 for schema stability; the already-float death deadline is Float32 as
listed above. `CellStore3D` uses float internally for all five values. Public
snapshots promote them to double, and the four Float64 datasets write that
promoted representation without inventing precision. Scheduled times are
rounded toward positive infinity, while last-update uses nearest rounding.
Checkpoint/resume therefore preserves the exact already-quantized state.

Rows are sorted by cell UID before writing, but cell slots and free-list order
are restored exactly. UID remains the persistent biological identity; slot is
also continuation state because deterministic conflict scheduling and the next
free-list reuse depend on it. The sparse chunk grid and density index are
rebuilt and checked for footprint conflicts. No dense 3D occupancy array is
checkpointed.

## Vasculature state

`/vasculature` has these UInt64 attributes and dataset:

```text
next_vessel_id next_node_uid next_tip_uid
perfused_vessel_ids                              UInt64[N]
```

`perfused_vessel_ids` is strictly sorted and unique. `/vasculature/process`
stores the complete dynamic Poisson-process state:

```text
eligible                                         UInt8 boolean
next_seed_time_hours eligibility_started_hours
accumulated_eligible_hours                       Float64
event_sequence attempted_events committed_roots
rejected_events                                  UInt64
schedule_generation                              UInt32
```

Each Poisson arrival is one attempted site and commits at most one root, so the
reader validates
`attempted_events = committed_roots + rejected_events`.

`/vasculature/nodes` is a typed structure-of-arrays table:

```text
x y z                                             Int32
uid parent_uid vessel_id                          UInt64
parent_node_slot                                  UInt32
role perfused                                     UInt8
diameter_voxels                                   Float32
created_time_hours                                Float64
```

`/vasculature/tips` contains:

```text
x y z bias_x bias_y bias_z target_x target_y target_z Int32
uid vessel_id current_node_uid event_sequence         UInt64
current_node_slot schedule_generation                 UInt32
role status perfused last_direction pending_direction UInt8
diameter_voxels speed_voxels_per_hour
max_length_voxels grown_length_voxels                  Float32
next_growth_time                                       Float64
```

Node and tip rows retain the stores' stable-slot order; they are never sorted
by UID. Consequently `parent_node_slot` and `current_node_slot` remain valid
and are cross-checked against their corresponding UIDs and vessel IDs. Terminal
tips are retained, including their inert last pending-direction/time fields,
because those fields participate in the deterministic state checksum.

Vessel voxel masks and the vascular influence field are derived state. On
resume they are rebuilt from the node segments, diameter, role, and perfusion
set, then validated against biological-cell occupancy. Scheduler entries for
active tips and the pending Poisson seed event are reconstructed from their
stored time and generation.

## Lineage

`/lineage` remains an append-edge table:

```text
birth_time                                       Float64
child_uid parent_uid                             UInt64
clone_id                                         UInt32
type                                             UInt8
```

The legacy fixed-width 150-ancestor matrix is not used.

## Validation and compatibility

The v2 reader rejects missing/unreadable HDF5 objects, unsupported version or
dimension, wrong fixed-width types, non-scalar attributes, unequal column
lengths, NaN/invalid times, invalid enums/directions/flags, duplicate or zero
IDs, invalid next-ID counters, invalid or incomplete cell-slot partitions,
broken node/tip references, inconsistent
perfusion or angiogenesis counters, configuration mismatch, and conflicting
reconstructed cell/vessel footprints. After those typed checks, the reader
restores a fresh `Simulation3D` and recomputes the complete state checksum; it
therefore also rejects schema-valid finite values that were changed without a
matching checksum update.

Schema v1 is explicitly rejected. It contains one aggregate cell scheduling
generation and no angiogenesis, node, tip, perfusion, or vessel scheduler
state, so silently treating it as v2 would not be a deterministic resume.

The HDF5 test suite covers cell-only and active bidirectional-vessel
checkpoints. In both cases, checkpoint/resume to a later time must have the
same final checksum and counters as an uninterrupted run. It also injects an
invalid schema version, a mismatched cell-column length, an invalid file, and
a dynamics mismatch and requires a clear nonzero failure. The state checksum
includes `slot_count`, each live row's slot, exact free-list order, cells
(including remaining division work), scheduler generations, every
simulation statistic, all lineage edge fields, and vasculature. Both the reader
and the resume application recompute it from restored state before continuing.

## Resuming an older checkpoint in an existing run directory

A checkpoint may be valid even when the same process wrote later visualization,
lineage, metrics, or checkpoint files before stopping. Output-enabled resume
therefore performs a checkpoint-prefix reconciliation on its first restored
snapshot. Catalog entries later than the restored clock and lineage edges after
the exact checkpoint lineage are not treated as part of the resumed history.
The common lineage prefix must match edge-for-edge; a mismatch remains a hard
error.

Superseded artifacts are preserved under
`recovery/checkpoint_<completed_events>[_NNNN]/`, including future and orphan
VTK-HDF frames, the lineage tail, future checkpoints, stale final metrics, and
recognized temporary files. `recovery.json` records what was retained and
quarantined. The live `.series` files and `lineage/edges.csv` are replaced
atomically with the checkpoint prefix before resumed output reuses any frame
index. This recovery directory is forensic data and is not read as simulation
state; only the explicitly selected HDF5 checkpoint supplies continuation
state.
