# ATCG3D checkpoint format

Checkpoint files are versioned HDF5 state files and are not visualization
frames. Schema v3 is the first schema that can resume independent, persistent
angiogenesis processes for multiple primary or metastatic lesions. Numeric
datasets and attributes use explicit little-endian fixed-width HDF5 types;
native types are used only for the in-memory transfer.

The writer creates `NAME.h5.tmp` in the destination directory, globally
flushes and closes it, reads it back through the strict v3 reader, and then
atomically renames it to `NAME.h5`. A failed write or verification removes the
temporary file. Existing final checkpoints are not overwritten.

Periodic files use
`checkpoint_<completed_events>_time_<float64-bits>.h5`. The hexadecimal
binary64 time suffix keeps checkpoints unique when several requested output
hours pass without a biological event and `completed_events` is unchanged.
The recovery reader also accepts the earlier
`checkpoint_<completed_events>.h5` name.

## Metadata and statistics

`/meta` contains scalar attributes:

```text
schema_version = 3             UInt32
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
rejected. `effective_config_json` is immutable provenance for the configuration
that created the checkpoint, not the override policy for the resumed run. The
reader nevertheless requires it to exist as a nonempty JSON/YAML mapping and
requires its scalar `schema_name`, `schema_version`, and `profile` identity to
match the requested configuration. A missing, malformed, sequence-valued, or
identity-incompatible provenance record is rejected.

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
stores a compatibility aggregate of all retained lesion processes:

```text
eligible                                         UInt8 boolean
next_seed_time_hours eligibility_started_hours
accumulated_eligible_hours                       Float64
event_sequence attempted_events committed_roots
rejected_events                                  UInt64
schedule_generation                              UInt32
```

The aggregate is derived, not an independent scheduler. Process rows are first
ordered by `lesion_id`, and all additions use checked arithmetic. Boolean
eligibility is the OR of lesion eligibility and the next time is the minimum
among eligible lesions. Each active row contributes
`accumulated_eligible_hours + checkpoint_time - eligibility_started_hours`;
inactive rows contribute their stored accumulated value. Those exact elapsed
values and all three event counters are summed in lesion-ID order, while event
sequence/generation are maxima. If any row is active, the derived aggregate's
eligibility start is the checkpoint time itself; otherwise it is zero. This
rebasing prevents the already-folded active intervals from being counted
twice. The aggregate is recomputed for every snapshot and is never advanced as
a scheduler. The reader requires it to match
`/vasculature/lesions/processes` exactly and rejects numeric/counter overflow.

`/vasculature/lesions` has scalar attributes:

```text
next_lesion_id                                   UInt64
last_refresh_time_hours next_refresh_time_hours  Float64
refresh_schedule_generation                      UInt32
```

`next_lesion_id` is nonzero and greater than every persisted lesion ID.
`last_refresh_time_hours` is not after the checkpoint clock. When angiogenesis
is enabled, `next_refresh_time_hours` is the event-driven sparse lesion-index
refresh deadline and must not precede the checkpoint clock; its generation makes
older queued refresh events stale after restore.

`/vasculature/lesions/core_identity` preserves stable lesion identity across a
sparse topology rebuild:

```text
block_x block_y block_z                          Int32
lesion_id                                        UInt64
```

Rows are strictly sorted by unique block coordinate. A lesion ID may occur in
multiple rows.

`/vasculature/lesions/dirty_blocks` preserves the old observation for each
block changed after the last lesion refresh, so checkpoint restore does not
prematurely expose those changes before the saved refresh deadline:

```text
block_x block_y block_z                          Int32
exists                                           UInt8 boolean
cell_count occupied_voxel_count                  UInt64
biological_volume                                Float64
cell_coordinate_sum_x cell_coordinate_sum_y
cell_coordinate_sum_z
occupied_coordinate_sum_x occupied_coordinate_sum_y
occupied_coordinate_sum_z                       Int64
```

Rows are strictly sorted by unique block coordinate. An absent row
(`exists=0`) has an all-zero payload. A present row has a nonempty observation;
counts, volume, occupancy capacity, and coordinate-sum bounds are validated.
The table is nonempty exactly when a lesion refresh event is pending.

`/vasculature/lesions/processes` contains one row per retained lesion process,
strictly sorted by unique `lesion_id`:

```text
lesion_id event_sequence attempted_events
committed_roots rejected_events                  UInt64
eligible                                         UInt8 boolean
next_seed_time_hours eligibility_started_hours
accumulated_eligible_hours                       Float64
schedule_generation                              UInt32
```

Each Poisson arrival is one attempted site and commits at most one root, so the
reader validates `attempted_events = committed_roots + rejected_events` for
every process. Every current core identity has a process. On lesion merge,
predecessor counters are accumulated into the retained result process before
the predecessor process is erased; therefore only the global sum of committed
roots is compared with total root-node count.

`/vasculature/lesions/source_ownership` is the sparse, versioned link between
immutable vessel provenance and current topology:

```text
source_lesion_id current_lesion_id               UInt64
```

Rows are strictly sorted by unique, nonzero `source_lesion_id`. Identity rows
are omitted: a source that is still a current lesion implicitly owns itself.
For a historical source, nonzero `current_lesion_id` must name a current lesion;
zero means that the source has no living lesion owner. Targets are direct (no
alias chains), and both IDs must be below `next_lesion_id`. On later merges or
removals, every historical row is retargeted to the new current owner or zero.
This table participates in the state checksum and is restored exactly. It does
not rewrite the `source_lesion_id` stored on vessel nodes or tips.

`/vasculature/nodes` is a typed structure-of-arrays table:

```text
x y z                                             Int32
uid parent_uid vessel_id source_lesion_id         UInt64
parent_node_slot                                  UInt32
role perfused                                     UInt8
diameter_voxels                                   Float32
created_time_hours                                Float64
```

`/vasculature/tips` contains:

```text
x y z bias_x bias_y bias_z target_x target_y target_z Int32
uid vessel_id source_lesion_id current_node_uid
event_sequence                                        UInt64
current_node_slot schedule_generation                 UInt32
role status perfused last_direction pending_direction UInt8
diameter_voxels speed_voxels_per_hour
max_length_voxels grown_length_voxels                  Float32
next_growth_time                                       Float64
```

Node and tip rows retain the stores' stable-slot order; they are never sorted
by UID. Consequently `parent_node_slot` and `current_node_slot` remain valid
and are cross-checked against their corresponding UIDs, vessel IDs, and source
lesion IDs. All nodes in a vessel have one source lesion; tips must match their
current node. A vessel keeps its historical nonzero source ID after a lesion
merge, so that ID need only be below `next_lesion_id` and need not name a
current process, but every such non-current source must have an explicit
`source_ownership` row (including owner zero after removal). A missing
historical ownership row is rejected. Terminal tips are retained, including their inert last
pending-direction/time fields,
because those fields participate in the deterministic state checksum.
Per-lesion active-tip limits resolve each tip's historical source through
`source_ownership`, so tips created before a merge still consume the retained
current lesion's budget.

Vessel voxel masks and the vascular influence field are derived state. On
resume they are rebuilt from the node segments, diameter, role, and perfusion
set, then validated against biological-cell occupancy. Scheduler entries for
active tips, pending per-lesion Poisson seed events, and the lesion refresh are
reconstructed from their stored times and generations.

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

The v3 reader rejects missing/unreadable HDF5 objects, unsupported version or
dimension, wrong fixed-width types, non-scalar attributes, unequal column
lengths, NaN/invalid times, invalid enums/directions/flags, duplicate or zero
IDs, invalid next-ID counters, invalid or incomplete cell-slot partitions,
broken node/tip/source-lesion references, unsorted or duplicate lesion tables,
invalid dirty-block snapshots, missing or invalid source-ownership tables,
inconsistent per-lesion/aggregate/root counters, invalid lesion refresh state,
inconsistent perfusion, missing/malformed effective-configuration provenance,
configuration mismatch, and conflicting
reconstructed cell/vessel footprints. After those typed checks, the reader
restores a fresh `Simulation3D` and recomputes the complete state checksum; it
therefore also rejects schema-valid finite values that were changed without a
matching checksum update.

Schemas v1 and v2 are explicitly rejected. Schema v1 contains one aggregate
cell scheduling generation and no complete vessel state. Schema v2 contains a
single tumour-wide angiogenesis process and has no lesion identity, per-lesion
processes, lesion refresh scheduler, or vessel source-lesion provenance.
Silently interpreting either as v3 would not be a deterministic resume.

The HDF5 test suite covers cell-only and active bidirectional-vessel
checkpoints. In both cases, checkpoint/resume to a later time must have the
same final checksum and counters as an uninterrupted run. It also checks the
lesion tables, dirty observations, and source IDs, then injects v1/v2/future
schema versions, a mismatched lesion aggregate, an invalid dirty-block boolean,
inconsistent node/tip source lesions, missing/corrupt source ownership,
missing/malformed effective-configuration provenance, a mismatched cell-column
length, an invalid file, and a dynamics mismatch and requires a clear nonzero
failure. The
state checksum
includes `slot_count`, each live row's slot, exact free-list order, cells
(including remaining division work), scheduler generations, every
simulation statistic, all lineage edge fields, lesion-source ownership, and
vasculature. Both the reader and the resume application recompute it from
restored state before continuing.

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
