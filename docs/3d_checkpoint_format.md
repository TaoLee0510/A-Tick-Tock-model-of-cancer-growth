# ATCG3D checkpoint format

Checkpoint files are versioned HDF5 state files and are not visualization
frames. A newly written self-contained base uses schema v6. The default
`journal_delta_hdf5_v2` child uses schema v8 and references its immediately
preceding checkpoint. Readers retain deterministic resume compatibility with
schema-v7 field deltas and the older v4/v5 row-delta chain. Numeric
datasets and attributes use explicit little-endian fixed-width HDF5 types;
native types are used only for the in-memory transfer.

All one-dimensional numeric datasets are chunked, protected by Fletcher32, and
use shuffle+deflate according to `output.storage.hdf5_compression_level`
(default 1). The chunk target is `hdf5_chunk_elements` (default 262144).

The writer creates `NAME.h5.tmp` in the destination directory, globally
flushes and closes it, locally verifies that temporary file's schema, kind,
and row counts, and then atomically renames it to `NAME.h5`. It does not
recursively materialize the parent chain during every write. The strict reader
performs complete chain and reconstructed-checksum validation on resume. A
failed write or verification removes the temporary file. Existing final
checkpoints are not overwritten.

Periodic files use
`checkpoint_<completed_events>_time_<float64-bits>.h5`. The hexadecimal
binary64 time suffix keeps checkpoints unique when several requested output
hours pass without a biological event and `completed_events` is unchanged.
The recovery reader also accepts the earlier
`checkpoint_<completed_events>.h5` name.

## Base and delta policy

`hdf5_base_v6_slot_journal_v8` writes a self-contained base first. During
normal evolution, `CellStore3D` records one final mutation for each changed
stable slot plus the exact ordered free-list push/pop operations. At a
checkpoint boundary the simulator moves only that journal, global clock/stats,
current vascular state, and the lineage suffix into the writer job. It does
not freeze or scan every live cell. A schema-v8 child contains:

- `/cells/changed`: complete current rows for slots alive at the boundary;
- `/cells/removed_slot`: sorted slots that are dead at the boundary;
- `/cells/free_list_mutations/{kind,slot}`: exact ordered free-list changes;
- full current stats and vascular state;
- lineage edges after `lineage_prefix_count`.

Repeated mutations of one slot within an interval are coalesced to its final
state. Slot reuse is represented by the final live row plus its free-list
operations, so UID remains independent from slot. The reader recursively
reconstructs the parent, applies removals/changes/free-list operations in
validated order, appends lineage, then verifies the complete state checksum.
A new base is forced when either of these holds:

- `checkpoint_base_every_hours` has elapsed (default 168 h);
- `checkpoint_max_delta_chain` is reached (default 168);

`delta_full_ratio` remains accepted in schema-v3 YAML for older run metadata,
but schema-v8 does not convert high scheduling churn into an hourly full
checkpoint.

Schema-v8 `/meta` adds `kind=slot_journal_v1`, `parent_file`,
`parent_state_checksum`, `parent_time_hours`, `chain_length`, and
`lineage_prefix_count`. Unsafe parent names, cycles, missing parents, excessive
chains, parent checksum/time mismatches, invalid stable slots, inconsistent
free-list operations, non-append lineage, and corrupt reconstructed checksums
are fatal.

The older `hdf5_base_v6_field_delta_v7` strategy remains readable. Each v7
child records births, sorted removed UIDs, and only the individual fields that
changed on surviving cells. It stores the exact current global state, stable
free-list, vascular state, and lineage suffix after the parent's prefix.

Schema-v7 `/meta` adds `kind=field_delta_v2`, `parent_file`,
`parent_state_checksum`, `parent_time_hours`, `chain_length`, and
`lineage_prefix_count`. `/cells/births` contains complete rows for new UIDs,
`/cells/removed_uid` contains removals, and `/cells/updates` contains sorted
`uid`, stable `slot`, and a 25-bit `field_mask`. Every field column under
`updates` is packed: it has one value only for rows whose mask contains that
field. Anchor uses one bit and three equally packed `x/y/z` columns. Parent
paths must be a filename in the same directory; cycles, missing parents,
excessive chains, parent checksum/time mismatches, non-append lineage,
stable-slot changes, unknown/empty masks, packed-column length mismatches, and
corrupt reconstructed checksums are fatal. Resume accepts a v4/v6 base or any
v5/v7/v8 chain tip and reconstructs the exact state before rebuilding sparse
runtime indexes. Derived sparse grids, density indexes, activation-class
caches, and indexed event-heap positions are rebuilt and validated; they are
not serialized as dense state.

## Metadata and statistics

`/meta` contains scalar attributes:

```text
schema_version = 6 or 8        UInt32
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
migration_attempts migration_commits
migration_swap_waits migration_swap_attempts
migration_swap_commits migration_swap_rejections
divisions deaths conflict_rejections
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
swap_wait_state pending_swap_direction           UInt8
inherent_growth_rate density_growth_rate
migration_rate normal_migration_rate division_work_remaining
death_deadline                                   Float32
next_migration_time migration_activation_end_time
next_division_time
last_update_time swap_ready_time                 Float64
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

`swap_wait_state=1` is persistent biological scheduling state for the
stage-1 singleton crowding exchange. It requires a nonzero
`pending_swap_direction`, a future `swap_ready_time`, and
`next_migration_time == swap_ready_time`. Inactive cells store zero for all
three fields. A thin-layer checkpoint cannot contain a pending z direction.

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
accumulated_eligible_hours remaining_hazard
hazard_last_update_hours hazard_not_before_hours
current_rate_sites_per_30_days
current_density_stress                           Float64
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
accumulated_eligible_hours remaining_hazard
hazard_last_update_hours hazard_not_before_hours
current_rate_sites_per_30_days
current_density_stress                           Float64
schedule_generation                              UInt32
```

`remaining_hazard` is the unconsumed unit-exponential target. The two hazard
times delimit the portion eligible for integration at the current piecewise
constant rate; an eligible zero-rate process has no absolute next event but
retains positive remaining hazard. Density stress is finite in `[0,1]`.
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
An active or transiting tip normally has a nonzero pending direction. The
intentional exception is a blocked-geometry retry: `pending_direction=0` with
a future `next_growth_time` means “wake and recompute feasibility”, not a
terminal tip or corrupt schedule. Schema-v6 validation and resume preserve this
state exactly.
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

The v6/v7 reader rejects missing/unreadable HDF5 objects, unsupported version or
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

Schemas v1, v2, and v3 are explicitly rejected. Schema v1 contains one aggregate
cell scheduling generation and no complete vessel state. Schema v2 contains a
single tumour-wide angiogenesis process and has no lesion identity, per-lesion
processes, lesion refresh scheduler, or vessel source-lesion provenance.
Schema v3 adds the per-lesion topology but lacks remaining hazard, current
rate, density stress, and integration times. Silently upgrading any of these
formats could resample the next density-modulated arrival and would not be a
deterministic resume.

Schema v4/v5 files remain readable. Their absent swap state and four swap
statistics are restored as zero. New writes always use v6/v7.

The HDF5 test suite covers cell-only and active bidirectional-vessel
checkpoints. In both cases, checkpoint/resume to a later time must have the
same final checksum and counters as an uninterrupted run. It also checks the
lesion tables, dirty observations, and source IDs, then injects v1/v2/v3/future
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
