# ATCG3D checkpoint format

Checkpoint files are versioned HDF5 state files, separate from VTK-HDF
visualization frames. They are written to a temporary file, globally flushed
and closed, then atomically renamed.

## Metadata

`/meta` attributes include:

- `schema_version=1`, `dimension=3`;
- completed event count and current model time;
- next UID and state checksum;
- migration/division/death/conflict statistics;
- complete effective configuration JSON;
- dynamics-only configuration JSON used for resume compatibility.

Changing end time, maximum events, thread count, or output policy is allowed on
resume. Any change to geometry, density, migration, division, death, RNG seed,
or other dynamics is rejected.

## Cell columns

`/cells` contains one-dimensional fixed-width datasets:

```text
x y z                         Int32
uid parent_uid event_sequence UInt64
clone_id schedule_generation  UInt32
type stage viability flags last_direction UInt8
inherent_growth_rate density_growth_rate migration_rate Float32
death_deadline Float32
next_migration_time next_division_time last_update_time Float64
```

These columns include all state that affects future evolution. Scheduler queue
entries are reconstructed from stored next-event times and generations;
generation checks invalidate stale events. The sparse chunk grid is rebuilt
from cell anchors/stages and validated for footprint conflicts. No dense 3D
grid is checkpointed.

## Lineage

`/lineage` is an append-edge representation with equal-length columns:

```text
birth_time Float64
child_uid parent_uid UInt64
clone_id UInt32
type UInt8
```

The run directory also exposes newly observed edges as
`lineage/edges.csv`. The legacy fixed-width 150-ancestor matrix is not used.

## Recovery guarantees

The reader rejects unreadable HDF5, unsupported version/dimension, unequal
column lengths, duplicate/zero UIDs, configuration mismatch, and conflicting
reconstructed footprints. Resume continues strictly after the stored completed
event state. The test suite compares a 5-hour checkpoint/resume-to-10-hour run
against an uninterrupted 10-hour run and requires identical checksum, event
count, and migration commits.
