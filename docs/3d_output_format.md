# ATCG3D visualization output

## Run directory

An output-enabled run uses:

```text
run.json
metrics/final.json
checkpoints/checkpoint_*.h5
lineage/edges.csv
viz/preview/frame_*.vtkhdf
viz/full/frame_*.vtkhdf
viz/vessels/frame_*.vtkhdf
preview.vtkhdf.series
full.vtkhdf.series
vessels.vtkhdf.series
```

The 3D executable never invokes legacy PNG code and never creates PNG frames.
Visualization snapshots are read-only views over `CellStore3D`; their layout is
not a checkpoint ABI.

## VTK-HDF point schema

Every preview or full frame is official VTK-HDF written by `vtkHDFWriter` from
a points-only `vtkPolyData`:

| Data | Type | Meaning |
|---|---|---|
| `Points` | Float32 `[N,3]` | one center per biological cell |
| `cell_id` | UInt64 | stable UID |
| `clone_id` | UInt32 | clone identity |
| `cell_type` | UInt8 | r=1, K=2 (legacy-compatible labels) |
| `stage` | UInt8 | large=0, small=1, ultrasmall=2 |
| `viability` | UInt8 | biological viability flag |
| `display_radius` | Float32 | rendering scale hint |

Each cell frame also has one dataset-level FieldData value:

| Field array | Type | Meaning |
|---|---|---|
| `total_cell_count` | UInt64 scalar | all live biological cells at this time, including cells omitted by preview sampling |

No Verts, Lines, Polys, sphere meshes, or footprint voxels are stored. The six
required arrays plus float coordinates use 31 uncompressed bytes per biological
cell before HDF5 metadata/alignment. Optional growth or migration arrays are not
enabled in the v1 profile.

Lattice anchors are voxel-corner coordinates. A stage-1/stage-2 cell at anchor
`p` and a vessel node occupying voxel `p` are both written at `p + (0.5,0.5,0.5)`.
A stage-0 cell is written once at the center of its 2x2x2 footprint,
`p + (1,1,1)`. This common center convention keeps the cell and vessel overlays
spatially aligned.

## VTK-HDF vessel schema

Vasculature is a separate official VTK-HDF `vtkPolyData`; it is never mixed
into the cell files. Every live centerline node is one Float32 point and every
non-root node contributes one VTK Line connecting its parent to itself. It has
no Verts, Polys, sphere meshes, or replicated tube surface in storage. ParaView
creates the display tube at render time.

| Point array | Type | Meaning |
|---|---|---|
| `node_id` | UInt64 | stable vessel-node UID |
| `vessel_id` | UInt64 | root/network identity |
| `branch_role` | UInt8 | root=0, inward=1, outward=2 |
| `perfused` | UInt8 | compatibility field; generated vessels are always 1 |
| `diameter_voxels` | Float32 | configured biological diameter in lattice voxels |
| `radius_voxels` | Float32 | `diameter_voxels / 2`; ParaView Tube absolute-radius scalar |

`run.json` schema version 3 declares the cell and vessel topology, both array
catalogs (including the cell FieldData catalog), and `vessel_series`. The cell
dataset remains strictly points-only:
one point per biological cell, including large cells whose occupancy footprint
contains eight voxels.

## Preview and full policy

Preview selects at most `output.preview_max_cells` UIDs with the smallest
stable 64-bit hashes, breaking ties by UID. If N≤K every cell is included.
Sampling is deterministic for the same seed/UID set and never calls simulation
RNG. Full contains every live biological cell and all required fields.

A due full frame always causes a preview frame at the exact same simulation
step/time before the full frame. Preview may contain additional times.
`time` values in `.series` are real model hours, not ordinal frame numbers.
Every preview frame has one vessel frame at exactly the same time, including an
empty Lines frame before angiogenesis begins. Full switches reuse that exact-time
vessel frame, so cell and vascular state cannot drift on the viewer timeline.

Each frame is written to `*.tmp`, closed and checked, then atomically renamed.
Only after rename is the corresponding standard ParaView `.series` JSON
atomically replaced. Manifests therefore never reference temporary or partial
frames. `run.json` embeds the complete effective configuration and schema
version. It also stores the normalized dynamics configuration used for exact
resume-output compatibility checks.

## Resume and append semantics

With `run.mode: resume` and output enabled, `output.directory` must name the
existing run directory. Before writing anything, the output manager strictly
loads `run.json` and all three `.series` catalogs. It rejects a schema or
dynamics mismatch, unknown/non-canonical frame paths, `.tmp` references,
missing frame files, non-increasing times, mismatched preview/vessel times, or
a full time without an exact preview time. On the first restored snapshot it
then reconciles the run directory to the selected checkpoint instead of simply
rejecting a run whose process had continued after that checkpoint.

The retained visualization prefix contains only catalog entries at or before
the checkpoint time. The retained lineage is verified edge-for-edge against
the checkpoint lineage; a mismatch in the common prefix is fatal. An existing
lineage shorter than the checkpoint is completed by normal append, while an
existing lineage tail beyond the checkpoint is rolled back. Preview, full,
vessel, and lineage catalogs are rewritten through temporary files and atomic
rename. Preview and vessel prefixes must remain paired, and every retained full
time must still have its exact preview.

No superseded data is deleted. Future catalog frames, canonical orphan frames,
future checkpoints, a stale final metrics file, lineage tail, and recognized
temporary artifacts are moved into a unique
`recovery/checkpoint_<completed_events>[_NNNN]/` directory. Its directory
layout mirrors the run directory and `recovery.json` records the checkpoint,
retained counts, and quarantined paths. Catalogs are published before referenced
future frames are moved, so readers never see a committed catalog pointing to a
quarantined frame. The next resumed frame reuses the first free canonical index;
the preserved recovery copy is not overwritten.

Periodic preview/full/checkpoint schedules restart at the first regular
boundary strictly after the restored time, so the checkpoint instant is not
emitted twice. Fresh `run.mode: new` runs continue to reject an output directory
that already has a `run.json`.

## Visualization clients

Native ParaView can open either series directly. From a ParaView Python shell
or `pvpython`:

```sh
pvpython visualization/paraview_load.py RUN_DIRECTORY --quality preview \
  --cell-radius-scale 1.0 --vessel-radius-scale 1.0
```

The independent trame application is:

```sh
pvpython visualization/viewer/app.py RUN_DIRECTORY --port 8080
```

It renders cells and Tube-filtered vascular Lines in the ParaView backend;
10-million-point full data is not sent to the browser. r cells use a fixed
green categorical color, K cells fixed red, and vessels fixed blue. Vessel
radius can be scaled without changing stored geometry. Whole, one-sided cut,
and adjustable X/Y/Z slab modes clip both cells and vessel tubes in the backend
while the camera remains freely rotatable. Drag/play loads preview cells plus the matching
vessel frame. The toolbar reports both the true live-cell count from
`total_cell_count` and the number of sampled points currently displayed. A
250 ms idle debounce then loads a full cell frame only when its
time exactly matches. Every slider input invalidates the prior full token, so
100 rapid inputs do not queue 100 full reads. Camera position, focal point, up
vector, and parallel scale are retained across cell/vessel frame switches. A
new slider event invalidates the preceding token both before and after a full
read, so a stale completion is never published and no further stale vessel/full
work is queued. ParaView's individual `UpdatePipeline()` call is synchronous
and cannot be interrupted safely once inside VTK; on the measured 10-million
frame that non-preemptible interval was 0.10 s. The viewer polls all three
atomic series files to discover frames from a live simulation.

Cell splats use a server-side Calculator array
`display_radius * point_size_scale`, with Point Gaussian `ScaleByArray`
explicitly enabled and transfer remapping disabled. Thus the YAML stage radii
and the interactive point-size control both affect the actual rendered radius;
the full point array remains in ParaView rather than the browser.

The viewer requires ParaView Python plus the packages pinned in
`visualization/viewer/requirements.txt`. They are optional viewer dependencies,
not simulation dependencies. Install them into a Python environment with the
same major/minor version as `pvpython`, then expose that environment's
`site-packages` to `pvpython` if ParaView does not bundle them.
