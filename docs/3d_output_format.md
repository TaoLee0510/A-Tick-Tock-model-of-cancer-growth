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
preview.vtkhdf.series
full.vtkhdf.series
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

No Verts, Lines, Polys, sphere meshes, or footprint voxels are stored. The six
required arrays plus float coordinates use 31 uncompressed bytes per biological
cell before HDF5 metadata/alignment. Optional growth or migration arrays are not
enabled in the v1 profile.

## Preview and full policy

Preview selects at most `output.preview_max_cells` UIDs with the smallest
stable 64-bit hashes, breaking ties by UID. If N≤K every cell is included.
Sampling is deterministic for the same seed/UID set and never calls simulation
RNG. Full contains every live biological cell and all required fields.

A due full frame always causes a preview frame at the exact same simulation
step/time before the full frame. Preview may contain additional times.
`time` values in `.series` are real model hours, not ordinal frame numbers.

Each frame is written to `*.tmp`, closed and checked, then atomically renamed.
Only after rename is the corresponding standard ParaView `.series` JSON
atomically replaced. Manifests therefore never reference temporary or partial
frames. `run.json` embeds the complete effective configuration and schema
version.

## Visualization clients

Native ParaView can open either series directly. From a ParaView Python shell
or `pvpython`:

```sh
pvpython visualization/paraview_load.py RUN_DIRECTORY --quality preview
```

The independent trame application is:

```sh
pvpython visualization/viewer/app.py RUN_DIRECTORY --port 8080
```

It renders in the ParaView backend; 10-million-point full data is not sent to
the browser. Drag/play loads preview only. A 250 ms idle debounce then loads a
full frame only when its time exactly matches. Every slider input invalidates
the prior full token, so 100 rapid inputs do not queue 100 full reads. Camera
position, focal point, up vector, and parallel scale are retained across frame
switches. The viewer polls atomic series files to discover frames from a live
simulation.

The viewer requires ParaView Python, trame, trame-vtk, and trame-vuetify. These
packages are optional and are not simulation dependencies.
