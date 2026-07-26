# ATCG3D trame viewer

This is a thin client over a ParaView backend. It supports arbitrary 3D camera
rotation/pan/zoom, time slider, play/pause, time and cell count, `cell_type`
coloring, Point Gaussian display, radius scaling, preview/full status, and live
series refresh. A synchronized server-side Tube overlay renders
`vessels.vtkhdf.series`. r cells are fixed green, K cells fixed red, and vessel
tubes fixed blue. Independent r/K/vessel switches allow any class to be hidden.
The top bar is reserved for playback, the time line, and current-frame status;
all camera, visibility, size, and slicing parameters are in the collapsible
side drawer. `Primary drag=Rotate` orbits with the left button, while
`Primary drag=Pan` translates both the camera and focal point with the left
button so the tumour can be placed anywhere in the viewport. Middle-button or
Shift-left drag remains a direct pan shortcut, and wheel/right drag zooms.

The slice panel provides whole, one-sided cut, and slab modes. World X/Y/Z are
available alongside three view-relative planes captured from the current
camera basis: axial is perpendicular to the viewing direction, sagittal is
normal to screen-right, and coronal is normal to screen-up. Selecting a
view-relative orientation captures it through the current camera focal point;
`Recapture current view` replaces it after a new camera pose. Subsequent camera
rotation does not rotate the captured plane, while the offset and thickness
sliders continue moving it along the captured normal. The same backend planes
clip cells, vessel tubes, and the optional
vascular-influence overlay. Vessel radius remains adjustable. The overlay reads
the exact cutoff from `run.json` and shows the approximate centreline reach
`radius_voxels + cutoff_radius_voxels` as a translucent blue tube. The
top status bar distinguishes the true live-cell count
stored in frame FieldData from the number of preview/full points currently
displayed.

Cell Point Gaussian radii are computed server-side as
`display_radius * Point size`; the viewer explicitly enables array scaling, so
the per-stage YAML radii and the side-panel multiplier are both applied.

Install the tested trame packages into a Python environment whose major/minor
version matches ParaView Python, then run:

```sh
python3 -m venv .venv-atcg3d-viewer
.venv-atcg3d-viewer/bin/pip install -r visualization/viewer/requirements.txt
PYTHONPATH="$(pwd)/.venv-atcg3d-viewer/lib/python3.X/site-packages" \
  pvpython visualization/viewer/app.py /path/to/run --port 8080
# If pvpython already sees the installed packages, simply use:
pvpython visualization/viewer/app.py /path/to/run --port 8080
```

The standalone viewer defaults to English. Use `--language zh-CN` for
Simplified Chinese; ATCG3D Studio passes its selected language automatically.

Replace `python3.X` with the exact version reported by
`pvpython -c 'import sys; print(f"python{sys.version_info.major}.{sys.version_info.minor}")'`.

During dragging or playback only preview is loaded. After 250 ms idle, the
newest slider token may load an exactly matching full VTK-HDF keyframe or an
exact checkpoint state reconstructed from its bounded delta chain. Full point data
and vessel centerlines stay server-side in ParaView; the browser receives only
rendered images and UI state. A newly selected time invalidates a stale full
token and its result; the synchronous ParaView `UpdatePipeline()` call itself
is not safely preemptible after it starts.

When Studio is attached, the catalog also merges the one-entry
`live.vtkhdf.series` and `live-vessels.vtkhdf.series` into the archived
timeline. Their files are atomically overwritten rather than accumulated. If
the slider is already at the newest time it follows each replacement; if the
user moved into history, live refresh does not move the selected time.

The dependency-free controller is tested with:

```sh
python3 tests/3d/viewer_controller_test.py
```

Incremental reconstruction requires `h5py` (pinned in `requirements.txt`). A
native, self-contained VTK-HDF export can be created with:

```sh
pvpython visualization/viewer/checkpoint_materializer.py \
  RUN_DIRECTORY CHECKPOINT.h5 output.vtkhdf --compression-level 1
```
