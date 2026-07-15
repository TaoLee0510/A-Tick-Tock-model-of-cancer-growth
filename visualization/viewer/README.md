# ATCG3D trame viewer

This is a thin client over a ParaView backend. It supports arbitrary 3D camera
rotation/pan/zoom, time slider, play/pause, time and cell count, `cell_type`
coloring, Point Gaussian display, radius scaling, preview/full status, and live
series refresh. A synchronized server-side Tube overlay renders
`vessels.vtkhdf.series`; the UI can color it by perfusion or inward/outward role
and adjust its radius scale. The toolbar distinguishes the true live-cell count
stored in frame FieldData from the number of preview/full points currently
displayed.

Cell Point Gaussian radii are computed server-side as
`display_radius * Point size`; the viewer explicitly enables array scaling, so
the per-stage YAML radii and the toolbar multiplier are both applied.

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

Replace `python3.X` with the exact version reported by
`pvpython -c 'import sys; print(f"python{sys.version_info.major}.{sys.version_info.minor}")'`.

During dragging or playback only preview is loaded. After 250 ms idle, the
newest slider token may load an exactly matching full frame. Full point data
and vessel centerlines stay server-side in ParaView; the browser receives only
rendered images and UI state. A newly selected time invalidates a stale full
token and its result; the synchronous ParaView `UpdatePipeline()` call itself
is not safely preemptible after it starts. The dependency-free controller is
tested with:

```sh
python3 tests/3d/viewer_controller_test.py
```
