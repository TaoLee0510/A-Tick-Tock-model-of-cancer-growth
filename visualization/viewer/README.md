# ATCG3D trame viewer

This is a thin client over a ParaView backend. It supports arbitrary 3D camera
rotation/pan/zoom, time slider, play/pause, time and cell count, `cell_type`
coloring, Point Gaussian display, radius scaling, preview/full status, and live
series refresh.

Install trame packages into the Python environment used by ParaView, then run:

```sh
pvpython visualization/viewer/app.py /path/to/run --port 8080
```

During dragging or playback only preview is loaded. After 250 ms idle, the
newest slider token may load an exactly matching full frame. Full point data
stays server-side in ParaView; the browser receives rendered images and UI
state. The dependency-free controller is tested with:

```sh
python3 tests/3d/viewer_controller_test.py
```
