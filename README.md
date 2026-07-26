# ![image](Tick-Tock_1.gif) A Tick-Tock model of cancer growth


Cancer growth model with phenotypic trade-offs between proliferation and survival; colonization and competition.

The historical two-dimensional simulator remains available as `atcg`. The
separate `atcg3d` target implements sparse, event-driven three-dimensional
growth with dynamic angiogenesis, HDF5 checkpoints, VTK-HDF time-series output,
and ParaView/trame visualization. Runtime density queries use exact incremental
70³ activation-block counts and per-cell 6³ growth-window counts. Asynchronous
VTK-HDF/checkpoint output writes directly from a frozen snapshot without
rebuilding a second simulation. All 3D model and run parameters are supplied
through strict schema-v3 YAML; see
[`docs/3d_model_spec.md`](docs/3d_model_spec.md),
[`docs/3d_angiogenesis_spec.md`](docs/3d_angiogenesis_spec.md), and
[`docs/3d_benchmark.md`](docs/3d_benchmark.md).

The production YAML uses `output.storage.mode: journal_delta_hdf5_v2`.
Hourly previews stay bounded, compressed, and self-contained; full VTK-HDF is
a configurable 24-hour series; hourly exact restart state is stored as an
HDF5 v6 base followed by schema-v8 stable-slot mutation journals, with a
weekly base by default. An hourly checkpoint therefore copies only changed
slots, exact free-list operations, global/vascular state, and new lineage
edges instead of freezing every live cell. All three sampling intervals remain
ordinary YAML values and may be set below one hour.

```sh
cmake -S . -B build-3d -DCMAKE_BUILD_TYPE=Release \
  -DATCG_BUILD_LEGACY_2D=OFF -DATCG_BUILD_3D=ON \
  -DATCG3D_ENABLE_HDF5_CHECKPOINT=ON -DATCG3D_ENABLE_VTKHDF=ON
cmake --build build-3d --parallel
build-3d/atcg3d --config configs/atcg3d_legacy_2d_mapped_v3.yaml --dry-run
```

For a long run that must survive the terminal or Codex task ending, use a
detached `screen` session. The foreground managed wrapper records the actual
simulation PID, forwards session termination, and atomically writes its
eventual exit code. On the older macOS `/usr/bin/screen`, configure its PTY log
in a small screenrc; `logfile flush 1` publishes output at one-second cadence:

```sh
CONTROL=/path/to/control/atcg3d_run
SCREENRC=/path/to/control/atcg3d.screenrc
printf '%s\n' \
  'deflog on' \
  'logfile /path/to/control/atcg3d_run.log' \
  'logfile flush 1' > "$SCREENRC"

/usr/bin/screen -c "$SCREENRC" -dmS atcg3d_run \
  /absolute/path/scripts/run_atcg3d_managed_foreground.sh \
  /absolute/path/build-3d/atcg3d \
  /absolute/path/config.yaml \
  "$CONTROL"

cat "${CONTROL}.pid"
tail -f "${CONTROL}.log"
# Present only after clean completion, failure, or a managed termination:
cat "${CONTROL}.exit"
```

Use `screen -ls` to verify the detached session. To stop it cleanly, send
`TERM` to the PID recorded in `${CONTROL}.pid`; the wrapper records exit code
143. A plain background `&` process is not sufficient in managed shells.

## ATCG3D Studio macOS app

`ATCG3D Studio` is the first native macOS control application. Its right-side
panel loads and edits the complete YAML tree, shows the C++-validated effective
configuration, starts the simulator through the detached launcher, publishes
pause/resume/checkpoint-stop/terminate requests, displays progress and ETA, and
opens an existing run directory. The embedded 3D area runs the ParaView/trame
backend, so full point data remains outside the web view.

Build the ad-hoc-signed local application with:

```sh
cd visualization/studio/desktop/src-tauri
cargo test --release
cargo tauri build --bundles app
open "target/release/bundle/macos/ATCG3D Studio.app"
```

The resulting application is
`visualization/studio/desktop/src-tauri/target/release/bundle/macos/ATCG3D Studio.app`.
It is locally ad-hoc signed but not Apple-notarized. ParaView 6.1.1 and the
Python/trame environment below remain external runtimes and their paths are
editable in the Run panel. Closing Studio stops its viewer backend but not a
detached simulation.

Studio defaults to English and provides an English/Simplified Chinese selector
in the top bar. The choice is persisted locally and also applies to the
embedded 3D viewer.

## ATCG3D trame viewer setup on macOS

ParaView 6.1.1 uses Python 3.12.7, so the viewer environment must use the same
Python major/minor version. Create a persistent pyenv environment and install
the repository-pinned viewer dependencies with:

```sh
# Skip this command when `pyenv versions` already lists Python 3.12.7.
pyenv install 3.12.7

pyenv virtualenv 3.12.7 atcg3d-paraview-3.12.7

$HOME/.pyenv/versions/atcg3d-paraview-3.12.7/bin/python -m pip install \
  -r visualization/viewer/requirements.txt
```

Activate the environment for ordinary Python commands with:

```sh
pyenv activate atcg3d-paraview-3.12.7
```

The viewer itself must run with ParaView's `pvpython`. Expose the pyenv
packages through `PYTHONPATH`, then provide an ATCG3D output run directory:

```sh
PYTHONPATH="$HOME/.pyenv/versions/atcg3d-paraview-3.12.7/lib/python3.12/site-packages" \
/Applications/ParaView-6.1.1.app/Contents/bin/pvpython \
visualization/viewer/app.py /path/to/run --port 8080
```

Open `http://localhost:8080` in a browser. The run directory must contain the
generated `preview.vtkhdf.series`, `full.vtkhdf.series`, and
`vessels.vtkhdf.series` catalogs. More viewer details are in
[`visualization/viewer/README.md`](visualization/viewer/README.md).
The viewer uses fixed green r cells, red K cells, and blue vessel tubes, with
independent visibility switches. It offers whole, cut, adjustable X/Y/Z slab,
and a frozen current-view-aligned plane; the camera remains freely rotatable.
An optional translucent overlay displays the configured vascular density-relief
range from `run.json`.

At non-keyframe hours the viewer reconstructs the matching checkpoint only in
the ParaView backend after the 250 ms idle debounce. For a standalone file that
native ParaView can open directly:

```sh
PYTHONPATH="$HOME/.pyenv/versions/atcg3d-paraview-3.12.7/lib/python3.12/site-packages" \
/Applications/ParaView-6.1.1.app/Contents/bin/pvpython \
visualization/viewer/checkpoint_materializer.py \
/path/to/run /path/to/run/checkpoints/checkpoint_....h5 output.vtkhdf
```

Graphical abstract of a Tick-Tock model of cancer growth: 

![image](model.gif)



Population initial growth with different parameters:

![image](visualization.gif)



Single cell initial growth with different parameters:

![image](visualization_SingleCell.gif)



Experimental observation of populationi initial growth:

![image](observation.jpg)
