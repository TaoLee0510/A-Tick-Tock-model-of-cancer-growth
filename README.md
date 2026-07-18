# ![image](Tick-Tock_1.gif) A Tick-Tock model of cancer growth


Cancer growth model with phenotypic trade-offs between proliferation and survival; colonization and competition.

The historical two-dimensional simulator remains available as `atcg`. The
separate `atcg3d` target implements sparse, event-driven three-dimensional
growth with dynamic angiogenesis, HDF5 checkpoints, VTK-HDF time-series output,
and ParaView/trame visualization. All 3D model and run parameters are supplied
through strict schema-v3 YAML; see
[`docs/3d_model_spec.md`](docs/3d_model_spec.md),
[`docs/3d_angiogenesis_spec.md`](docs/3d_angiogenesis_spec.md), and
[`docs/3d_benchmark.md`](docs/3d_benchmark.md).

```sh
cmake -S . -B build-3d -DCMAKE_BUILD_TYPE=Release \
  -DATCG_BUILD_LEGACY_2D=OFF -DATCG_BUILD_3D=ON \
  -DATCG3D_ENABLE_HDF5_CHECKPOINT=ON -DATCG3D_ENABLE_VTKHDF=ON
cmake --build build-3d --parallel
build-3d/atcg3d --config configs/atcg3d_legacy_2d_mapped_v3.yaml --dry-run
```

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
The viewer uses fixed green r cells, red K cells, blue vessel tubes, and offers
whole, cut, and adjustable X/Y/Z slab views with free 3D camera rotation.

Graphical abstract of a Tick-Tock model of cancer growth: 

![image](model.gif)



Population initial growth with different parameters:

![image](visualization.gif)



Single cell initial growth with different parameters:

![image](visualization_SingleCell.gif)



Experimental observation of populationi initial growth:

![image](observation.jpg)
