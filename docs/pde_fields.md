# PDE VTK-HDF fields

HDF5 builds accept `--vtkhdf-fields` on continuum/structured executables to
write fields beside existing CSV snapshots. `--field-vtkhdf FILE` writes one
final field independently of configured periodic output. Both options work
through `atcg_sim --model pde`. Existing output is refused.

```sh
build-codex/atcg_sim --model pde --config ATCG3D_SharedRules/config/angiogenesis_v8.yaml --vtkhdf-fields --output-root field_runs
```

Files use VTK-HDF ImageData v2 with sample coordinates at PDE voxel centres.
PointData contains normal/active/refractory r by stage, K, total r/K, nutrient,
occupied volume fraction and perfused vessel fraction. Vascular models also
include VEGF and vessel tips. Continuum files contain its four population
fields. Scalar dataset axes run z/y/x with singleton axes omitted. Attributes
carry extent, origin, spacing and identity direction; FieldData carries time.
The layout follows the [VTK-HDF specification](https://docs.vtk.org/en/v9.6.0/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html).

The writer streams rows through compressed HDF5 chunks, closes a temporary file
and then renames it. `fields.vtkhdf.series` is committed after each file and
uses relative paths. Resume appends to that series. File enumeration counts CSV
snapshots when field export is enabled, so paired files do not skip indices.

The existing viewer/Studio embedded viewer discovers this field timeline when
an ABM point timeline is absent. It renders density surfaces and provides a
field selector for nutrient, VEGF and the other arrays. The count label becomes
Grid samples; point radii and r/K individual visibility switches are hidden.
World/view clipping and the time slider remain available. Open the PDE run
directory in Studio's viewer panel or run the viewer as described in
[viewer.md](viewer.md). For 3D fields, clipping exposes the interior.

Validation includes direct HDF5 shape/value assertions, independent ParaView
6.1.1 reading of 2D and 3D fixtures, and offscreen viewer rendering/field
selection. Rendering QA uses `pvpython --no-mpi --force-offscreen-rendering`;
the ordinary viewer uses its configured ParaView runtime.
