# Simulator entry and verification

Build from the repository root with C++20, CMake, yaml-cpp and OpenMP:

```sh
cmake -S . -B build-codex -DCMAKE_BUILD_TYPE=Release -DATCG_BUILD_LEGACY_2D=OFF
cmake --build build-codex --parallel 4
ctest --test-dir build-codex --output-on-failure
```

On macOS, Apple clang uses Homebrew `libomp`. Install `yaml-cpp`, `libomp` and,
for checkpoints/field export, `hdf5`. Add
`-DATCG3D_ENABLE_HDF5_CHECKPOINT=ON` for HDF5 ABM checkpoints and PDE VTK-HDF.
The default build still supports native ODE, PDE and hybrid checkpoints.
The separate `ATCG3D_ENABLE_VTKHDF` option enables the existing ABM point writer
through VTK; PDE ImageData output only needs HDF5.

`atcg_sim --model abm|ode|pde|hybrid --config YAML` dispatches to the model's
native executable without a shell. The executables stay available individually.
`--dry-run` validates configuration; `--output-root` and `ATCG_OUTPUT_ROOT`
mount relative output directories under a local external root. Each supplied
configuration has a distinct default directory. CLI options beyond these depend
on the selected native model; invoke its `--help` for the complete list.

| Model | Configuration | Contract and output |
| --- | --- | --- |
| ABM | Native or nutrient YAML; shared structured YAML | Individual event clocks, directional migration and angiogenesis. Shared mode writes requested JSON reports and paired checkpoints. |
| ODE | `ATCG3D_ODE/config/ode_smoke_v1.yaml` | Well-mixed mean-clock reaction, adaptive RK45, CSV and checkpoint. |
| PDE | Continuum or structured YAML | Density fields, resource/vascular fields, CSV and native checkpoint. |
| Hybrid | `ATCG3D_Hybrid/config/hybrid_smoke_v1.yaml` | Dense PDE core and sparse/front ABM, conservative exchanges, CSV and paired checkpoint. |

Examples:

```sh
build-codex/atcg_sim --model abm --config ATCG3D_SharedRules/config/validation_v7.yaml --report abm_summary.json
build-codex/atcg_sim --model pde --config ATCG3D_SharedRules/config/validation_v7.yaml
build-codex/atcg_sim --model ode --config ATCG3D_ODE/config/ode_smoke_v1.yaml
build-codex/atcg_sim --model hybrid --config ATCG3D_Hybrid/config/hybrid_smoke_v1.yaml --output-root run_outputs
ctest --test-dir build-codex -L validation --output-on-failure
```

The validation label runs 16 paired seeds on a 256-square, 48-hour growth
example and a 48-square, eight-hour vascular example. Reports include the
predeclared tolerances and sampling uncertainty. Unit tests separately cover
published checksum fixtures, exact growth windows, high-rate transport,
nonnegativity, conservation, ODE/PDE agreement and thread-independent restart.
GitHub Actions configures Linux and macOS HDF5 builds, runs all CTests and
uploads both validation reports. The workflow is checked in; a hosted run is
only confirmed once GitHub executes it.

Version and approximation boundaries are explicit in
[shared_rule_contract.md](shared_rule_contract.md),
[abm_pde_alignment.md](abm_pde_alignment.md),
[ode_model.md](ode_model.md), [angiogenesis_fields.md](angiogenesis_fields.md)
and [hybrid_model.md](hybrid_model.md). Aggregate agreement on the included
smokes is not a biological calibration. Published schemas keep their previous
arithmetic; new rules and storage use new schema/model identifiers.

Use [configuration_migration.md](configuration_migration.md) to create a new
configuration graph, [pde_fields.md](pde_fields.md) for viewer/Studio field
output, and [sparse_pde_storage.md](sparse_pde_storage.md) for 10000-square runs.
Model manuals, specifications, viewer and Studio instructions live in `docs/`;
old module README paths are navigation pointers. `cmake --install build-codex
--prefix install-root` installs the executables, YAML graph and documentation.
