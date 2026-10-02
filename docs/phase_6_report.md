# Phase 6 verification

Added `atcg_sim` dispatch for ABM/ODE/PDE/hybrid, an installed executable/config
layout, and `atcg_config_migrate` for supported version graphs. Migration
rewrites relative references, validates native loaders, records rule changes
and requires an explicit per-cell uptake rate for old nutrient v1 inputs.
Hybrid CLI now emits CSV metrics under the shared local output-root override.

PDE VTK-HDF ImageData output streams rows into compressed fields and appends
an atomic relative series. The existing viewer/Studio viewer discovers PDE
series, renders surfaces and selects fields. Substantive model, viewer and
Studio manuals/specifications now live in `docs/`; old paths are link stubs.
The GitHub Actions workflow builds/tests HDF5 on Linux and macOS and publishes
both small ensemble reports. Hosted execution is not claimed by local checks.
Installation and five installed configuration graphs pass native dry-run
validation, including the nutrient-chemotaxis configuration references.

Structured v9 sparse zero-page storage preserves native arithmetic. New
assert-based tests cover dense/sparse per-field bitwise equivalence, four-thread
restart, full zero remapping, active-region budget rejection, CLI routing,
configuration graph migration and HDF5 field dimensions/values. Existing viewer
controller tests add PDE timeline discovery; concurrent control tests isolate
their request directories by process. Hybrid CLI tests run the complete smoke
and verify metrics under the output-root override. The opt-in 10000-square
benchmark loads the supplied configuration, runs 48 hours, computes final
diagnostics/checksum and measures peak resident memory: 4,573,052,928 bytes
(4.57 GB). Final mass is 668.387; checksum is 12240940797807697984.
Default 36/36 and HDF5 38/38 CTests pass, including all legacy checksums and
both 16-seed validation targets. ParaView 6.1.1 independently reads 2D/3D files;
the viewer's offscreen density/nutrient rendering and field selection pass.

Limits are explicit: sparse v1 covers thin-layer static-vascular sparse support,
not arbitrarily filled grids or dynamic angiogenesis at 10000-square scale.
Early PDE schemas without shared-resource calibration cannot be upgraded
without parameter choices. Hybrid cohort/individual conversion and the ODE
mean-clock closure retain the scientific approximations listed in their model
manuals. Quantitative smoke agreement is not comprehensive biological calibration.
