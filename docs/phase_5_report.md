# Phase 5 verification

Added `ATCG3D_Hybrid`, `atcg3d_hybrid` and hybrid schema/model v1. Native ABM
and structured PDE share nutrient, VEGF, vessels, density/activation counts
and occupancy. Lie splitting advances the PDE before each ABM interval;
fixed exchanges use conservative regional dependent rounding. Checkpoints
retain both engines, clocks, core classification, lineage and ABM vasculature.
Published-model hooks are dormant by default and legacy checksum fixtures pass.

New assert-based `atcg3d_hybrid_test` checks frozen-biology mass conservation,
agent-to-density and density-to-agent conversion, fractional refractory mass,
3D large-cell activation clocks, four-thread bitwise mixed restart, native
ABM/PDE limiting checksums, coupled vascular restart, and all-ABM hypoxic
vascular checkpoint continuation without requiring HDF5.
Default 34/34 and HDF5 35/35 CTests pass, including both 16-seed validation
suites. The 48-square coupled example completes eight hours with conversion
in both directions. Full results and implementation limits are described in
`hybrid_model.md`: restored individuals receive fresh division/rate samples,
cohort clocks are mean clocks, and exchange redistributes within 16-voxel blocks.
Hybrid ensemble calibration and convergence studies remain future scientific
validation, beyond conservation/restart/limiting-engine correctness checks.
