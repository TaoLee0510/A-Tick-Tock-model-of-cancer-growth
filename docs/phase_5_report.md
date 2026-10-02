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
These v1 approximations are retained for published runs.

## Distribution and refinement follow-up

Hybrid schema/model v2 transfers actual remaining division work and the
structured duration/rate/direction marginals. Reconstructed cells keep sampled
remaining work. Sub-cell fronts remain at their original sites; residuals after
conversion retain spatial weights subject to vacant volume. This fixes
artificial block-scale tail diffusion caused by repeated exchange.

The new assert-based `atcg3d_hybrid_distribution_test` covers aged-work
conversion, work-moment preservation on immediate round trip, unchanged
sub-cell tails after repeated exchange, 3D large-cell duration/rate transfer,
mixed four-thread bitwise restart and native limiting checksums. Product CLI
tests add v2 dry-run/migration and interval/seed/report dispatch.

`atcg3d_hybrid_validation` runs seven scenarios for each of 16 paired seeds on
a 256-square grid for 48 hours. The regular growing case passes the unchanged
ABM/PDE ensemble tolerances: adaptive/native ABM mass error 0.029583 and radial
L2 0.295661; adaptive/native PDE mass error 0.015309 and radial L2 0.024348.
Both representations remain present, with mean 9 density-to-agent and 256
agent-to-density conversions. Step and exchange refinement checks both pass.
Middle/fine step mass drift is 0.002294 and exchange drift is 0.000417. Reports
retain all raw runs and paired uncertainty; refinement thresholds were fixed
before the ensemble was run.

The refinement verifies reduced interval sensitivity within sampling error,
not a formal convergence-order proof. Joint spatial/work/rate correlations,
reconstructed growth-rate heterogeneity and mean refractory clocks remain
explicit closures. General activated invasive-front calibration remains
separate from this regular growing mixed-case validation.
