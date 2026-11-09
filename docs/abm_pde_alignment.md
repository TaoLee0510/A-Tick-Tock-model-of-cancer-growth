# ABM/PDE alignment and version boundaries

Continuum schema v5 and structured schema v7 introduce the corrected spatial
contract. Continuum v1-v4 and structured v1-v6 retain their published equations
and checkpoint layouts. `atcg3d_legacy_pde_regression_test` compares six pairs
of state checksums with fixtures generated from the phase-0 commit, including
transport, reaction, clock expiry and nutrient updates.

The growth window contains exactly `edge` lattice sites per dimension:
`lower = (edge - 1) / 2`, `upper = edge - lower - 1`. The inclusive interval is
`[anchor - lower, anchor + upper]`. Thus edge 70 covers offsets -34 through 35.
Finite grids clip this interval. At unit spacing, small-cell densities are
literal ABM anchor counts; large-cell densities remain conservative footprint
averages and consequently smooth an anchor at a window boundary.

Continuum v5 subdivides migration until `2*d*D_max*dt/h^2 <= 0.45`, including
the configured active-r multiplier. Ordinary structured diffusion has the
same preflight CFL check; its discrete active transport keeps its independent
move-probability substeps. Unstable ordinary diffusion is rejected before field
allocation. The continuum r200 companion can therefore run independently.

## Shared static sources

Nutrient schema v3 has a required `vascular_geometry` mapping with grid shape,
origin, unit spacing, source mode, synthetic axis/centre/radius and optional
integer `static_sources`. Continuum v5 carries the same settings under `grid`
and `vascular`. Source modes include `synthetic_central_line`, `static_voxels`,
`abm_perfusion`, and `abm_plus_synthetic_line`. Explicit sources supplement
the selected mode. Cylindrical distances are evaluated at voxel centres.

The ABM excludes the shared sources and out-of-grid sites before initialization
and during migration and division. A large cell is admitted only if its whole
footprint is admissible; a thin-layer footprint projects its second z layer
onto the same source mask. ABM imports therefore preserve their initial cell
number in both PDEs. New PDE array initializers reject overlapping populations
instead of silently discarding supplied mass. Legacy structured initializers
keep their original vessel-clearing behavior. Static nutrient sources are
Dirichlet values; the v1/v2 quasi-steady exchange equations remain unchanged.

Use `nutrient_smoke_static_guided_v3.yaml` or
`structured_smoke_2d_256_v7.yaml` as executable examples. Nutrient v3 requires
`halo_voxels >= direction.density_radius`; legacy wrappers retain their original
halo contract. Output-root overrides are described in [output_paths.md](output_paths.md).

## Direction guidance and refractory mass

The new ABM model `low_density_high_resource_bounded_v2` excludes out-of-domain
sites from both density and resource-average denominators. An initial active
direction uses density eligibility and density/resource weights. During
persistence, `persistence_uses_density: true` keeps those gates and weights;
`false` uses resource and distance weights alone. Feasibility and the configured
turn cone always apply. The published `low_density_high_resource_v1` always
uses its original density gates, including during persistence, regardless of
this flag. Structured v5+ requires a resource guidance model; structured v7
also excludes off-grid samples from directional averages.

Structured v7 uses `cohort_clock_refractory_hysteresis_v3`. Each ordinary-r
stage contains an eligible subset and a refractory subset. The refractory
subset is included in ordinary-r totals exactly once. A second field stores
mass times remaining cooldown hours. Both fields undergo the same conservative
ordinary migration operator. Mixing produces a mass-weighted mean cooldown.
Clock expiry of active r adds its mass and cooldown clock to that subset.
The subset becomes eligible when its mean clock reaches zero and local density
is below the off threshold. A different eligible cohort at the same voxel can
activate independently. This mean-clock closure approximates a distribution of
individual cooldowns; it does not resolve separate age bins.

Death scales refractory mass and its clock equally. Newly produced ordinary-r
mass is eligible. Restart format 4 stores both fields, while v5/v6 continue
using the site-local refractory state and restart format 3. The moving-front
refresh preserves its previous update bounds even when the clipped new box is
empty, so previously depleted host voxels are restored.

`atcg3d_alignment_test` checks window endpoints in 2D and 3D, initial counts,
source masks, bounded averages, persistence, r200 conservation, CFL rejection,
cohort transport/mixing/release, thread independence and binary restart.
