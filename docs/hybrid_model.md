# Hybrid model

The growing-cycle configuration
`ATCG3D_Hybrid/config/hybrid_regular_cycle_v2.yaml` selects schema 2,
`hybrid_distributions_v2`. It carries the structured PDE's transported
division-work distributions through both interface directions. Structured
v11/v12 activation-duration and ordinary-rate distributions are also supported.
Schema/model v1 retains the closure described below.

`atcg3d_hybrid --config ATCG3D_Hybrid/config/hybrid_smoke_v1.yaml` loads
hybrid schema 1, `hybrid_shared_grid_v1`, and a structured v7/v8 rules file.
The new wrapper uses the native sparse ABM event engine and the native
structured PDE, including its nutrient, VEGF, tip and vessel equations.
`--mode all_abm|all_pde` delegates to the corresponding native engine; their
cell-state checksums match the corresponding standalone run exactly.

Adaptive mode classifies smoothed occupancy with on/off hysteresis. Dense
locations carry PDE populations; agents represent locations below the
smoothed-occupancy thresholds. Classification currently uses occupancy alone,
so front distance, active-r presence and nutrient gradient are not yet
representation guarantees. At each macro step the PDE advances first, including nutrient and
vascular fields. The ABM then runs to the same time with those fields held
fixed. Agent occupancy blocks PDE transport/growth; density occupancy blocks
agent placement, migration and division. Both populations contribute to the
same growth/activation counts, nutrient consumption, VEGF production and
moving-front mask. New perfused vessels delete overlapping cells in either representation.
They do not conservatively relocate those cells.

The fixed exchange interval must be an integer multiple of the PDE step.
Agent-to-density conversion transfers one cell, stage, active direction and
remaining activation/refractory time. Large cells spread one unit over their
4/8-voxel footprint. Density-to-agent conversion uses deterministic dependent
rounding within fixed 16-voxel exchange blocks. Each phenotype/stage/activity
budget pays exactly one unit per agent; a seed/exchange-index permutation
chooses sites. Separate ordinary-r and refractory-r budgets prevent cooling
mass from being rounded into armed cells. Fractional remainders stay in PDE
compartments, sharing available volume. This bounds redistribution distance
and conserves total cell mass, without creating a cell from a fractional budget.

`--checkpoint FILE` writes a manifest plus `.pde.bin` (except all-ABM) and
`.resource.bin`. Keep the set together. The manifest is committed last and
existing checkpoint names are refused. `--resume-checkpoint FILE` restores
cell-slot/free-list layout, UID/event counters, lineage, full ABM vasculature,
PDE buckets/clocks/resources, the core mask, refractory map and exchange
counters. Derived coupling arrays are reconstructed. All modes support
checkpoints without HDF5. Binary v1 requires a compatible C++ ABI; this is not
an architecture-independent interchange format.

Limitations of the new mixed closure: reconstituted agents draw a fresh
individual division cycle and migration-rate heterogeneity. PDE-to-ABM active
cohorts retain a mass-weighted remaining clock and restart direction history
in the stay bucket; ordinary-r refractory cohorts retain their mean clock.
Regional conversion redistributes residual spatial density inside its block.
Ultrasmall colocations remain ABM because the structured model has two stages.
These are approximation choices in hybrid v1, not changes to published models.
Mixed active mass retains the native PDE float precision; integer exchange
and normal/K mass budgets use doubles. No ensemble calibration of hybrid
front statistics is claimed by the unit tests.

## Distribution interface v2

An agent entering the core deposits its actual remaining division work after
the lazy growth-clock update, remaining activation time, ordinary migration
rate and direction. A large agent spreads those distributions over its 4/8
footprint with the same cell-mass weights. A departing density cohort supplies
mass-weighted work/time/rate/direction marginals. Each new UID samples these
marginals with independent deterministic seed/UID/exchange-index domains;
its division cycle is retained instead of being redrawn. The regional budget
debits an exact unit and proportionally debits its distribution mixtures.
An expired sampled activation time creates an ordinary refractory agent.

Sub-cell regions that cannot produce an individual retain their fields and
distributions exactly. After a conversion, residual mass is restored in
proportion to each cohort's original spatial density, subject to free volume;
bounded water filling redistributes only when a new agent occupies a site.
This prevents repeated exchange from concentrating fractional tails at a
permuted site and creating block-scale numerical diffusion.

The v2 manifest uses a distinct native format identifier. Checkpoints retain
the work/time/rate banks, reconstructed agents and exchange counters. The
same ABI limitation as v1 applies. Native `all_abm` and `all_pde` modes still
delegate exactly to their standalone engine, including state checksums.

The interface retains marginal distributions in expectation, not their joint
correlations. Newly reconstructed individual inherent growth rates are sampled
from the configured law. Ordinary refractory cohorts still use a weighted
remaining clock. Fractional activity fields retain native float roundoff.
Those approximations are explicit and independent of the exact interface mass
budget. Ultrasmall colocations remain individual agents.

## Ensemble and interval refinement

```sh
python3 scripts/validate_hybrid.py --exe build-codex/atcg3d_hybrid --seeds 16 --output hybrid-validation
ctest --test-dir build-codex -R atcg3d_hybrid_validation --output-on-failure
```

The validation uses 16 paired seeds, a 256-square grid and 48 hours of regular
shifted-geometric growth. Seven runs per seed compare native ABM/PDE limits
and adaptive hybrid; split steps are 0.25, 0.125 and 0.0625 hours, and exchange
intervals are 2, 1 and 0.5 hours. Every adaptive ensemble must contain both
representations and conversions in both directions.

Native comparisons use the existing ABM/PDE tolerances: mass and r/K ratio
0.35, active fraction 0.05 absolute, radius quantiles 0.25/0.25/0.30 and
normalized radial-profile L2 0.35. Refinement uses scalar relative tolerance
0.10, active fraction 0.02 absolute and profile L2 0.15. Middle/fine drift must
pass these fixed tolerances and decrease relative to coarse/fine drift within
paired sampling uncertainty; profile error uses delete-one jackknife
uncertainty. Raw runs, uncertainty, coverage and pass/fail are retained in JSON
and Markdown. This is a growing mixed-case verification, not calibration of
arbitrary invasive fronts or proof of asymptotic convergence order.

## Volume coupling v3

`ATCG3D_Hybrid/config/hybrid_regular_cycle_v3.yaml` selects
`hybrid_volume_coupling_v3` and structured v13. Published v1/v2 arithmetic,
fingerprints, placement thresholds and checkpoint identifiers remain intact.
The v3 manifest has a separate native format identifier.

PDE anchor counts contribute to a large agent's activation density as
`count / (edge^d / large_cell_volume)`, matching the native ABM denominator.
Growth and activation counts remain anchor counts. Nutrient consumption and
VEGF production use a separate consumer field: each large agent contributes
one cell divided equally over its four planar or eight spatial voxels.
Total demand is one biological cell. This is the same footprint weighting as
the shared ABM environment. Occupancy remains one unit per footprint voxel.

An entering footprint voxel must satisfy
`PDE volume + live agent volume + 1 <= maximum_occupied_fraction`.
`minimum_density` no longer defines an empty destination. The native thin-layer
ABM stores the upper half of a large footprint at z=1; v3 checks that storage
against the corresponding planar PDE voxel at z=0. The candidate agent is
checked only at entering voxels during migration, as in the native grid.

Structured v13 records deleted normal r, active r and K mass by small/large
stage in metrics and checkpoint state. Adaptive v3 also records agents deleted
by the shared vessel exclusion rule in those counters. This loss accounting
is cumulative and is not a conservative vessel-displacement algorithm.

The v3 ensemble repeats the existing seven native/mixed/refinement scenarios
with 16 seeds and the exact previously declared tolerances. It remains a smoke
comparison until the separate statistical-equivalence work is completed.
