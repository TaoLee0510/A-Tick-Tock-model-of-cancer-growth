# Hybrid model

`atcg3d_hybrid --config ATCG3D_Hybrid/config/hybrid_smoke_v1.yaml` loads
hybrid schema 1, `hybrid_shared_grid_v1`, and a structured v7/v8 rules file.
The new wrapper uses the native sparse ABM event engine and the native
structured PDE, including its nutrient, VEGF, tip and vessel equations.
`--mode all_abm|all_pde` delegates to the corresponding native engine; their
cell-state checksums match the corresponding standalone run exactly.

Adaptive mode classifies smoothed occupancy with on/off hysteresis. Dense
locations carry PDE populations; agents represent sparse locations and the
front. At each macro step the PDE advances first, including nutrient and
vascular fields. The ABM then runs to the same time with those fields held
fixed. Agent occupancy blocks PDE transport/growth; density occupancy blocks
agent placement, migration and division. Both populations contribute to the
same growth/activation counts, nutrient consumption, VEGF production and
moving-front mask. New perfused vessels displace either representation.

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
