# Configuration migration

`atcg_config_migrate` copies a configuration and its referenced base,
continuum and structured YAML into a new directory. It rewrites references and
output directories to relative names, validates the resulting graph with the
native loaders, and writes `MIGRATION.md`. Source files remain intact; existing
output directories are refused. After failure, a draft graph may remain for
inspection; choose a new directory for the next attempt.

```sh
build-codex/atcg_config_migrate --input ATCG3D_SharedRules/config/validation_v7.yaml --output-directory migrated_run
build-codex/atcg_sim --model pde --config migrated_run/model.yaml --dry-run
```

Supported upgrades are ABM v1-v3 to v3, nutrient v1-v3 to v3, continuum v3-v6
to v6, and shared-resource structured v5-v15 to v15. ODE v1 and hybrid v1/v2/v3
wrappers are copied with their references upgraded. Hybrid v1 retains its
mean-clock closure; v2 retains its transported-work selection. Structured v1-v4 and continuum v1-v2 lack
shared-resource calibration parameters; the tool rejects them rather than
inventing scientific parameters.

Structured v5/v6 upgrades select transported cohort refractory clocks;
new continuum versions use the exact ABM edge window. These changes have new
fingerprints. Migration creates a fresh run and clears old checkpoint imports;
it does not convert a biological checkpoint or promise identical new trajectories.
The numeric time interval and initial parameters remain visible for review.
Schema v13 retains the mean-rate division closure unless a configuration
explicitly selects transported shifted-geometric division work. Migration
preserves that explicit selection and its work-grid parameters. Duration and
rate distributions also require explicit selections; upgrading a schema alone
keeps its previous mean-clock/rate closures.

Older nutrient configurations need `--grid-edge N` to define the new finite
resource domain. Nutrient v1 also needs `--nutrient-K-per-cell-hour RATE`: its
old rates used occupied voxels and the new contract uses cells. The caller
selects the new K rate; the tool carries forward the old r/K rate ratio and a
common saturation. Unequal old saturations or nonpositive rates require manual
calibration. The tool expands the halo to at least the direction radius.

```sh
build-codex/atcg_config_migrate --input ATCG3D_Nutrient/config/nutrient_smoke_v1.yaml --output-directory migrated_nutrient --grid-edge 128 --nutrient-K-per-cell-hour 0.01
```

Review the new domain and uptake parameters for the intended experiment before
running it. The explicit rate above is an example choice, not a conversion of
all large/small-cell voxel uptake into identical individual uptake.

Structured v13 adds cumulative vascular-deletion accounting without changing
the vascular deletion operator. Hybrid v3 must be selected explicitly; the
migrator preserves a wrapper's existing model instead of silently changing its
volume, activation-density or footprint-resource coupling.

Structured v14 migration keeps `axial_diffusion_v1` and
`clipped_location_mean_v1` defaults to preserve the numerical operators of a
v13 source. The named `feasible_fixed_lattice_jump_v3` and
`truncated_normal_expectation_v2` corrections are explicit choices in the new
validation configurations. Operator interventions require
`shared_operator_switches_v1`; they are not enabled by migration.

Structured v15 adds the explicit `shared_vegf_lattice_v2` vascular model.
Migration preserves `vegf_tip_density_v1` or `disabled`; selecting v15 alone
does not replace a published vascular mechanism. See the
[shared angiogenesis contract](shared_angiogenesis_contract.md) for the
individual tip law, centerline diagnostics and common perfusion assumptions.
