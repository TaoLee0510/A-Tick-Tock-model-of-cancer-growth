# ATCG3D ABM-aligned structured PDE

`ATCG3D_StructuredPDE` is the density model for comparisons in which the r-cell
fast-migration rule must retain its ABM memory. It separates ordinary and
activated r density, carries a remaining-activation clock, and transports the
activated density in the same fixed directions used by the ABM.

The legacy profiles below retain their published mean-clock/rate closures.
Schema 7 carries refractory mass with transport, schema 8 adds vascular fields,
schema 9 adds sparse storage, and schema 10 supports a remaining division-work
distribution. Schema 11/12 add activation-duration and speed distributions. Schema 13
records cumulative vascular-deletion mass by phenotype, activity and stage;
its native checkpoint format version is 9. Older schemas retain their state
and output columns unchanged.
Use [the shared contract](shared_rule_contract.md),
[division renewal](renewal_validation.md), and
[activation distributions](activation_distribution.md) for these newer models
and their ensemble verification. All refinements require explicit model choices.

The supplied two-dimensional production profile uses a `2000 x 2000 x 1`
unit-spaced grid and runs to 2160 hours. It is separate from
`ATCG3D_Continuum`, which remains the simpler instantaneous-mobility baseline.

## Shared ABM contract

The aligned base profile is
`configs/atcg2d_legacy_native_r20_aligned_v3.yaml`. It defines activated r
migration as

```text
lambda_active(cell) = 20 * lambda_normal(cell).
```

This is implemented in the ABM itself. Existing configurations that use an
independent activated-rate Beta distribution retain their old behavior.

The structured PDE uses the matching ensemble rate
`20 * E[lambda_normal]`, the ABM `70 x 70` density window, `32 x 32` anchor
blocks, threshold `0.9`, eight in-plane directions (26 in 3D), the 45-degree
initial-direction filter, and the 0.9 direction-persistence probability.
Density only starts the active state. Falling below 0.9 does not turn it off;
the stored remaining clock does.

The v2 profile additionally uses per-cell resource demand: large and small
cells of the same type consume the same total resource, and r consumes 1.2
times K. Activated r selects among density-eligible directions using both lower
density and higher nutrient supply. The low-density gate remains primary.

The v3 profile keeps the v2 resource contract and adds four requested spatial
rules. R-to-K daughter conversion uses the same quantized `70 x 70` density
field as activation; active-r direction scores use 45-degree sectors clipped
by an exact `70 x 70` boundary; crowded small active-r can conservatively
exchange with adjacent small K; and every perfused vessel voxel has zero cell
capacity. V1/v2 remain available unchanged for reproduction of earlier runs.

The v5 profile is isolated in
`ATCG3D_StructuredPDE_NutrientChemotaxis`. It replaces the quasi-steady resource
closure with one transient shared nutrient, fixes all planar edges and vessel
voxels at maximum supply, gives every biological cell identical demand, and
uses common r/K density limits and carrying capacity. Activated r direction
selection depends only on the 70-voxel nutrient field. Clock expiry is followed
by a refractory interval and density hysteresis, so an expired cohort cannot
immediately reactivate. Schemas v1-v4 retain their previous equations.

## Build and run

```sh
cmake -S . -B build-structured -DCMAKE_BUILD_TYPE=Release \
  -DATCG_BUILD_LEGACY_2D=OFF \
  -DATCG_BUILD_3D=ON \
  -DATCG_BUILD_3D_CONTINUUM=ON \
  -DATCG_BUILD_3D_STRUCTURED_PDE=ON
cmake --build build-structured --parallel

build-structured/atcg3d_structured_pde \
  --config build-structured/configs/structured_smoke_r20_resource_guided_v2.yaml

build-structured/atcg3d_structured_pde \
  --config build-structured/configs/structured_legacy_2d_2000_r20_exchange_v3.yaml
```

`structured_benchmark_2d_2000_24h_resource_guided_v2.yaml` exercises the real
production grid and v2 direction rule without writing output. It is intended
for machine-specific runtime and memory checks before the 2160-hour run.

## State and output

For each small/large stage, the model stores ordinary r density, activated r
density by persistent direction, activated remaining-clock mass, and K
density. It also evolves effective nutrient and the vascular source field.

The output directory contains:

```text
metrics.csv                         masses, active fraction, radii, nutrient
fields/field_XXXXXXXX.csv           normal/active r, K, nutrient, vessel
checkpoints/*.structured.bin        sparse exact restart state
config/requested.yaml               structured wrapper
config/continuum_requested.yaml     resolved continuum wrapper
config/effective.json               full resolved configuration
final.json                          final summary and checksum
```

Fields report `r_normal_*` and `r_active_*` separately, as well as total r/K,
activation density, occupancy, nutrient and vessel fraction. Checkpoints store
directional density only inside its active support and accept a changed end
time or output schedule on restart; a changed dynamics fingerprint is rejected.
Metrics also report the assembled r and K consumption rates, making the 1.2
ratio directly auditable.

## Interpretation boundary

This is a deterministic ensemble closure, not a cell-by-cell replay. It aligns
the activation trigger, persistent active state, mean rate multiplier,
direction set/filter and exclusion rule. An ABM checkpoint import preserves
the exact remaining clock and last direction of every imported active cell
after coarse-graining.

The legacy mean-clock model uses the mean of the ABM
`Beta(0.005, 0.011666...) * remaining division time` duration. Consequently,
the model does not reproduce the full duration distribution, UID lineage,
finite-number fluctuations, event ordering, or individual hard-placement
conflicts. Those require ABM ensembles or a PDE-individual hybrid. See
[the model specification](structured_pde_model_spec.md) for the equations and validation
contract.
