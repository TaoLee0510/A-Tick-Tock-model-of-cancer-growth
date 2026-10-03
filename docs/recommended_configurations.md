# Recommended configuration registry

Each model has one current starting point below. A recommendation applies to
the stated experiment and evidence, rather than to every tumor, duration or
vascular regime. The graphs were generated with `atcg_config_migrate` and
validated through `atcg_sim --dry-run`. All references and output directories
are relative. Generated output namespaces include the model directory so
ABM and PDE recommendations can run under one external output root.

| Model | Current configuration | Applicable evidence and limits |
| --- | --- | --- |
| ABM | [abm/model.yaml](recommended/abm/model.yaml) | Shared r20, 256-square, 8-hour invasion reference. Phase B uses two independent 16-seed groups. Native individual clocks remain the spatial reference. |
| Structured PDE | [pde/model.yaml](recommended/pde/model.yaml) | Same r20 experiment and explicit feasible jumps/true normal mean. Phase B TOST passes; this does not establish long-time or r200 equivalence. |
| Hybrid | [hybrid/model.yaml](recommended/hybrid/model.yaml) | Version 5 resources with the v4 front/core policy, r retained as agents, K core conversion. Phase D compares invasive mixed representation and four splitting/exchange refinements. See its report for the final statistical result. |
| ODE | [ode/model.yaml](recommended/ode/model.yaml) | Version 1 well-mixed mean-clock reaction and adaptive RK45. Uniform periodic PDE agreement, conservation and nonnegativity tests apply; no spatial-front equivalence is claimed. |

The ABM and PDE graphs select structured schema 17 with
`published_sector_sums_v1`. Schema migration preserves their explicitly
validated numerical operators; it does not silently enable the new FFT sector
cache. The hybrid graph explicitly selects v5 resource initialization and full
workspace persistence. Its unchanged front policy is from v4. ODE retains
its mean-clock closure. Native checkpoint fingerprints include graph versions,
so these migrated recommendations start new runs.

`configuration_registry.json` lists every repository YAML configuration and
its status. The four top-level inputs above are `recommended`; their copied
references are `recommended_dependency`. Named production and 3D inputs are
`benchmark`: they define workloads and remain conditional on recorded
performance, boundary and restart results. Prior statistical experiments,
smokes, benchmarks and calibration inputs are `reproduction`. Published
native base profiles and prior hybrid policies are `legacy` unless they are
dependencies of a reproduction. These categories do not change loader
behavior or delete any input. Temporary YAML under build directories is not
part of the registry.

Deprecation requires a new recommendation, a migration note that identifies
changed rules and evidence, and at least one release retaining the old loader
and reproduction fixtures. Published arithmetic and checkpoint readers stay
available. A schema/model may be removed only in an explicitly announced
breaking release with its reproducible source revision documented. A failed
TOST or boundary-contact run remains a failed reproduction result.

The production ABM/PDE/hybrid graphs use 2000-square grids, r200, transient
nutrient and shared VEGF angiogenesis for a configured 2160 hours. They are
workload definitions; short performance prefixes do not make them the current
scientific recommendation. See [the production protocol](production_benchmark_contract.md)
and [the 3D sector protocol](three_dimensional_sector_contract.md).

Regenerate a graph in a new directory, then review and dry-run it:

```sh
build-codex/atcg_config_migrate \
  --input ATCG3D_SharedRules/config/active_r20_v14.yaml \
  --output-directory recommended_pde
build-codex/atcg_sim --model pde --config recommended_pde/model.yaml --dry-run
```
