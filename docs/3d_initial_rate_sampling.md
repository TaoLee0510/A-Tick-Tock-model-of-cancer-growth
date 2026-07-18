# 3D initial growth and migration rate sampling

## Mapping rule

Initial per-cell kinetic rates are biological parameters, not spatial counts.
Unlike the 6x6 to 6x6x6 density carrying-capacity mapping, they are not
multiplied by six when the lattice becomes three-dimensional. The production
profile keeps the legacy growth laws and the same activated-r beta family, but
uses the new configurable multiplier default of 1 instead of the historical
hard-coded 200. Each draw is deterministic for an immutable cell UID.

The production growth model is `legacy_truncated_normal_v1`:

| Type | Mean | Standard deviation | Minimum | Maximum |
|---|---:|---:|---:|---:|
| r | 1.1832 | 0.2441 | 1.0722619 | 1.3171805 |
| K | 0.6832 | 0.3764 | 0.33963482 | 0.99505180 |

Sampling uses inverse-CDF sampling over the configured truncated interval. The
values in the `mean` column are parameters of the untruncated normal, so the
mean of generated values is generally different after truncation.

The r-cell inherent migration rate used during density activation is configured
under `migration.activated_r_rate`, independently of initialization. Its
schema-v3 and supplied-profile default scale is `1.0`:

| Type | Alpha | Beta | Scale | Post-processing |
|---|---:|---:|---:|---|
| r (activated, schema and supplied-profile default) | 0.01 | 0.0566666667 | 1 | values `<=0.5` become `0.25` |
| K | 5 | 5 | 0.25 | none |

The historical hard-coded r multiplier was 200. It remains reproducible by
setting `migration.activated_r_rate.scale: 200` in a YAML file, but no supplied
schema-v3 profile uses it by default.

The beta sampler uses two deterministic gamma variates and evaluates their
ratio in log space. Log-space evaluation is necessary for the very small r-cell
shape parameters and avoids underflow-driven bias. K keeps its initialization
model because it does not use the r-cell density-activation state machine.

## YAML schema

The activated r distribution is nested below `migration`; only K migration is
nested below `initial`. Models are discriminated and unknown, unused, missing,
non-finite, or inconsistent fields cause configuration loading to fail. Schema
v2 files are rejected rather than silently reinterpreting the old initial-r
field as the new activated-r field.

```yaml
migration:
  activated_r_rate:
    model: beta
    alpha: 0.01
    beta: 0.0566666667
    scale: 1.0
    lower_clamp:
      enabled: true
      threshold: 0.5
      value: 0.25

initial:
  growth_rate:
    model: legacy_truncated_normal_v1
    truncated_normal:
      r:
        mean: 1.1832
        standard_deviation: 0.2441
        minimum: 1.0722619
        maximum: 1.3171805
      K:
        mean: 0.6832
        standard_deviation: 0.3764
        minimum: 0.33963482
        maximum: 0.99505180
  migration_rate:
    K:
      model: legacy_beta_v1
      beta:
        alpha: 5.0
        beta: 5.0
        scale: 0.25
        lower_clamp:
          enabled: false
```

## Reproducibility contract

`sample_initial_cell_rates()` and its component functions are pure. The initial
r migration field is sampled from `migration.activated_r_rate`; subsequent
committed division cycles use the same configured law with the division event
sequence. Draws use separate stateless RNG domains, so calling growth before
migration, retrying initialization, changing thread count, or enabling output
does not alter a cell's sampled rates.

The same sampler initializes a daughter that converts from r to K after a
successful division. That K draw is keyed by the immutable daughter UID and is
therefore independent of retry count and thread order. The sampled inherent
growth rate is then limited by
`division.timing.K_max_inherent_growth_rate`; an unconverted daughter inherits
the divided parent rate subject to the corresponding r or K ceiling.

The effective configuration and dynamics JSON include the activated-r model,
shape, scale, and clamp. Checkpoint/config compatibility checks therefore
reject any change that would alter initialization or future division cycles.
