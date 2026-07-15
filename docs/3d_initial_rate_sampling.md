# 3D initial growth and migration rate sampling

## Mapping rule

Initial per-cell kinetic rates are biological parameters, not spatial counts.
Unlike the 6x6 to 6x6x6 density carrying-capacity mapping, they are not
multiplied by six when the lattice becomes three-dimensional. The production
profile preserves the 2D probability laws and assigns one deterministic sample
to each immutable cell UID.

The production growth model is `legacy_truncated_normal_v1`:

| Type | Mean | Standard deviation | Minimum | Maximum |
|---|---:|---:|---:|---:|
| r | 1.1832 | 0.2441 | 1.0722619 | 1.3171805 |
| K | 0.6832 | 0.3764 | 0.33963482 | 0.99505180 |

Sampling uses inverse-CDF sampling over the configured truncated interval. The
values in the `mean` column are parameters of the untruncated normal, so the
mean of generated values is generally different after truncation.

The production migration model is `legacy_beta_v1`:

| Type | Alpha | Beta | Scale | Post-processing |
|---|---:|---:|---:|---|
| r | 0.01 | 0.0566666667 | 200 | values `<=0.5` become `0.25` |
| K | 5 | 5 | 0.25 | none |

The beta sampler uses two deterministic gamma variates and evaluates their
ratio in log space. Log-space evaluation is necessary for the very small r-cell
shape parameters and avoids underflow-driven bias. The smoke profile uses the
`fixed` model for both rates.

## YAML schema

The rate configuration is nested below `initial`. Models are discriminated:
`fixed` accepts only `fixed`, while the legacy models accept only their named
parameter block. Unknown, unused, missing, non-finite, or inconsistent fields
cause configuration loading to fail.

```yaml
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
    model: legacy_beta_v1
    beta:
      r:
        alpha: 0.01
        beta: 0.0566666667
        scale: 200.0
        lower_clamp:
          enabled: true
          threshold: 0.5
          value: 0.25
      K:
        alpha: 5.0
        beta: 5.0
        scale: 0.25
        lower_clamp:
          enabled: false
```

## Reproducibility contract

`sample_initial_cell_rates()` and its component functions are pure. Each draw
is keyed by `(seed, uid, rate kind, cell type)` using separate stateless RNG
domains. Calling growth before migration, migration before growth, retrying an
initialization, changing thread count, or enabling output does not alter a
cell's sampled rates. These functions do not read or increment the cell event
sequence used for migration, division, death, or angiogenesis.

The same sampler initializes a daughter that converts from r to K after a
successful division. That K draw is keyed by the immutable daughter UID and is
therefore independent of retry count and thread order. The sampled inherent
growth rate is then limited by
`division.timing.K_max_inherent_growth_rate`; an unconverted daughter inherits
the divided parent rate subject to the corresponding r or K ceiling.

The effective configuration JSON includes every model and distribution
parameter so checkpoint/config compatibility checks detect any change that
would alter future initialization.
