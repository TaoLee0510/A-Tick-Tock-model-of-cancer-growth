# Transported division work and regular-cycle validation

Structured schema 10 adds `division_clock.model:
transported_shifted_geometric_v1`. Published schemas and the default
`mean_rate_v1` model retain their arithmetic. The new model follows the ABM
division law: minimum work is `minimum_fraction * base_cycle_hours`, followed
by a geometric number of work quanta. A quantum is inherent growth times
`stochastic_time_quantum_hours`; its success probability is the quantum
divided by the larger of the configured tail work and the quantum.

Four sparse distributions carry remaining work for small/large r and K mass.
ABM initialization imports the actual remaining work. Diffusion, directed
migration and swaps transport the source distribution with accepted mass.
Positive growth advances the work coordinate; completed mass supplies the
division event rate, with newborns and completed parents receiving the same
renewal law. Failed small-K division retains the configured retry work.
Division also resets the completed active-r parent fraction to ordinary r.

`work_bin_width` controls interpolation error. `maximum_work` must cover the
geometric tail to a residual probability below 1e-9; configuration validation
rejects inadequate support. The supplied example uses width 0.5 and maximum
128. Work and activity are separate local marginal distributions. Newborn
work uses the phenotype mean inherent rate, so inherited growth-rate
heterogeneity and age/direction correlations remain closure approximations.
Native checkpoint format 6 includes the work distributions and checksum.
ODE v1 and hybrid v1 reject this model because their current state cannot
carry these distributions.

Run the unchanged, regular-cycle biological parameters with:

```sh
python3 scripts/validate_abm_vs_pde.py --exe build-codex/atcg3d_shared_abm --config ATCG3D_SharedRules/config/regular_cycle_v10.yaml --seeds 16 --output build-codex/validation-regular
```

The 256-square, 48-hour ensemble passes the previously declared tolerances:

| Metric | ABM mean | PDE mean | Error |
| --- | ---: | ---: | ---: |
| Total cell mass | 468.375 | 489.728 | 0.04559 |
| r/K ratio | 1.42338 | 1.32850 | 0.06666 |
| Active fraction | 0 | 0 | 0 |
| r50 | 14.4375 | 13.1875 | 0.08658 |
| r90 | 22.375 | 19 | 0.15084 |
| r99 | 26.1875 | 22.8125 | 0.12888 |
| Normalized radial density L2 | | | 0.27459 |

Scalar errors use the ABM mean as denominator, except absolute active-fraction
error. The original mass/ratio tolerances are 0.35, radius tolerances
0.25/0.25/0.30, and radial L2 tolerance 0.35. Reports retain the paired sampling
allowance and every seed. This regular-cycle case does not validate nonzero
activation.

`atcg3d_division_renewal_test` checks the literal shifted-geometric law,
conservative transport, serialization, the independent renewal equation
for expected pure-birth population, and native restart across thread counts.
`atcg3d_regular_cycle_validation` registers the ensemble under the
`validation` CTest label.

Release checks pass: default 38/38 CTests and HDF5 40/40 CTests, including
the regular-cycle ensemble and all published checksum fixtures.
