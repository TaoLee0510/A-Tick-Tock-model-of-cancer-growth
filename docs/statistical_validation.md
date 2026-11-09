# Statistical validation protocol

This protocol is declared before running the phase B ensembles. The old
`error <= tolerance + sampling_allowance` reports remain smoke checks. They
cannot establish equivalence. New statistical reports use the following
unchanged, previously published practical equivalence margins:

| Quantity | Margin | Scale and rationale |
| --- | ---: | --- |
| Total cell-number mass | 0.35 | Relative to the ABM mean; existing closure tolerance |
| r/K ratio | 0.35 | Relative to the ABM mean; existing phenotype closure tolerance |
| Active fraction of r | 0.05 | Absolute fraction; existing activation tolerance |
| r50 and r90 | 0.25 | Relative to the ABM mean; existing spatial closure tolerance |
| r99 | 0.30 | Relative to the ABM mean; existing tail closure tolerance |
| Radial profile L2 | 0.35 | Existing normalized, measure-weighted profile tolerance |
| Vascular length and perfused volume | 0.75 | Existing exploratory vascular smoke tolerance |
| Lesion perfused fraction | 0.15 | Absolute fraction; existing vascular smoke tolerance |

These margins express the original numerical closure acceptance budget, not
biological or clinical interchangeability. In particular, the broad vascular
margins do not validate the mechanisms that phase C will align. A failed
statistical result must be reported as inconclusive or outside the margin;
neither margins nor sampling uncertainty may be added to the equivalence bound.
Command-line tolerance overrides may only tighten these margins.

## Independent realization and baseline

For N >= 16, ABM group A uses seeds 1 through N, ABM group B uses seeds N+1
through 2N, and PDE uses seeds 1 through N. The A/B seed sets are disjoint.
The baseline is the signed difference of the two ABM ensemble means, with
paired standard error and RMS realization difference also reported. The
PDE-minus-ABM error is divided by the absolute baseline mean error. A zero
baseline yields a null ratio with an explicit reason, not an arbitrary floor
or an infinite JSON value. The RMS-normalized difference remains informative
when cancellation makes the baseline mean small. Seed pairing is declared in
advance; it does not imply that the two model RNG trajectories remain coupled.

## Equivalence tests

At alpha = 0.05, scalar endpoints use two one-sided paired Student t tests
(TOST). Absolute fraction contrasts use P-A with bounds +/- margin. Relative
contrasts test mean(P-(1-margin)A) > 0 and mean(P-(1+margin)A) < 0. This
linear formulation includes uncertainty in the ABM reference mean. The
two-sided paired t statistic and p value for P-A = 0 are descriptive only:
a small difference can be significant and still equivalent. Each one-sided
p value must be strictly below 0.05. Zero standard error is handled by the
strict deterministic inequalities; equality at a bound fails.
Relative margins at a zero reference mean are undefined. Such contrasts are
reported with null errors/p values and an explicit reason, and fail the global
equivalence claim. No absolute floor or extra tolerance is introduced to make
an uncovered vascular quantity pass.

The profile uses four-voxel annuli in 2D, or spherical shells in 3D, and the
published measure-weighted relative L2 distance between normalized ensemble
profiles. A delete-one-pair jackknife supplies its standard error. TOST bounds
are +/- the existing L2 margin, applied to the estimated nonnegative distance.
This nonlinear jackknife t approximation is explicitly approximate at small N
and near zero distance. Identical paired profiles pass deterministically;
heterogeneous paired profiles are not treated as zero uncertainty. Profile
intervals and all delete-one estimates are retained for review.
For the descriptive paired profile t statistic, the contrast is each seed's
L2(P,A) minus L2(B,A). It tests excess profile discrepancy over the independent
ABM variability baseline, rather than a signed field difference. This is an
exploratory diagnostic and does not affect the global-distance TOST criterion.

Student t probabilities are evaluated with a regularized incomplete beta
function; no normal approximation or fixed t(15) critical value is used for
the new tests. The implementation is checked against analytic df=1 and df=2
distributions and the [NIST t table](https://www.itl.nist.gov/div898/handbook/eda/section3/eda3672.htm).
The TOST interpretation follows [Lakens (2017)](https://doi.org/10.1177/1948550617697177):
equivalence at alpha 0.05 corresponds to containment of a 90% confidence
interval within prespecified equivalence bounds, rather than failure to reject
a zero difference.

## Time series, boundaries and attribution

Sampling is fixed before a run: every four hours by default, with the initial
state and final endpoint retained. Intervals from 4 through 24 hours are
accepted. Runs shorter than four hours are explicitly endpoint-only statistical
diagnostics, without a resolved time-series claim.
Every sampled time compares mass, r/K, active fraction, all three radius
quantiles and the profile. The report retains every time table and maximum
absolute error. Global trajectory equivalence requires all prespecified
metric/time tests to pass (an intersection-union test); repeated samples from
one seed are not independent replicates. No multiplicity adjustment is needed
for this joint claim, and passing selected times cannot imply global success.

Reports include grid shape, radius to the nearest domain face, r99/half-width
and cell mass at the outer grid faces. r99 >= 0.8 of the minimum planar
half-width is flagged as boundary influenced. The existing 128-square, two-hour
r200 reproduction is retained with that warning; it is not free-space evidence.
Long validation rejects boundary-influenced runs instead of enlarging margins.
Nonzero mass at an outer face is separately flagged, including numerical
diffusion tails. The declared interior acceptance guard uses r99/half-width,
not an arbitrary cutoff on trace amounts at the face.
The published shared ABM base uses an expandable sparse domain; its resource
grid does not itself confine individual positions. The PDE has reflecting
population transport at its grid faces. Boundary-exposed reproductions thus
also compare these different domain closures. The radius guard limits this
effect in accepted ensembles; it does not prove exactly zero boundary flux.

Mechanism attribution uses matched-seed one-at-a-time interventions for
migration, activation, division, crowding exchange and boundary geometry.
Intervention effects and their ABM/PDE interaction contrasts are reported;
effects are not assumed additive. Boundary changes must preserve the initial
physical population and sources. A numerical correction is only accepted on
a new schema/model, with its mechanism tests and old-schema bitwise checks.

## Opt-in runs

Long ensembles are enabled by `ATCG_ENABLE_LONG_VALIDATION=ON`, registered
with the `validation_long` label and excluded from default CI. The declared
cases are 512-square/360 hours and 2000-square/720 hours for r20 and r200.
The r200 cases must satisfy the radius/boundary guard; otherwise they fail and
require a larger domain. Registering these cases is not evidence that they
have completed. Reports always distinguish executed runs from unrun workloads.
The r200 grids are enlarged to 4096-square in both long configurations; their
filenames retain the requested nominal workload for comparison. The actual
shape is always written in realization reports. The existing one-million
active-voxel budget remains a hard limit, so these opt-in workloads may also
fail that budget before phase E adds another storage implementation.

## Version 14 corrections and interventions

Schema 14 adds named normal transport and growth closures. The published
axial-diffusion and clipped-location-mean closures remain available and are
unchanged in schemas 1 through 13. `truncated_normal_expectation_v2` uses the
mean of the actual truncated ABM growth distribution, rather than its location
parameter. Integrating x times the normal density gives
`mean + sigma * (phi(lower)-phi(upper))/(Phi(upper)-Phi(lower))`.
This fixes the parameter-to-mean mapping; it does not represent cell-specific
inherited growth rates or their correlations with age and phenotype.

`feasible_fixed_lattice_jump_v3` includes all eight planar or 26 spatial unit
lattice directions. It also models the ABM choice among feasible directions,
which the old vacancy-throttled axial flux does not. For independent direction
feasibility probabilities a_d, the unconditional chance to choose d is
`a_d * integral_0^1 product_{e!=d}(1-a_e+a_e*t) dt`. Their sum is
`1-product_d(1-a_d)`: only a cell with no feasible direction stays. Large cells
use the number of newly entered footprint sites (volume minus footprint
intersection), rather than testing the overlapping footprint again.
Feasibility still assumes independent destination occupancy; overlapping
footprints and quenched per-cell migration rates remain closure limitations.
The model requires unit spacing and zero distance/crowding exponents. An
intermediate `fixed_lattice_jump_v2` retains unconditional vacancy throttling
but corrects the stencil alone, allowing the two effects to be distinguished.

`operators.model: shared_operator_switches_v1` allows one of migration,
activation, division or exchange to be disabled without changing other rate
parameters. No-migration retains active-clock expiry but disables displacement
and swapping. No-division disables birth, shape reduction and division-triggered
conversion while retaining survival and resource consumption. ABM lazy division
work can accumulate, but division events never commit; PDE division work does
not advance. No-activation starts and remains inactive in the fresh-run
intervention. No-exchange disables both native swapping mechanisms. Enlarging
both physical grid extents around the same origin-centered initialization is
the boundary intervention; it preserves sources and population, rather than
switching nutrient boundary physics. These are interventions, not independent
additive terms in a fitted decomposition. All interaction contrasts, including
null effects in uncovered regimes, must be retained.

The no-migration intervention exposed an additional division closure error:
ordinary PDE births were added at the mother's occupied voxel, whereas a small
ABM daughter is placed in a free adjacent voxel. Schema 14's optional
`small_daughter_placement: feasible_neighbor_birth_v2` estimates neighboring
feasibility before reaction, retains the mother at its original location and
places prospective daughters using the same uniform-feasible kernel. Newborns
receive fresh work and ordinary-r state. A deterministic recipient-volume cap
can reject competing prospective births; it never deletes existing mass.
Large-cell division and reduction retain their published local closure.
Independent neighbor occupancy and simultaneous birth competition remain
approximations. This model is tested separately from the preceding v14
transport/growth corrections and uses the same predeclared margins.

## Reproducing the reports

The default build registers two statistical ensembles alongside the eight
historical smoke suites. The validation label runs both types; each report
states its acceptance mode explicitly.

```sh
ctest --test-dir build-codex --output-on-failure -L '^validation$'
python3 scripts/validate_abm_vs_pde.py --exe build-codex/atcg3d_shared_abm \
  --config ATCG3D_SharedRules/config/regular_cycle_spatial_birth_v14.yaml \
  --seeds 16 --require-interior --output build-codex/validation-statistical-spatial-birth
python3 scripts/ablate_abm_vs_pde.py --exe build-codex/atcg3d_shared_abm \
  --config ATCG3D_SharedRules/config/regular_cycle_v14.yaml \
  --output build-codex/ablation-regular
python3 scripts/ablate_abm_vs_pde.py --exe build-codex/atcg3d_shared_abm \
  --config ATCG3D_SharedRules/config/active_r20_v14.yaml \
  --output build-codex/ablation-active-r20
python3 scripts/ablate_abm_vs_pde.py --exe build-codex/atcg3d_shared_abm \
  --config ATCG3D_SharedRules/config/regular_cycle_spatial_birth_v14.yaml \
  --interventions full no_migration --output build-codex/ablation-spatial-birth
```

After completing the regular full reference, `--correction-study` reuses that
reference and runs the three declared closure reversions. `--analysis-only`
recomputes statistics from complete existing seed/time-series reports and
rejects missing or misaligned records. Ablation command success means the
diagnostics completed; its JSON reports equivalence separately for every
intervention. The 48-hour `activation_ablation_v14.yaml` is an additional
prepared workload, not an executed result in the phase B report.

```sh
cmake -S . -B build-long -DCMAKE_BUILD_TYPE=Release \
  -DATCG_BUILD_LEGACY_2D=OFF -DATCG_ENABLE_LONG_VALIDATION=ON
cmake --build build-long -j 4
ctest --test-dir build-long --output-on-failure -L validation_long
```
