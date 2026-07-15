# ATCG 2D to 3D parameter mapping v2

## Scope

This document defines the numerical mapping used by the
`legacy_2d_mapped_v2` profile. The input configuration keeps
the source 2D values and the mapping revision.  Startup resolves those values
once into an effective 3D configuration; simulation hot loops never multiply
an already-resolved value again.

## Unit rules

| Quantity | Mapping |
| --- | --- |
| lattice length, radius, and window edge | unchanged when lattice spacing is unchanged |
| hours, rates per hour | unchanged after confirming rate versus interval semantics |
| probabilities, fractions, and angles | unchanged |
| counts tied to the 6 x 6 growth window | multiplied by the extrusion depth 6 |
| disk and annulus initialization | regenerated as a sphere and spherical shell at the same lattice radii |
| estimator/chunk block edges | numerical controls; selected by convergence and performance tests, not biological scaling |

The density-dependent legacy source has two carrying-capacity branches:

| Ultrasmall | 2D r/K capacity | 2D r/K limit | Effective 3D r/K capacity | Effective 3D r/K limit |
| --- | --- | --- | --- | --- |
| disabled | 31 / 36 | 15.5 / 18 | 186 / 216 | 93 / 108 |
| enabled | 36 / 42 | 18 / 21 | 216 / 252 | 108 / 126 |

The effective values above are final. They must not be multiplied by six a
second time. Accordingly, effective metadata reports
`density_count_scale_2d_to_3d=6` for the documented mapping and the downstream
growth-rule multiplier `carrying_capacity_scale_2d_to_3d=1`; the latter is one
because the stored 93/108 and 186/216 values are already resolved.

Other legacy mappings are:

- growth window: `6 x 6` to `6 x 6 x 6`;
- directional density radius and threshold: `5` and `0.60`;
- persistence and turn half-angle: `0.90` and `45 degrees`;
- migration activation: `70 x 70` to `70 x 70 x 70`, threshold `0.90`;
- normal r migration: `Beta(5,5) * 0.5`; K uses its own base rate;
- finite activation duration: `Beta(0.005, 0.011666666666666667)` times the
  remaining division-cycle time (mean fraction `0.30`);
- r-to-K daughter conversion density: `70 x 70` to `70 x 70 x 70`, with
  the legacy stage-aware capacity divisor, threshold `0.50`, and probability
  `0.05` per successful r-cell division;
- dimensionless interaction defaults: `alpha=2.2`, `beta=0`;
- division timing: base `24 hours`, deterministic fraction `0.9`, stochastic
  fraction `0.1`;
- post-division inherited growth-rate ceilings: r `1.3171805`, K `0.99505180`;
- density-death geometric means: r `48 hours`, K `120 hours`;
- death eligibility threshold: density-adjusted growth rate `0.005`, with the
  legacy positive-support geometric waiting-time model;
- initial geometry source: `R0=60`, `R1=50`, and r fraction `0.5` for the
  density-dependent profile.

Initialization counts are derived from the generated 3D sphere/shell
geometry.  They are not obtained by multiplying the 2D cell count by an
arbitrary constant.

The production profile enables r-to-K conversion. Its 70³ query uses the
incremental anchor block index (`query_block_edge=32`), not a voxel scan. The
smoke profile disables conversion so small deterministic smoke runs do not
silently change their cell-type composition, while retaining all parameters
as explicit schema-v2 fields.

## Values without a 2D biological counterpart

The distance-weight exponent remains zero to preserve equal lattice-direction
weighting.  Density block edge, activation cache block edge, and scheduler
bucket size are engineering controls.  They require estimator-error and
performance tests rather than biological calibration.

Angiogenesis threshold, seed-event rate, vessel speeds, diameter, stopping
lengths, and vascular relief parameters have no legacy 2D source.  A YAML
configuration that enables angiogenesis must state them explicitly.
