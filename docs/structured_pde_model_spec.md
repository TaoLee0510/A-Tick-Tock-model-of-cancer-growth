# Structured migration PDE specification

## Population state

For stage `j` (small or large), the state is

```text
n_j(x,t)       ordinary r density
a_j,d(x,t)     activated r density with direction d
q_j,d(x,t)     remaining activation-time mass
k_j(x,t)       K density
N(x,t)         effective nutrient
V(x)           perfused vessel fraction
```

Direction bucket `d=0` means that the cohort has not yet selected a persistent
direction. The other buckets correspond exactly to the eight in-plane or 26
three-dimensional ABM offsets. Total r density and occupied fraction are

```text
r_j = n_j + sum_d a_j,d
phi = r_small + k_small + v_large (r_large + k_large),
```

where `v_large=4` in 2D and `8` in 3D.

`q/a` is the cohort's mean remaining active time. The same accepted transport
fraction is applied to `a` and `q`, so migration cannot alter that clock.

## Activation and return to baseline

The trigger field reproduces the ABM block-anchor query. Total cells are
summed in a `70 x 70` window centred on each quantized `32 x 32` query block.
For a normal-r cohort at a site whose stage-specific density is at least 0.9,

```text
n_j -> 0
a_j,0 -> a_j,0 + n_j
q_j,0 -> q_j,0 + tau_bar n_j.
```

For newly activated continuum mass,

```text
tau_bar = E[Beta(0.005, 0.011666...)]
          * T_base / E[growth_r]
        = 0.30 * T_base / E[growth_r].
```

The density field is not consulted to deactivate a cohort. During each active
transport substep, `q` decreases by `a Delta t`; when the remaining clock is
exhausted, `a` is returned to `n`. This is the required hysteresis missing from
the instantaneous-mobility continuum baseline.

When initialized from an ABM checkpoint, the imported contribution to `q`
uses each cell's exact scheduled activation end time. New activations after the
PDE starts use the mean-duration closure above, rather than the full Beta clock
coordinate.

In v5, clock expiry also disarms the local activation state for 24 hours. It
can re-arm only after the cooldown reaches zero and the same 70-by-70 density
field is at most 0.80. A later rise to the 0.90 on-threshold may start a new
episode; the expired episode itself cannot immediately restart.

## Directional active transport

The aligned ABM configuration defines

```text
lambda_active(cell) = 20 lambda_normal(cell).
```

The continuum event intensity is therefore

```text
Lambda_active = 20 E[lambda_normal].
```

For the v2 profile, candidate directions use the legacy radius-five cone. In
v3, each 45-degree forward sector is clipped by an exact `70 x 70` boundary,
so direction choice and activation observe the same spatial scale. Let `rho_d`
and `N_d` be the mean cell density and normalized nutrient in that sector.
Candidate weights are

```text
w_d = |d|^-eta
      [epsilon + (1-epsilon)(1-rho_d/0.6)]^alpha_rho
      [epsilon + (1-epsilon)N_d]^alpha_N.
```

The density gate is applied before nutrient weighting, so high nutrient cannot
pull activated r back into a direction whose density exceeds 0.6. With no
previous direction, mass is split in proportion to `w_d`. With a previous
direction, the 0.9 forward/0.1 turn prior is multiplied by these weights and
renormalized. Uniform density and nutrient therefore recover the original
persistence probabilities. The nutrient-coupled ABM samples the same law.

An attempted directional jump is accepted with destination volume-filling
factor

```text
A_j(x+d) = (1 - phi(x+d)/phi_max)^v_j,
```

with `v_j=1` for small cells and `v_large` for large cells. Rejected mass stays
in its original direction bucket. Boundaries are no-flux. The explicit jump
operator is substepped so its event probability does not exceed the configured
0.15 limit.

Ordinary r and K retain the symmetric exclusion-diffusion closure used by
`ATCG3D_Continuum`. Thus only activated r has a persistent velocity state.

In v5, active-r direction selection contains no density gate or density score.
For the normalized local nutrient `N_0` and mean nutrient `N_d` in each exact
70-voxel directional sector,

```text
w_d = |d|^-eta exp(chi (N_d-N_0)).
```

Every geometrically feasible, nonvascular direction remains available even
when all sectors exceed density 0.60. Density is retained only for activation
and r-to-K event gates. The fixed-direction persistence and conservative
exchange operators are unchanged.

### Crowding exchange and vessel exclusion (v3)

When every nonvascular neighbour of a small active-r cohort is occupied, v3
allows a swap with adjacent small K. For an accepted exchange amount `S`,

```text
a_small(x) -= S       K_small(x) += S
a_small(y) += S       K_small(y) -= S.
```

The active clock moves with `a_small`. Target demand is limited by available K
mass, so the exchange is nonnegative and conserves r mass, K mass, and occupied
fraction at both endpoints. Large cohorts do not exchange. This is a
deterministic mean-field closure of the ABM singleton transaction; it does not
reproduce individual wait/retry conflicts.

For v3, `V(x)>0` implies zero cell capacity. Initialization removes any cell
mass overlapping a supplied or synthetic vessel, and all migration, exchange,
and reaction operators preserve zero cell density on those voxels. The vessel
continues to act as a nutrient source.

## Reaction and nutrient operators

Growth uses the source ABM's continuous density-growth function over the
configured local box. In v4, division intensity is `max(g,0)/T_cycle`, matching
the ABM work clock: the inherent rate is already present in `g` and is not
divided out a second time. For a small cell, successful placement means that at
least one possible neighbour is free, represented by
`P_success = 1 - phi^m` with `m=8` in a thin layer and `m=26` in 3D. Earlier
profiles retain their original reaction closure for reproducibility.

The reaction closure includes small/large stages, volume-limited daughter
placement, large-to-small shape reduction, delayed death mean hazards and
density-gated r-to-K daughter conversion. New r daughter mass begins in the
ordinary compartment.

In v3, the conversion gate is stage-specific and uses the same quantized
`70 x 70`/`32 x 32` anchor-box density field as migration activation, with the
configured conversion threshold `0.5`. V1/v2 retain their local-voxel gate for
reproducibility.

Cell occupancy remains volume-weighted, but the v2 nutrient sink is counted by
cell number:

```text
c_r = r_small + r_large
c_K = K_small + K_large
q_K = 0.010
q_r = 1.2 q_K = 0.012.
```

Nutrient solves

```text
0 = D_N Laplacian(N) - lambda_N N
    - q_r c_r N/(K_N+N)
    - q_K c_K N/(K_N+N)
    + kappa_v V (N_v-N)
```

by deterministic relaxed Jacobi iteration at exact refresh times. Nutrient
changes carrying capacity through the same saturating multiplier used by the
continuum baseline.

Thus one large and one small cell of the same type have equal integrated
resource demand. The large-cell sink is spatially distributed over its
footprint; its total is not multiplied by four in 2D or eight in 3D. Legacy v1
profiles retain their original occupied-voxel sink for reproducibility.

V5 selects continuum nutrient schema v3: `q_r=q_K`, the four planar edges and
vessel voxels are fixed at `N_max`, and `N(x,0)=N_max`. Both phenotypes use the
same density limit and carrying capacity. Nutrient multiplies positive growth
by `N/(K_g+N)` but never changes the hard occupied-area limit. Vessel voxels
remain cell-free.

## Numerics, restart and validation

Population, activation-density and growth-density work is restricted to exact
support bounds; this changes traversal cost, not the equations. Direction
transition sets are precomputed. Values below `minimum_density` are pruned as
the declared finite-support tolerance.

The version-2 binary checkpoint stores clock state, step counters, nutrient
schedule, population support, directional active support, resource fields and
a dynamics fingerprint. Directional arrays are sparse by bounding box. V5 uses
checkpoint version 3, which additionally stores cooldown and re-arm fields.

The registered test verifies:

- the ABM assigns every initial r cell `active_rate = 20 * normal_rate`;
- a dilute active cohort stays active until its stored clock expires;
- the quantized 70-by-70 density trigger activates ordinary r;
- the ABM directional filter creates a positive outward centroid shift;
- equal-density resource gradients bias activated r toward higher nutrient;
- a resource source 30 voxels away biases v3 transport through its 70x70 sector;
- r-to-K conversion uses the quantized 70x70 density rather than voxel occupancy;
- packed small active-r and K exchange conservatively;
- vessel voxels remain cell-free after initialization and transport;
- small and large cells have equal v2 integrated consumption and r/K is 1.2;
- sparse checkpoint continuation is bit-identical;
- occupied fraction remains within capacity.

These are rule-level and deterministic-restart checks. Quantitative agreement
with the stochastic ABM should be evaluated against multi-seed ensemble
statistics (mass, active fraction, radial quantiles and spatial density), not
against one ABM realization.
