# ATCG3D continuum reduction specification

## State variables

The four cell-number densities are

```text
u_rs(x,t), u_rl(x,t), u_Ks(x,t), u_Kl(x,t),
```

where `s` and `l` mean small and large. Optional ABM ultrasmall cells are
mapped into the small class in v1; supplied reference configurations disable
that stage. The effective occupied fraction is

```text
phi = u_rs + u_Ks + v_l (u_rl + u_Kl),
```

with `v_l=8` in 3D and `v_l=4` in a thin layer. Density fields have units of
cells per ABM voxel measure, so integrating `u` returns cell-equivalent number.

## Population equations

For each field `u_a`,

```text
partial_t u_a = -div(J_a) + R_a(u,N).
```

The finite-volume face flux is the exclusion random-walk reduction

```text
F_a(L,R) = D_a A(phi_face) / h^2
           * [u_a(L) V_a(R) - u_a(R) V_a(L)],
```

where `V_a` is destination vacancy and `A` is the configured crowding factor.
The low-density fixed-26 mapping is

```text
D_a = (9/26) lambda_a,
```

in voxel-squared/hour, with `lambda_a` obtained from the mean configured ABM
migration clock. Large-cell and density-activated-r multipliers are explicit
wrapper parameters. Boundaries are no-flux. No nutrient chemotaxis is included
because the source ABM does not contain that rule.

## Density growth and reactions

Local r/K counts are box averages over the mapped ATCG3D growth-density window.
The same `calculate_density_growth_rate_continuous` function used by the ABM is
evaluated with

```text
r_effective = r_count / M(N),
K_effective = K_count / M(N).
```

Mean division intensity is derived from the ABM work clock. A newly initialized
ABM cycle stores mean work `T_cycle`; density growth depletes that work directly,
so the v2 reaction is

```text
b_i = max(g_i,0) / T_cycle.
```

The legacy v1 reaction retains the earlier normalized-rate closure for exact
reproduction of old continuum runs.

The deterministic reaction operator contains same-size daughters when space is
available, large-to-two-small shape reduction when it is not, density-gated
r-to-K daughter conversion, and the configured delayed-death mean hazard.
For a small cell, the ABM succeeds when at least one neighbour is free. The v2
mean-field closure is therefore `P_success = 1 - phi^m`, with `m=8` in a thin
layer and `m=26` in 3D. Positive reaction increments are capacity limited so
`phi <= phi_max`.

This reaction closure preserves the average rule structure, not the exact
distribution of individual event times or geometric division conflicts.

## Nutrient equation

Both schema versions use the same vascular-surplus convention as
`ATCG3D_Nutrient`. In v2, the consumer fields are cell-number densities,

```text
c_r = u_rs + u_rl,
c_K = u_Ks + u_Kl,
```

and the nutrient equation is

```text
0 = D_N Laplacian(N) - lambda_N N
    - q_r c_r N/(K_N+N)
    - q_K c_K N/(K_N+N)
    + kappa_v V(x,t)(N_v-N).
```

`V` is the fraction of an image voxel occupied by a perfused ABM vessel or a
configured synthetic line. The field uses zero exterior values and fixed-count
relaxed Jacobi iterations. It is refreshed at exact configured times.

The supplied v2 profiles set `q_K=0.010` and `q_r=1.2 q_K=0.012`. Because
`u_l` integrates to large-cell number, large and small cells of the same type
have equal total demand. The legacy v1 equation instead uses volume-weighted
consumer density and remains available for reproducibility.

The capacity multiplier is

```text
M(N) = 1 + (M_max-1)
       * [N/(K_M+N)] / [N_v/(K_M+N_v)].
```

Thus `N=0` reproduces baseline carrying capacity and `M_max=2` matches the old
maximum density relief of `0.5`.

### Transient shared nutrient (schema v3)

Schema v3 instead advances

```text
partial_t N = D_N Laplacian(N) - lambda_N N
              - q C N/(K_N+N),
C = u_rs + u_rl + u_Ks + u_Kl.
```

The supplied v3 profile has `q_r=q_K=q`; large and small fields already
integrate biological cell number, so phenotype and footprint do not change
integrated demand. The initial value is explicit (`N(x,0)=N_v` in v5).
Configured planar-edge and vessel sources obey `N=N_v`; a disabled planar
source leaves a zero-flux outer boundary. Diffusion and decay are explicit and
CFL checked, while the local Michaelis--Menten sink is solved implicitly to
preserve `0 <= N <= N_v`.

Schema v3 fixes `M(N)=1`, uses one common r/K density limit and carrying
capacity, and scales only positive growth by `N/(K_g+N)`. Thus nutrient supply
does not alter hard occupancy capacity. Schemas v1/v2 remain unchanged.

## Conservative ABM coarse-graining

Small agents contribute one cell divided by continuum-voxel measure. A large
agent contributes `1/v_l` at each site in its exact footprint before spatial
binning. Consequently:

```text
integral(u_a) = number of agents in class a,
integral(phi) = ABM occupied biological volume.
```

Perfused vessel sites are binned in the same way into `V`. Importing a
checkpoint first verifies the base ATCG3D checksum.

## Determinism and restart

Traversal order, prefix sums, face fluxes, reaction updates, and nutrient
iterations are fixed. The checkpoint contains all population arrays, nutrient,
vascular fraction, clock, step count, next refresh, solve count, dimensions,
and a dynamics fingerprint. Restart tests compare the full state checksum
across a nutrient refresh boundary.

## Validation contract

The registered continuum test covers strict wrapper loading, exact r/K mass
mapping from the ABM initialization, occupied-fraction bounds, a bounded
single-source nutrient gradient, and bit-identical checkpoint continuation.
The common ATCG3D test suite remains authoritative for the source ABM behavior.
