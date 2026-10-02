# ATCG3D structured PDE nutrient chemotaxis

This directory contains v5-v7 configurations for `atcg3d_structured_pde`.
The `continuum_*.yaml` files are companion continuum configurations referenced
by the `structured_*.yaml` entry points. Structured v5 uses continuum schema
v3; structured v6 uses continuum schema v4. Older structured schemas retain
their v1-v4 equations.

The v5 model has a transient shared nutrient field. Its configured planar edges
and perfused/synthetic vessel voxels are fixed sources at `N_max`. All cells have the
same integrated per-cell demand, independent of phenotype and footprint.
Nutrient scales division and survival but never changes the hard occupied-area
limit. Activated r transport uses only the 70-voxel directional nutrient field;
cell density is not a direction gate. The activation clock still expires, and
a site-local cooldown plus density hysteresis prevents immediate reactivation.
This refractory state is attached to the grid site in both v5 and v6; it does
not follow ordinary-r mass during migration.

V7 references continuum schema v5 and replaces the site-local cooldown with an
ordinary-r refractory mass subset and its transported mass-weighted cooldown
clock. It also uses the exact ABM growth-window endpoints and the shared static
source mask during ABM initialization. Published v5/v6 runs remain unchanged.
See [the alignment contract](../docs/abm_pde_alignment.md) for version boundaries
and the mean-clock closure, and use `config/structured_smoke_2d_256_v7.yaml`
for the new 48-hour smoke profile.

V6 replaces planar-edge supply with a moving tumour-front boundary. A smoothed,
thresholded occupied field selects the largest connected tumour component and
fills its holes. Host voxels outside that component are clamped to `N_max`;
the configured vessel source remains fixed inside the tumour. The v6
equivalence configurations compare reference and optimized active-region
execution; the model rules and checkpoint state are the same in both modes.

The transient initial condition is explicit: `N(x, 0) = N_max` throughout the
plane. Consumption then depletes the interior while planar edges and vessel
voxels remain fixed at `N_max`. This avoids interpreting an initially empty
nutrient field as an acute starvation experiment.

Configurations:

- `config/structured_smoke_2d_256_v5.yaml`: 24-hour implementation smoke test.
- `config/structured_control_edges_2d_256_24h_v5.yaml` and
  `config/structured_control_vessels_2d_256_24h_v5.yaml`: matched source-geometry
  controls.
- `config/structured_calibration_2d_512_360h_v5.yaml`: first calibration run.
- `config/structured_2d_2000_720h_v5.yaml`: completed production calibration
  after the 512-grid acceptance checks passed.
- `config/structured_2d_2000_resume_720_to_2160h_v5.yaml`: exact continuation
  from the accepted 720-hour checkpoint to the completed 2160-hour result.
- `config/structured_smoke_2d_256_24h_v6.yaml`: moving-front smoke run.
- `config/structured_calibration_2d_512_360h_v6.yaml`: moving-front calibration.
- `config/structured_2d_2000_2160h_v6.yaml` and
  `config/structured_2d_10000_2160h_v6.yaml`: larger moving-front configurations.
- `config/structured_equivalence_2d_256_24h_v6_reference.yaml` and
  `config/structured_equivalence_2d_256_24h_v6_optimized.yaml`: deterministic
  active-region equivalence configurations.

Each entry point has its own relative output directory. See
[`output path overrides`](../docs/output_paths.md) for external-volume runs.

Completed numerical checks and the 512-grid results are recorded in
[`CALIBRATION.md`](CALIBRATION.md).

The old v4 1440-hour checkpoint and partial 1584-hour metrics are intentionally
left in their original run directory.

## Governing nutrient equation

For nutrient `N` and total biological-cell density `C`, v5 advances

```text
∂N/∂t = D_N ΔN - λN - q C N/(H_N + N).
```

Here `C` is the sum of r and K cell numbers, not occupied volume. Consequently
a small cell and a large cell each contribute exactly one unit to `C`, and r
and K share the same `q`. The initial condition is `N(x,0)=N_max`; the four
planar edges and vessel voxels obey `N=N_max` for all later times. The
Michaelis--Menten sink is solved locally and implicitly after an explicit,
CFL-checked diffusion step, so `0 <= N <= N_max` is preserved.

With the supplied calibration values (`D_N=0.1`, `q=0.01`, `H_N=0.25`), the
reference Damkohler number at length 35 and one cell per voxel is
`q L²/[D_N(H_N+N_max)] = 98`. It is therefore a consumption-dominated regime,
as requested; this dimensionless check is more meaningful than directly
comparing a diffusion coefficient with a consumption rate.

If a source-control configuration disables planar Dirichlet supply, the
remaining outer boundary is zero-flux. The production configuration fixes both
planar edges and vessels at `N_max`.

## Cell equations and constraints

Normal r and K use the existing conservative diffusion/reaction closure. Both
types now share `common_density_limit=16.75` and
`common_carrying_capacity=33.5`; nutrient never multiplies the hard capacity.
It only scales positive division by

```text
S(N) = N/(growth_half_saturation + N).
```

Activated r remains a finite set of directional density fields in v5 and v6. For direction
`d`, the 70-voxel sector mean nutrient `N_bar[d]` gives the jump weight

```text
w[d] = |d|^(-distance_exponent)
       exp(chemotaxis_strength * (N_bar[d]/N_max - N_local/N_max)).
```

No density value appears in this direction score and no high-density sector is
discarded. A move is still rejected at the outer domain or into a vessel voxel.
The conservative active-r/K exchange closure handles occupied target sites.

The 70-voxel density field is retained only for event rules: it triggers r
activation at 0.90 and r-to-K conversion at the configured threshold. Each
activated cohort carries a remaining-time clock. On expiry it returns to normal
r, disarms the local site for a 24-hour refractory interval, and can re-arm only
when the same 70-voxel density falls to 0.80 or below.

## Acceptance checks

The automated suite checks fixed edge/vessel nutrient values, transient inward
propagation, equal integrated demand across phenotype and cell size,
nutrient-directed movement even when every sector exceeds the former 0.60
density threshold, activation expiry/refractory behavior, vessel exclusion,
and deterministic checkpoint/restart including the new refractory state.
