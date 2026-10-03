# Well-mixed ODE

Build target `atcg3d_ode` loads `ATCG3D_ODE/config/ode_smoke_v1.yaml`.
Its new schema/model is `atcg3d.ode_config` v1 / `well_mixed_structured_v1`.
The `structured_config` reference supplies all base biological parameters and
the v7 common resource/refractory contract. Population entries are small/large
pairs of cell-number densities per unit voxel volume, rather than total cells.

The state contains ordinary r, active r, and K in both stages; active time-mass;
mean nutrient N; optional vascular capacity V; refractory ordinary-r subsets
and their time-mass. Mean clock is time-mass divided by its compartment mass.
The density-growth window counts are the uniform density times its exact edge
volume. The shared density rule sets growth, N/(H+N) scales positive growth,
and vacancy controls division and large-to-small shape reduction. New r
births enter ordinary compartments; configured r-to-K conversion affects
new daughters. Threshold mortality removes all associated clocks in proportion
to their mass.

N follows decay, per-cell Michaelis-Menten consumption and V-weighted perfusion
exchange with the configured vessel value. V is a constant capacity in v1.
It is optional by setting V=0. This spatial mean supply is a volumetric exchange
closure of Dirichlet vessels, not an exact average of a heterogeneous boundary
layer. Activation transfers eligible ordinary mass into a mean-clock cohort;
expiry and refractory release are explicit conservative transitions.

An embedded Dormand-Prince RK45 integrates between transitions. Component-wise
absolute/relative error and positivity reject candidate steps. The solver clips
intervals to pending mean-clock expiries and reports step underflow rather than
silently proceeding. `--checkpoint FILE` and `--resume-checkpoint FILE` preserve
state and time, validate dynamics/solver settings, and support extending the
end time. Checkpoints at output-step boundaries restart bitwise on the same
platform, including a changed configured thread count.

`--model periodic-pde` runs a periodic finite-volume reference on the configured
unit-spaced grid. Its compartments and clocks share conservative fluxes; an
explicit CFL selects substeps. Each grid site uses the same adaptive reaction
kernel and periodic exact-edge density counts. This reference uses isotropic
active diffusion and is intended to verify the well-mixed limit, rather than
reproduce structured direction correlations. Large grids and wide counting
windows are expensive in this reference implementation. Production invasion
remains the structured PDE executable. Both modes emit metrics.csv and a JSON
summary; resumed ODE output uses metrics_resumed.csv.
