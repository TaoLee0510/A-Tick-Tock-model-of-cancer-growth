# Phase 1 verification

Continuum v5, nutrient v3 and structured v7 add exact growth-window endpoints,
continuum active-r migration substeps, shared static/synthetic vessel geometry,
bounded guidance averages and transported refractory mass/clock fields.
Published continuum v1-v4 and structured v1-v6 retain their arithmetic and
checkpoint formats. The ordinary structured diffusion CFL is checked before
allocation, v5+ guidance is validated, and empty moving-front updates preserve
the previous nutrient update bounds.

New CTest targets are `atcg3d_alignment_test` and
`atcg3d_legacy_pde_regression_test`. The latter fixtures were generated from
phase-0 commit `e9c3782`, compiled independently, and compare all six published
structured schemas and their continuum companions after eight transport and
reaction steps. The alignment target covers 2D/3D edge windows, initial mass,
source masks, bounded resource averages, persistence, r200 conservation,
ordinary CFL rejection, cohort activation/expiry/transport/mixing/hysteresis,
four-thread restart and HDF5 ABM restart with static exclusion. A mortality
restart test covers canonical reconstruction of rounded active-direction totals;
this fixes a cache discrepancy discovered by the full 48-hour restart experiment.

Release validation: 28/28 CTests pass with the default build; 29/29 pass with
HDF5 checkpoints enabled. Both new 256-square configurations complete 48 hours
through their command-line executables. Structured restart from the 24-hour
checkpoint reaches the same final checksum, `2845072932649579256`, as the
uninterrupted run. Continuum r200 completes independently with checksum
`3218589719640665088`. The nutrient-v3 example passes command-line dry-run.

Remaining modeling limits: large-cell field imports distribute their mass over
footprints, so window-boundary counts are smoothed relative to literal ABM
anchors. V7 refractory clock mixing is a mean-clock closure rather than a full
age distribution. Ensemble ABM/PDE quantitative agreement and the transient
ABM resource contract remain phase-2 work. The historical diagnostic artifacts
remain outside the tracked build dependency set.
