# Shared ABM and structured PDE rules

`atcg3d_shared_abm` accepts the same structured YAML as the PDE.
Use `--model abm` or `--model pde`, `--seed N`, `--threads N` and
`--report summary.json`. The ABM adapter selects the new
`nutrient_gradient_shared_resource_v3` model and a finite unit-spaced domain.
Existing ABM models keep their original rules. Sources and initial vessel
exclusion use the shared static geometry introduced in phase 1.

The resource solve advances at the continuum time interval with the same
explicit diffusion, exponential decay, positive implicit Michaelis-Menten
uptake, vessel Dirichlet values and optional moving tumour front. Each cell
consumes once, distributed across its footprint. Both phenotypes use the
configured common density limit and carrying capacity. Positive growth is
multiplied by N/(H+N); negative density growth retains the existing death rule.

Active-r directions use exp(chi*(sector_mean_N-local_N)), with exponents
clipped to [-40,40]. The 70-square sector uses offsets -34 through 35 and the
configured cone angle. Only finite-domain sites enter its average. Two-dimensional
row spans and nutrient row prefixes avoid repeated square scans; three-dimensional
sector offsets are precomputed once. Guidance applies to persistent choices and
swaps as well as initial direction selection. For this model
`persistence_uses_density` does not add a density gate: vacancy and vessel
exclusion determine legal destinations, and nutrient determines their weights.

Each ABM UID carries its own refractory deadline and hysteresis arm. Expiry
starts cooldown; the cell re-arms only after the deadline and a density at or
below the off threshold. Its activation duration distribution scales with its
inherent full cycle. PDE v7 transports refractory mass and its mean clock.
These are different distribution closures and are compared statistically.

The new `shared_fixed_lattice_means_v2` diffusion mapping uses D=3*lambda/8
in a thin layer and D=9*lambda/26 in 3D, matching the coordinate second moment
of eight or 26 equiprobable jumps. The validation configuration uses unit large
mobility and no additional crowding power; target vacancy already suppresses
entry. The old mapping remains unchanged.

HDF5 builds support `--checkpoint state.h5` and `--resume-checkpoint state.h5`.
The paired `state.h5.resource.bin` stores resource time, field, front mask and
UID refractory state, checks the ABM checksum and dynamics fingerprint, and is
required for restart. Both files must be retained together. Changing thread
count and extending the end time preserve the trajectory.

## Ensemble validation

Run `python3 scripts/validate_abm_vs_pde.py --exe build-codex/atcg3d_shared_abm
--config ATCG3D_SharedRules/config/validation_v7.yaml --seeds 16 --output results`
(on one shell line), or `ctest --test-dir build-codex -L validation`.
The default is 256 square voxels and 48 hours. Reports contain total mass,
r/K ratio, active fraction, r50/r90/r99 and normalized radial density L2.
Scalar tolerances are 35% mass/ratio, 5 percentage points active fraction,
25% r50/r90 and 30% r99, plus an explicitly reported paired sampling allowance.
The area-weighted radial L2 tolerance is 0.35 with no sampling allowance.
Reports use deterministic seed ordering and retain each seed's raw metrics.

The original CI case uses near-memoryless division work clocks (minimum fraction zero,
stochastic fraction one), sparse small cells and no angiogenesis. It does not
establish dense invasive agreement. The original regular-cycle case (minimum
fraction 0.9) failed the mass tolerance: 468.69 ABM versus 650.14 PDE mean cells.
This records a model closure limitation, not a numerical tolerance adjustment.
Structured schema 10 now supplies a transported division-work distribution;
the regular-cycle ensemble passes with the same tolerances. See
[the renewal validation report](renewal_validation.md) for the law, numerical
checks and closure limits. Individual growth and activation clocks, footprint correlations, direction
correlations and mixed refractory ages need broader regime validation.
