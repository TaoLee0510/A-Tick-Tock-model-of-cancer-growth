# Shared ABM and structured PDE rules

`atcg3d_shared_abm` accepts the same structured YAML as the PDE.
Use `--model abm` or `--model pde`, `--seed N`, `--threads N` and
`--report summary.json`. The ABM adapter selects the new
`nutrient_gradient_shared_resource_v3` model and a finite unit-spaced resource
grid. Native individual positions remain expandable outside that grid; the PDE
reflects at its grid faces. Boundary-exposed comparisons therefore mix different
domain closures and cannot establish free-space equivalence.
An explicitly selected v4 model keeps the same contract and corrects the first
jump time of initially activated cells to use their active rate.
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
--config ATCG3D_SharedRules/config/regular_cycle_spatial_birth_v14.yaml
--seeds 16 --output results`
(on one shell line), or `ctest --test-dir build-codex -L validation`.
The default is 256 square voxels and 48 hours. Reports contain total mass,
r/K ratio, active fraction, r50/r90/r99 and normalized radial density L2.
Equivalence margins are 35% mass/ratio, 5 percentage points active fraction,
25% r50/r90 and 30% r99. The area-weighted radial L2 margin is 0.35.
Every metric at every four-hour sample must pass paired TOST at alpha 0.05.
Two disjoint 16-seed ABM groups quantify within-model variation; the paired
PDE group uses the first group's seeds. Reports retain each realization,
baseline ratios, descriptive paired t tests, both one-sided equivalence tests
and the maximum error over time. See [the statistical protocol](statistical_validation.md)
for the prespecified hypotheses, boundary guard and profile approximation.

`--smoke-only` selects the historical endpoint criterion, which adds an
explicit paired sampling allowance to scalar tolerances. Its success is a
smoke result and does not establish statistical equivalence. Published smoke
CTest cases retain this criterion without changing their tolerances.

The original CI case uses near-memoryless division work clocks (minimum fraction zero,
stochastic fraction one), sparse small cells and no angiogenesis. It does not
establish dense invasive agreement. The original regular-cycle case (minimum
fraction 0.9) failed the mass tolerance: 468.69 ABM versus 650.14 PDE mean cells.
This records a model closure limitation, not a numerical tolerance adjustment.
Structured schema 10 now supplies a transported division-work distribution;
the regular-cycle ensemble passes with the same tolerances. See
[the renewal validation report](renewal_validation.md) for the law, numerical
checks and closure limits. [Activation distributions](activation_distribution.md)
describe the opt-in schema 11/12 duration and speed closures and the new active
ensemble coverage. Footprint, direction and mixed refractory-age correlations
remain approximate.

Schema 14 optionally uses the true truncated-normal growth expectation,
uniform feasible eight/26-direction normal jumps, and neighbor placement of
small daughters. Each operator can be disabled independently for paired
mechanism attribution. Published operators remain the default. See
[the phase B report](phase_b_report.md) for all raw statistical and intervention
tables, improvements in front position, remaining r/K bias, and the published
r200 and vascular cases that do not establish equivalence.
