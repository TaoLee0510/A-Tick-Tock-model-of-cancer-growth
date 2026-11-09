# Activation duration and speed distributions

The shared ABM model `nutrient_gradient_shared_resource_v4` schedules an
initially activated cell's first jump with its active rate. Model v3 retains
its original initial ordinary-rate scheduling for checkpoint compatibility.
Later migration, growth and per-UID refractory rules are shared by v3 and v4.

Structured schema 11 adds `beta_duration_distribution_v2`. Each stage carries
a sparse remaining-activation-time distribution in addition to directional
mass. New activation deposits the configured Beta law scaled by the full
inherent division cycle. Initial ABM import deposits actual remaining times.
Migration and r/K swaps transfer the local frozen distribution with their
accepted mass. Time advection expires only the completed cohort, transfers it
to ordinary refractory r, and preserves the mean clock of surviving mass.
This prevents a short-lived cohort from extending or prematurely expiring a
long-lived cohort through an averaged deadline.

The kernel uses the regularized incomplete Beta function, evaluated with the
continued fraction and symmetry identities in
[NIST DLMF 8.17](https://dlmf.nist.gov/8.17). Interval probabilities and first
moments are deposited on adjacent time nodes. Both total probability and the
configured mean are preserved. Tests independently check elementary Beta
CDFs, extreme-parameter symmetry, kernel moments and mixed-cohort expiry.
The configured time grid must cover the longest initial full inherent cycle.

Structured schema 12 adds `beta_rate_distribution_v2` with 4 to 64 rate bins.
The configured unclamped ordinary-r Beta law determines velocity nodes; the
active multiplier and stage mobility scale these nodes. Each migration substep
weights outgoing mass by the node's Poisson jump probability, so fast cells
preferentially reach the front. The rate supremum sets the substep bound.
Accepted flux debits the same weighted rate distribution, conserves mass,
and leaves rejected flux at its source. Ordinary refractory r continues to
use the existing diffusion closure; this is an active-population refinement.

The new `nutrient_gradient_feasible_direction_jump_exchange_v5` model averages
the literal persistent ABM choice over independent Bernoulli occupancies of
forward and turn destinations. A blocked forward site transfers its prior
to feasible turns. Nutrient guidance weights the choices within each feasible
subset. Target-volume budgets jointly limit incoming flux; rejected proposals
remain at their sources. An attempted jump with no feasible cone clears its
direction history; its next attempt can select any legal direction, as in the
ABM. The cone is restricted to at most 45 degrees to bound
subset enumeration in 3D. Earlier transport strings retain their arithmetic.

Division-work, duration and speed distributions are spatial marginals. Their
correlations with direction, phenotype stage and one another are approximated.
The r/K swap operator retains its existing mean-rate event closure. Active
direction densities and moment fields retain native float precision. Numerical
time/rate discretization and aggregate ensemble comparisons do not establish
pathwise equivalence to individual deterministic ABM migration clocks.

Native schema-11/12 checkpoints store duration/rate distributions and validate
their configuration identity. Unit tests cover conservative selective rate
transport and native PDE restart with four threads. The ODE and hybrid v1
wrappers reject duration distributions because their mean-clock state does
not carry these distributions.

The active validation configurations retain the regular division law
(minimum fraction 0.9, tail fraction 0.1), 500 initial cells per phenotype,
and a nonzero activation threshold. The r20 case uses 256-square voxels for
eight hours; the bounded r200 CI case uses 128-square voxels for two hours.
A 256-square, eight-hour r200 configuration is supplied for longer runs.
`atcg3d_active_r20_validation` and `atcg3d_active_r200_validation` run 16 paired
seeds with the same declared scalar/profile tolerances as the original ensemble
and require at least 10% active r in both models. Per-seed reports and JSON/Markdown
aggregate reports are retained under the build directory. These early invasive
cases complement the separate 256-square, 48-hour regular-cycle growth case.
