# Activated invasion and nonlinear vascular verification

The follow-up corrects shared ABM initial active-jump scheduling in model v4,
adds transported Beta activation durations in structured schema 11 and active
Beta speed distributions in schema 12, and adds feasible-subset persistent
transport with direction reset after an unsuccessful choice. The native PDE
entry initializes v4 agents through the same shared environment as the ensemble
adapter. Published model strings and checkpoint arithmetic remain unchanged.
See [activation distributions](activation_distribution.md) for equations,
numerical grids, checkpoint formats and approximation boundaries.

The old mean-duration/mean-speed closure failed activated comparisons. Retaining
direction after a blocked cone also trapped PDE mass in the core. These defects
were corrected without enlarging the previously declared tolerances. All active
cases use 16 paired seeds and require at least 10% active r in both models.

| Metric | r20 ABM / PDE, 256-square, 8 h | r200 ABM / PDE, 128-square, 2 h |
| --- | ---: | ---: |
| Total cell mass | 930.688 / 970.377 | 982.375 / 989.702 |
| r/K ratio | 0.894787 / 0.915090 | 0.979858 / 0.981857 |
| Active r fraction | 0.315805 / 0.274721 | 0.311173 / 0.295504 |
| r50 | 13.125 / 14.000 | 13.000 / 14.000 |
| r90 | 22.938 / 22.500 | 33.250 / 26.688 |
| r99 | 63.500 / 57.688 | 71.562 / 69.250 |
| Normalized radial L2 error | 0.055547 | 0.081327 |

Both ensembles pass the original tolerances, including the active-fraction
coverage check. Reports retain every realization and explicitly state the
paired sampling allowance. The r200 CI case is an early invasion check; it
does not replace the separate 256-square, 48-hour regular-cycle growth case
or establish long-duration invasive calibration.

The nonlinear vascular comparison enables PDE branching at 0.1/hour and
anastomosis at 0.2/hour. Sixteen ABM realizations commit 58 roots and 2 actual
anastomoses. PDE endpoint branching and anastomosis rates are 1.522818 and
0.817897/hour. Volume-derived length and perfused volume are 26.6875 ABM
versus 44.8097 PDE, relative error 0.6791. Lesion perfusion is 0.037531 versus
0.120273, absolute error 0.082742. This passes the existing broad early-smoke
tolerances (0.75 relative, 0.15 absolute, plus reported sampling uncertainty).
The previously supplied comparison with these terms disabled is retained.
Continuous tip branching and discrete hypoxic roots are different closures;
the length proxy is not a centreline-length calibration.

New `atcg3d_activation_distribution_test` checks independent Beta identities
and moments, selective rate transport, cohort expiry and four-thread native
restart. Extended structured tests check doorway flux, the analytic Bernoulli
choice average, direction reset, target capacity and conservation. Shared ABM
tests check initial jump scheduling and event-time resource restart. CLI tests
check identical final native/ensemble PDE checksums for a resource-initialized
active case. The new active and nonlinear vascular CTests carry the validation
label; the workflow retains their reports.

Full Release verification passes 42/42 default CTests and 44/44 HDF5 CTests,
including all six 16-seed validation suites and published checksum fixtures.
The added native-entry CLI comparison also passes in both builds. Remaining
limits are distribution correlations, the mean-rate swap closure and broad
early vascular tolerances. Mixed-model exchange and convergence verification
are handled separately.
