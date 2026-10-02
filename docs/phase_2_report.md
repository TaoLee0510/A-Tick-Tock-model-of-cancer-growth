# Phase 2 verification

The new shared-resource ABM environment loads structured YAML and applies
transient nutrient, common r/K growth limits, resource-scaled positive growth,
70-square sector chemotaxis, UID-specific cooldown and hysteresis. Coupled
resource checkpoints bind to the HDF5 ABM checksum. The new lattice-mean
mapping corrects the thin-layer diffusion second moment without changing
published mapping strings. Published PDE checksum fixtures remain identical.

New tests are `atcg3d_shared_rules_test` and `atcg3d_abm_pde_validation` (label
`validation`). Unit coverage includes identical transient resource updates,
common growth rates, brute-force versus cached sector averages, separate UID
cooldowns, resource/ABM event-time restart with four threads, 3D cached sectors
and the eight-direction lattice second moment. The ensemble uses 16 paired
seeds, 256-square voxels and 48 hours. It retains per-seed data and emits JSON
and Markdown reports.

Release checks: default 30/30 CTests and HDF5 31/31 CTests pass. A command-line
24-hour HDF5 checkpoint resumed with four threads reaches the same complete
48-hour report as the one-thread uninterrupted seed-5 run. ABM checksum is
10529943340798383128; resource checksum is 11395177022692243872.

The ensemble passes the predeclared tolerances: ABM/PDE mean mass 743.31/674.39,
mass error 0.0927, ratio error 0.1887, r50/r90/r99 errors 0.0593/0.1837/0.1983,
and radial L2 0.3270. Both active fractions are zero in this sparse smoke case.
The mass tolerance is 0.35 plus paired sampling uncertainty; radial tolerance
is a fixed 0.35. These tolerances were not enlarged to make the run pass.

Remaining limits: the regular-cycle baseline failed the mass comparison and
is documented in the contract. The passing baseline uses near-memoryless work
clocks matching a PDE mean-rate closure. Dense activation, correlations among
cell footprints and directions, mixed refractory ages, and angiogenesis require
additional validation regimes. The sparse test's zero active fraction does not
validate activated invasion quantitatively. ODE, dynamic vessels and hybrid
coupling remain subsequent phases.

Follow-up: [transported division work in schema 10](renewal_validation.md)
addresses the regular-cycle mass bias with an independently checked renewal
law. The new 16-seed case passes the original tolerances. This does not resolve
the separate activated-invasion and age/direction-correlation limits above.

The next [activated verification follow-up](activated_vascular_validation.md)
adds duration/rate distributions, corrects blocked-direction reset and initial
active scheduling, and passes nonzero r20/r200 ensembles at the original
tolerances. Native dispatch and the ensemble adapter also agree bitwise.
