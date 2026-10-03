# Ver7 simulator validation and shared resource models

The phase A-F commits and their evidence are included in PR #67. Following
the instruction to avoid hosted execution, active hosted runs were cancelled;
subsequent builds, tests and benchmarks ran locally on macOS. The delivery
commit includes `[skip ci]` to suppress automatic push and pull-request
workflow runs. Hosted Linux/macOS results for phases D-F are unavailable.

## Resulting behavior

The simulator exposes individual ABM, well-mixed ODE, continuum/structured
PDE and adaptive hybrid models through `atcg_sim`. Relative configuration
graphs, versioned rule contracts, migration, checkpoints and field output
make the models reproducible and independently testable.

The follow-up correctness work fixes hybrid large-cell activation density,
volume-based placement and distributed resource footprints. New vascular
schemas diagnose removed r/K mass by stage. New hybrid front policies retain
every r individual and classify K cores using volume, front distance,
active-r protection, nutrient gradients and hysteresis. Core transport
reflects at the interface; placement respects native swaps and self-overlap.
Regional conversion conserves mass and uses cached regional membership.

New shared VEGF angiogenesis gives individual tips and the conditional PDE
tip field the same source, hop, branch and connection rate definitions.
Diagnostics use centerline edge length, perfused volume, lesion perfusion,
branches and anastomoses. Deterministic parallel vascular work uses active
bounds and reusable buffers. The choice and stochastic closure are documented
in `shared_angiogenesis_contract.md`.

New prepared 3D sector means use cone row segments, summed-volume corner
kernels and tiled spectral convolution. Queries have constant cost after
resource preparation. A new sparse population storage policy supports 3D
and dynamic vasculature. Resource and vascular arrays remain dense.

Six immutable structured PDE policies select activation, transport, exchange,
reaction, nutrient and vascular contracts once. Published arithmetic and
traversal order are preserved. Resource-limited restart checks discovered an
older omission of positive field tails, active totals and resource work/front
buffers. Structured schema 17 persists the complete state; hybrid v5 also
initializes resources and individual clocks through the standalone ABM path.
Older schema formats and behavior remain available.

## Evidence and interpretation

Statistical validation uses at least 16 paired realizations, an independent
16-seed ABM baseline, time-series comparisons, paired t tests and TOST with
prespecified margins. Existing smoke allowances remain separate. All scalar
contrasts and observation times must pass; the radial-profile interval is a
nonlinear jackknife approximation. Passing a broad equivalence margin can
coexist with a statistically significant model bias.

The regular 48-hour and early r20 8-hour ABM/PDE cases pass the declared
equivalence screen. Neighbor birth placement improves the regular case's
maximum radial-profile error from 0.27458671035633075 to 0.22627052907292164.
The r20 PDE endpoint r90 is 21.125 versus ABM 22.9375, and r99 is 60.9375
versus 63.5. These residual front biases remain visible in the raw reports.
Mechanism ablations measure migration, activation, exchange and boundary
effects; the short no-division case contains no births and cannot establish
division accuracy.

The 256-square, 240-hour angiogenesis experiment passes the prespecified
controlled screen with frozen hypoxic consumers. Endpoint ABM/PDE centerline
lengths are 530 and 684.54572450997489; perfused volumes are
256.30050099315639 and 368.65860492795275. Relative errors are
0.29159570662259404 and 0.43838425402764436. Lesion-perfusion error is
0.13599071061585738, near its 0.15 absolute margin. The rate contract is
shared, while spatial correlation and independent-arrival closures explain
material differences. This experiment does not establish growing-tumor
vascular equivalence. ABM ending-population and PDE starting-population
source sampling also needs a coupled splitting study.

The invasive hybrid case passes the ABM-primary, PDE-secondary and step/
exchange refinement screens. Its minimum ABM mass fraction is
0.6714706209799789 and minimum active fraction is 0.24675324675324675, above
the declared 0.20 and 0.10 requirements. Maximum hybrid/ABM errors are
0.0018419799339300113 for mass, 0.003222270324416123 for r/K,
0.00073267345559853425 for active fraction, 0.1008174386920981 for r90 and
0.040654172338689773 for profile L2. The corrected resource-initialized
hybrid repeats this high-resource case under the same statistical margins.
Low-nutrient initialization and dynamic-vascular restarts have separate
assert-based coverage; their growing-tumor equivalence is unmeasured.

The older 128-square r200 case has boundary contact and fails the statistical
screen despite passing smoke checks. Its r90 TOST lower p is
0.05393702365741495. The earlier nonlinear vascular smoke case also fails
the new statistical screen. Both failed results and successful experiments
are retained; no declared tolerance has been widened.

## Verification and performance

| Phase | Release CTests | HDF5 CTests | Validation suites per build |
| --- | ---: | ---: | ---: |
| A | 47/47 | 49/49 | 8 |
| B | 52/52 | 55/55 | 10 |
| C | 56/56 | 59/59 | 11 |
| D, full run plus documented retry | 60/60 | 63/63 | 12 |
| E | 62/62 | 65/65 | 12 |
| F | 65/65 | 68/68 | 13 |

Phase F's full-suite times are 4895.79 seconds for Release and 4914.24
seconds for HDF5. An additional expanded schema-17 restart check passes
in both builds, including a migrated hybrid-v4 wrapper; it takes 6.83 and
6.99 seconds respectively. These builds and tests ran locally on macOS.

Each phase has a separate commit and full local CTest evidence, including
the validation label. Phase D's complete run had a compatibility-scanner
failure caused by later-phase configurations in the shared source tree.
Recompiling only that unchanged test against immutable phase-D sources
passed the retry; its report records the original failure and combined
test-set result. Simulator libraries and completed ensembles were unchanged.

The final maintenance comparisons check 416 published native reports,
96 schema-14 reports, 48 schema-15 vascular reports and 176 hybrid-v4 reports
exactly against preceding fixtures in both builds. C++ assert tests compile
with `-UNDEBUG`; checked compiler and linker logs contain no new warnings.
The largest measured direct 3D sector-mean error is
1.2086454059812013e-10, within its separately prespecified arithmetic budget.

The phase E/F reports record the full test counts and isolated 1/4/8-thread
256-cube measurements, plus production-prefix wall time, peak resident memory,
checkpoint size and exact cross-process continuation. Production inputs use
2000-square r200 transient nutrient and VEGF angiogenesis with a configured
2160-hour horizon. Measured prefixes preserve that horizon; they do not
establish feasibility over the full duration.

| 256-cube model | Threads | Median seconds | Speedup | Maximum resident bytes |
| --- | ---: | ---: | ---: | ---: |
| ABM | 1 | 82.891666458 | 1 | 4435918848 |
| ABM | 4 | 24.356374708 | 3.4032842511153243 | 4738252800 |
| ABM | 8 | 13.043498458 | 6.355017921373683 | 5141200896 |
| PDE | 1 | 1065.031624458 | 1 | 13064814592 |
| PDE | 4 | 677.232565375 | 1.572623170990406 | 13369917440 |
| PDE | 8 | 602.862700916 | 1.766623847917233 | 13775110144 |

All 18 3D measurements have exact biological diagnostics and checksums
across threads and repetitions within each model. ABM/PDE endpoint masses
are 475 and 511.99999007859481, active masses are 146 and
79.319266422425045, and centerline growth values are 223 and
278.57877937237703. These single-trajectory model differences need a separate
3D ensemble and mechanism study. Peak PDE memory reaches 13775110144 bytes;
the sparse policy has not established an 8-GB bound for this 3D workload.

| Production model | Scope | Threads | Seconds | Resident bytes | Checkpoint bytes | Mass | Active mass | r99 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| abm | reference | 1 | 0.65745583399999996 | 796393472 | 0 | 905 | 46 | 18 |
| abm | checkpoint_prefix | 1 | 0.877671333 | 799506432 | 228169845 | 990 | 42 | 18 |
| abm | resumed | 8 | 0.45349208299999999 | 773865472 | 0 | 905 | 46 | 18 |
| pde | reference | 1 | 1598.1177935840001 | 8387461120 | 0 | 971.02563189167199 | 153.90012883284106 | 73 |
| pde | checkpoint_prefix | 1 | 405.01444720799998 | 4942020608 | 2001976223 | 995.34471763408965 | 155.43032215426302 | 47 |
| pde | resumed | 8 | 1199.6755417500001 | 8939520000 | 0 | 971.02563189167199 | 153.90012883284106 | 73 |
| hybrid | reference | 1 | 2.6449263749999998 | 2724380672 | 0 | 984.47457234368994 | 44 | 18 |
| hybrid | checkpoint_prefix | 1 | 4.1789416250000002 | 2714533888 | 1236775335 | 998.73465277777746 | 43 | 18 |
| hybrid | resumed | 8 | 2.8618167909999999 | 2689449984 | 0 | 984.47457234368994 | 44 | 18 |

All three production-prefix restarts are exact. Each observation remains
inside the prespecified r99/half-width guard. The one-seed endpoint masses
are ABM 905, PDE 971.02563189167199 and hybrid 984.47457234368994.
Endpoint r99 values are 18, 73 and 18 respectively; active masses are
46, 153.90012883284106 and 44. These material r200 model differences have
not passed an ensemble equivalence study. The inputs remain benchmark
workloads, rather than validated scientific recommendations. The eight-thread
PDE restore peaks at 8939520000 resident bytes; the workload has not
established a decimal 8-GB bound. Checkpoint and resource sidecar sizes are
reported together. All raw state/field checksums are in the phase F report.

## Scope still requiring evidence

- Full 2160-hour production runs and the opt-in `validation_long` ensembles.
- Long-time 3D biological equivalence and growing-tumor vascular equivalence.
- Stronger vascular accuracy than the controlled screen's broad margins.
- A memory bound for a fully occupied 3D domain; dense fields and caches
  contribute substantial memory even with sparse population storage.
- Linux/macOS hosted verification of phases D-F, intentionally not run.

Native binary checkpoint continuation requires a compatible ABI. Current
recommended configuration graphs and migration/deprecation policy live in
`recommended_configurations.md`; every other input is registered as a
dependency, benchmark, legacy or reproduction configuration.

Per-phase changes, added tests, original numerical tables, failed pilots and
source digests are retained in `phase_a_report.md` through `phase_f_report.md`
and their corresponding numeric JSON reports.
