# Phase B: statistical equivalence and mechanism attribution

The two new default statistical validation suites use 16 paired ABM/PDE
realizations and 16 additional disjoint ABM baseline seeds. They require
every metric at every prescribed time to pass paired TOST at alpha 0.05.
The historical smoke checks retain their original sampling allowance;
their success is not a statistical equivalence claim. All previously
declared margins are unchanged. See [the protocol](statistical_validation.md)
for the prespecified margins, hypotheses, profile jackknife approximation
and boundary rule.

Schema 14 adds independent operator switches, the true truncated-normal
growth-rate expectation, an eight/26-direction uniform-feasible normal
jump closure and optional neighbor placement of small daughters. Their
published alternatives remain unchanged. The ordinary daughter stays
at its maternal voxel only under the retained local-growth closure.
Newborns in the neighbor model receive fresh division work; competing
prospective births are limited by recipient volume. These are mean-field
closures, not exact representations of correlated individual occupancy.

## Numerical findings

Regular-cycle runs cover 0 through 48 hours every four hours. The latest
neighbor-birth closure passes all prescribed contrasts. Its maximum r90
error is 0.12285714285714286, r99 error 0.10739856801909307, and radial
profile error 0.22627052907292164. The published closure gives
0.15083798882681565, 0.1288782816229117 and 0.27458671035633075.
The normal transport correction reduces the front lag; a systematic
tail lag remains. The maximum r/K error is 0.15075424357325423, slightly
worse than the published 0.13288175690618612. Correcting a distribution
mean does not remove inherited growth-rate and work/age correlations.

The no-migration experiment exposed the small-daughter location error.
Its former local-birth mass error is 0.44602414040553423 and fails
equivalence. Neighbor birth reduces that error to 0.05698273276152335
(ABM 442.9375; PDE 468.1772891925572) and passes every sampled contrast.
The old failure is retained below rather than replaced by its correction.

The r20 v14 eight-hour activation case passes every contrast at 0, 4 and
8 hours, with endpoint ABM/PDE active fractions 0.3158049654902044 and
0.30992378821835553. Its maximum active-fraction error is
0.005881177271848892. It uses transported division work rather than
the old mean-rate birth approximation. Activation/exchange interventions
use this active case; the short duration leaves division unexercised.

At eight hours disabling crowding exchange changes the paired r99
model difference by +8.75 voxels, versus +2.5625 without migration
and +2.25 without activation. These intervention effects interact
and are not additive components of one error. Doubling the domain
changes all three radius quantiles by exactly zero; its mass interaction
is -9.237055564881302e-14. This case is interior under the r99 guard,
while the original PDE still has a tiny nonzero outer-face tail.

The published 128-square/two-hour r200 case remains a smoke reproduction.
Its r90 TOST lower p value is 0.05393702365741495, above 0.05, so its
statistical result is inconclusive. Its r99 reaches/exceeds the boundary
guard (maximum ratios 1.203125 ABM, 1.21875 independent ABM, 1.09375 PDE).
It cannot support a free-space equivalence claim. The published shared
ABM has expandable sparse positions while the PDE reflects at its grid
faces; the resource grid does not constrain native individual positions.
These boundary-exposed results also mix different domain closures.

The published nonlinear vascular smoke also fails the new statistical
criterion: its endpoint perfused-volume proxy has relative error
0.6790534066033419 and upper TOST p 0.37057439858497371. Its initial
vascular reference is zero, making relative margins undefined. These
failures remain visible. Its reported vascular length is still the old
volume/cross-section proxy, not the centerline quantity required by
phase C.

## Verification and remaining work

Release CTest passes 52/52 tests in 1784.25 seconds; the HDF5 build
passes 55/55 in 1760.42 seconds. Both runs include all ten tests labeled
`validation`. Four additional checks after the exact cache optimization
pass, including the compatibility and published PDE regressions. Neither
build emits a new compiler or linker warning. In each build all 416
published per-seed native reports equal the phase A baseline exactly;
all 96 saved schema 14 traces also remain exactly equal. The two new
statistical datasets are identical between Release and HDF5 builds.

Added tests cover analytic Student probabilities and NIST quantiles,
strict deterministic TOST boundaries, noisy non-equivalence despite a
smoke pass, baseline cancellation, nonlinear profile uncertainty,
trajectory failure despite endpoint success and misaligned interventions.
Native trace tests verify unchanged endpoints, one/four-thread identity
and HDF5 trace/restart identity. Operator tests cover enumerated feasible
choice probabilities, 2D/3D transport and daughter placement, exact clock
expiry without migration, ABM event disabling and checkpoint versioning.

Four opt-in `validation_long` workloads are registered behind
`ATCG_ENABLE_LONG_VALIDATION=ON` and successfully dry-run. They have not
been executed. The r20 domains are 512-square/360 hours and
2000-square/720 hours; both r200 domains are enlarged to 4096-square
with the same durations and a strict r99/half-width < 0.8 guard.
They retain the one-million active-voxel budget. Registration is not
production-scale evidence. These opt-in configurations retain local
birth/mean-rate division and therefore do not inherit the demonstrated
neighbor-birth/renewal accuracy without further validation.

Vascular mechanism equivalence, centerline statistics, true hybrid-front
coverage, 3D/production benchmarks and operator-strategy refactoring
remain phases C through F. Existing hybrid and vascular passes below
are only historical smoke evidence. The eight-hour active ablation
does not characterize division-driven phenotype changes at 240 hours.

All table values use binary64 round-trip precision. Complete profile
arrays, every scalar realization, baseline diagnostics, TOST contrasts
and intervention arrays are retained in
[phase_b_numeric_results.json](phase_b_numeric_results.json). Per-seed
native reports are retained in build validation directories and CI
artifacts. The unchanged smoke allowance is displayed explicitly.

CTest statistics use Python 3.11; manual diagnostic studies use Python
3.14. Their aggregate summation can differ in the last binary64 digits.
Tables retain their generating report values, and raw native traces
are identical. No acceptance result changes between those interpreters.

## Historical smoke suites

### validation

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 743.3125 | 674.38546078935178 | 9.8053030358293007 | 0.092729557501923104 | 0.34999999999999998 | 0.028110788893436121 | PASS |
| r_K_ratio | 1.2542428737808808 | 1.0176043954212508 | 0.03112628731277876 | 0.1886703790042592 | 0.34999999999999998 | 0.052884588503645399 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.050000000000000003 | 0 | PASS |
| r50 | 15.8125 | 14.875 | 0.11063265039459795 | 0.059288537549407112 | 0.25 | 0.014909608094285421 | PASS |
| r90 | 24.5 | 20 | 0.18257418583505536 | 0.18367346938775511 | 0.25 | 0.015880228163857261 | PASS |
| r99 | 29.3125 | 23.5 | 0.27716947282604315 | 0.19829424307036247 | 0.29999999999999999 | 0.020150043380547478 | PASS |
| radial_profile_L2 | null | null | null | 0.32698891172682604 | 0.34999999999999998 | 0 | PASS |

Coverage: `{"active_fraction": {"abm": 0, "pde": 0}}`.

Boundary: `{"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": null, "maximum_r99_to_half_width": 0.25, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": null, "maximum_r99_to_half_width": 0.1875, "outer_face_mass_present": false}}`.

### validation-regular

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 468.375 | 489.72804166911277 | 1.8941015470146889 | 0.045589627262583976 | 0.34999999999999998 | 0.0086177323654941049 | PASS |
| r_K_ratio | 1.4233832881828432 | 1.3285012486394163 | 0.011453123436082308 | 0.066659514925567015 | 0.34999999999999998 | 0.01714689658429951 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.050000000000000003 | 0 | PASS |
| r50 | 14.4375 | 13.1875 | 0.11180339887498948 | 0.086580086580086577 | 0.25 | 0.016502375272907537 | PASS |
| r90 | 22.375 | 19 | 0.15478479684172258 | 0.15083798882681565 | 0.25 | 0.014741738639987075 | PASS |
| r99 | 26.1875 | 22.8125 | 0.125 | 0.12887828162291171 | 0.29999999999999999 | 0.010171837708830548 | PASS |
| radial_profile_L2 | null | null | null | 0.27458671035633064 | 0.34999999999999998 | 0 | PASS |

Coverage: `{"active_fraction": {"abm": 0, "pde": 0}}`.

Boundary: `{"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": null, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": null, "maximum_r99_to_half_width": 0.1796875, "outer_face_mass_present": false}}`.

### validation-vascular

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| vascular_length | 26.6875 | 33.682986743131785 | 3.1983212011413391 | 0.26212596695575774 | 0.75 | 0.25538632242181519 | PASS |
| perfused_volume | 26.6875 | 33.682986743131785 | 3.1983212011413391 | 0.26212596695575774 | 0.75 | 0.25538632242181519 | PASS |
| lesion_perfused_fraction | 0.037530736013472225 | 0.090234588360324455 | 0.0039945825957257361 | 0.05270385234685223 | 0.14999999999999999 | 0.0085124555114915422 | PASS |

Coverage: `{"active_fraction": {"abm": 0, "pde": 0}}`.

Boundary: `{"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": null, "maximum_r99_to_half_width": 0.4583333333333333, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": null, "maximum_r99_to_half_width": 0.375, "outer_face_mass_present": false}}`.

### validation-active-r20

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 930.6875 | 970.37748642098722 | 2.1778784515893821 | 0.042645878902410554 | 0.34999999999999998 | 0.0049866995960910323 | PASS |
| r_K_ratio | 0.89478658768318697 | 0.91508983096094598 | 0.0037481771122353726 | 0.02269059858209194 | 0.34999999999999998 | 0.0089265591774847081 | PASS |
| active_fraction | 0.31580496549020443 | 0.27472069115986003 | 0.0024065440219043521 | 0.041084274330344395 | 0.050000000000000003 | 0.0051283453106781736 | PASS |
| r50 | 13.125 | 14 | 0.085391256382996647 | 0.066666666666666666 | 0.25 | 0.013864287036355491 | PASS |
| r90 | 22.9375 | 22.5 | 0.61892884620662714 | 0.019073569482288829 | 0.25 | 0.057501356785452748 | PASS |
| r99 | 63.5 | 57.6875 | 0.59314662324476453 | 0.091535433070866146 | 0.29999999999999999 | 0.019905440222592018 | PASS |
| radial_profile_L2 | null | null | null | 0.055546805567441175 | 0.34999999999999998 | 0 | PASS |

Coverage: `{"active_fraction": {"abm": 0.3158049654902044, "pde": 0.27472069115986003}}`.

Boundary: `{"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": null, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": null, "maximum_r99_to_half_width": 0.4765625, "outer_face_mass_present": false}}`.

### validation-active-r200

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 982.375 | 989.7017335674069 | 1.4657306528003706 | 0.0074581840614906776 | 0.34999999999999998 | 0.0031795109007431883 | PASS |
| r_K_ratio | 0.97985804784240171 | 0.98185735548383302 | 0.0024945664867440129 | 0.0020404053891619075 | 0.34999999999999998 | 0.0054251952055268437 | PASS |
| active_fraction | 0.31117281779932982 | 0.29550364174646931 | 0.00085367000780809125 | 0.015669176052860501 | 0.050000000000000003 | 0.0018191707866390423 | PASS |
| r50 | 13 | 14 | 0 | 0.076923076923076927 | 0.25 | 0 | PASS |
| r90 | 33.25 | 26.6875 | 1.1760412053438718 | 0.19736842105263158 | 0.25 | 0.075372746122941078 | PASS |
| r99 | 71.5625 | 69.25 | 0.58251430597597054 | 0.032314410480349345 | 0.29999999999999999 | 0.017346207665115011 | PASS |
| radial_profile_L2 | null | null | null | 0.081326884302303726 | 0.34999999999999998 | 0 | PASS |

Coverage: `{"active_fraction": {"abm": 0.3111728177993298, "pde": 0.2955036417464693}}`.

Boundary: `{"abm": {"annotated": true, "boundary_influenced": true, "maximum_boundary_mass": null, "maximum_r99_to_half_width": 1.203125, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": true, "maximum_boundary_mass": null, "maximum_r99_to_half_width": 1.09375, "outer_face_mass_present": false}}`.

### validation-vascular-nonlinear

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| vascular_length | 26.6875 | 44.809737788726686 | 3.2142268696439849 | 0.6790534066033419 | 0.75 | 0.25665639191424189 | PASS |
| perfused_volume | 26.6875 | 44.809737788726686 | 3.2142268696439849 | 0.6790534066033419 | 0.75 | 0.25665639191424189 | PASS |
| lesion_perfused_fraction | 0.037530736013472225 | 0.12027296915844937 | 0.004029409737840264 | 0.08274223314497714 | 0.14999999999999999 | 0.0085866721513376022 | PASS |

Coverage: `{"active_fraction": {"abm": 0, "pde": 0}, "vascular_activity": {"abm_anastomoses": 2, "abm_roots": 58, "pde_anastomosis_rate": 0.8178972356693133, "pde_branching_rate": 1.5228178080433576}}`.

Boundary: `{"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": null, "maximum_r99_to_half_width": 0.4583333333333333, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": null, "maximum_r99_to_half_width": 0.375, "outer_face_mass_present": false}}`.

### validation-hybrid: adaptive/abm

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 468.375 | 482.23076250929751 | 1.9805545246701453 | 0.029582626120731263 | 0.34999999999999998 | 0.009011073802128804 | PASS |
| r_K_ratio | 1.4233832881828432 | 1.2942909243521372 | 0.013278079148936769 | 0.090694028026359133 | 0.34999999999999998 | 0.019879105579852424 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.050000000000000003 | 0 | PASS |
| r50 | 14.4375 | 13.0625 | 0.125 | 0.095238095238095233 | 0.25 | 0.018450216450216449 | PASS |
| r90 | 22.375 | 19 | 0.15478479684172258 | 0.15083798882681565 | 0.25 | 0.014741738639987075 | PASS |
| r99 | 26.1875 | 25.9375 | 0.89209491273817576 | 0.0095465393794749408 | 0.29999999999999999 | 0.072593957385968591 | PASS |
| radial_profile_L2 | null | null | null | 0.29566126697015688 | 0.34999999999999998 | 0 | PASS |

### validation-hybrid: adaptive/pde

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 489.72804166911277 | 482.23076250929751 | 0.66320581017523117 | 0.015309066506101438 | 0.34999999999999998 | 0.0028858702406882289 | PASS |
| r_K_ratio | 1.3285012486394163 | 1.2942909243521372 | 0.0036041601668488845 | 0.025751066716960622 | 0.34999999999999998 | 0.0057813007879525258 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.050000000000000003 | 0 | PASS |
| r50 | 13.1875 | 13.0625 | 0.085391256382996647 | 0.0094786729857819912 | 0.25 | 0.013798579514856177 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.25 | 0 | PASS |
| r99 | 22.8125 | 25.9375 | 0.88447253584645957 | 0.13698630136986301 | 0.29999999999999999 | 0.082621850910194208 | PASS |
| radial_profile_L2 | null | null | null | 0.024348409394435524 | 0.34999999999999998 | 0 | PASS |

### validation-hybrid: step/coarse_vs_fine

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 486.9802930790737 | 482.23076250929751 | 0.49424241999384083 | 0.0097530241721814192 | 0.10000000000000001 | 0.0021627786831937691 | PASS |
| r_K_ratio | 1.340664755305145 | 1.2942909243521372 | 0.0025524528120272143 | 0.034590176827951825 | 0.10000000000000001 | 0.0040571492022194429 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r99 | 26.25 | 25.9375 | 0.70544755297612305 | 0.011904761904761904 | 0.10000000000000001 | 0.057268904205414015 | PASS |
| radial_profile_L2 | null | null | null | 0.0043438312838042529 | 0.14999999999999999 | 0 | PASS |

### validation-hybrid: step/middle_vs_fine

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 486.9802930790737 | 485.863310751124 | 0.26169794670775864 | 0.0022936910257440927 | 0.10000000000000001 | 0.001145176370296529 | PASS |
| r_K_ratio | 1.340664755305145 | 1.3279632083777151 | 0.0012631369495516359 | 0.0094740664115831839 | 0.10000000000000001 | 0.0020077687795125748 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r99 | 26.25 | 26.375 | 0.42695628191498325 | 0.0047619047619047623 | 0.10000000000000001 | 0.034660717590888734 | PASS |
| radial_profile_L2 | null | null | null | 0.0011602356935058633 | 0.14999999999999999 | 0 | PASS |

### validation-hybrid: exchange/coarse_vs_fine

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 485.6608560281644 | 485.85028715879059 | 0.61231255817774666 | 0.00039004817513068067 | 0.10000000000000001 | 0.002686726849159752 | PASS |
| r_K_ratio | 1.3272983106838141 | 1.3282066995729778 | 0.0028352873678673916 | 0.00068438939600222513 | 0.10000000000000001 | 0.0045521020649929248 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 18.9375 | 0.0625 | 0.0032894736842105261 | 0.10000000000000001 | 0.007009868421052631 | PASS |
| r99 | 24.6875 | 25.6875 | 0.92195444572928875 | 0.040506329113924051 | 0.10000000000000001 | 0.079582174130597025 | PASS |
| radial_profile_L2 | null | null | null | 0.0015720915166280418 | 0.14999999999999999 | 0 | PASS |

### validation-hybrid: exchange/middle_vs_fine

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 485.6608560281644 | 485.863310751124 | 0.69645645116231436 | 0.00041686440331081625 | 0.10000000000000001 | 0.0030559364194276822 | PASS |
| r_K_ratio | 1.3272983106838141 | 1.3279632083777151 | 0.003351661756234164 | 0.00050094066160488032 | 0.10000000000000001 | 0.0053811499231512594 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r99 | 24.6875 | 26.375 | 1.2803604635674544 | 0.068354430379746839 | 0.10000000000000001 | 0.11051941864758461 | PASS |
| radial_profile_L2 | null | null | null | 0.0019474915248959784 | 0.14999999999999999 | 0 | PASS |

### validation-hybrid: step refinement order

| Metric | Coarse error | Middle error | Sampling allowance | Result |
| --- | ---: | ---: | ---: | --- |
| total_mass | 0.0097530241721814192 | 0.0022936910257440927 | 0.0022898115650144776 | PASS |
| r_K_ratio | 0.034590176827951825 | 0.0094740664115831839 | 0.00414680985892969 | PASS |
| active_fraction | 0 | 0 | 0 | PASS |
| r50 | 0 | 0 | 0 | PASS |
| r90 | 0 | 0 | 0 | PASS |
| r99 | 0.011904761904761904 | 0.0047619047619047623 | 0.057388651130111393 | PASS |
| radial_profile_L2 | 0.0043438312838042529 | 0.0011602356935058633 | 0.0013372693805159147 | PASS |

### validation-hybrid: exchange refinement order

| Metric | Coarse error | Middle error | Sampling allowance | Result |
| --- | ---: | ---: | ---: | --- |
| total_mass | 0.00039004817513068067 | 0.00041686440331081625 | 0.0028729545688062816 | PASS |
| r_K_ratio | 0.00068438939600222513 | 0.00050094066160488032 | 0.0050094686622969902 | PASS |
| active_fraction | 0 | 0 | 0 | PASS |
| r50 | 0 | 0 | 0 | PASS |
| r90 | 0.0032894736842105261 | 0 | 0.0070098684210526301 | PASS |
| r99 | 0.040506329113924051 | 0.068354430379746839 | 0.11079996974126083 | PASS |
| radial_profile_L2 | 0.0015720915166280418 | 0.0019474915248959784 | 0.0023188905601691846 | PASS |

Mixed coverage: `{"abm_mass": 6.0625, "pde_mass": 476.1682625092975, "to_abm": 9, "to_pde": 256}`.

### validation-hybrid-volume: adaptive/abm

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 468.375 | 482.33519851790538 | 1.9742094932395273 | 0.029805601319253552 | 0.34999999999999998 | 0.0089822053484781041 | PASS |
| r_K_ratio | 1.4233832881828432 | 1.2947931434296158 | 0.013217386631744866 | 0.090341193282795593 | 0.34999999999999998 | 0.01978824055761301 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.050000000000000003 | 0 | PASS |
| r50 | 14.4375 | 13.0625 | 0.125 | 0.095238095238095233 | 0.25 | 0.018450216450216449 | PASS |
| r90 | 22.375 | 19 | 0.15478479684172258 | 0.15083798882681565 | 0.25 | 0.014741738639987075 | PASS |
| r99 | 26.1875 | 26.0625 | 0.89849411053532602 | 0.0047732696897374704 | 0.29999999999999999 | 0.073114690197643134 | PASS |
| radial_profile_L2 | null | null | null | 0.29532787173799768 | 0.34999999999999998 | 0 | PASS |

### validation-hybrid-volume: adaptive/pde

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 489.72804166911277 | 482.33519851790538 | 0.6657033356306552 | 0.015095813435576963 | 0.34999999999999998 | 0.0028967379596927792 | PASS |
| r_K_ratio | 1.3285012486394163 | 1.2947931434296158 | 0.0036029438539854127 | 0.025373032388432176 | 0.34999999999999998 | 0.0057793497452156726 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.050000000000000003 | 0 | PASS |
| r50 | 13.1875 | 13.0625 | 0.085391256382996647 | 0.0094786729857819912 | 0.25 | 0.013798579514856177 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.25 | 0 | PASS |
| r99 | 22.8125 | 26.0625 | 0.88741196746494244 | 0.14246575342465753 | 0.29999999999999999 | 0.082896434089547055 | PASS |
| radial_profile_L2 | null | null | null | 0.024220026542160759 | 0.34999999999999998 | 0 | PASS |

### validation-hybrid-volume: step/coarse_vs_fine

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 487.07775221479494 | 482.33519851790538 | 0.4009039965359234 | 0.0097367487538173453 | 0.10000000000000001 | 0.0017539836560658715 | PASS |
| r_K_ratio | 1.3411611329909079 | 1.2947931434296158 | 0.0020658388295670457 | 0.034573019170252385 | 0.10000000000000001 | 0.0032824561027874778 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r99 | 25.8125 | 26.0625 | 0.48733971724044822 | 0.0096852300242130755 | 0.10000000000000001 | 0.040233256656247753 | PASS |
| radial_profile_L2 | null | null | null | 0.0045289386677212573 | 0.14999999999999999 | 0 | PASS |

### validation-hybrid-volume: step/middle_vs_fine

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 487.07775221479494 | 485.85691882733988 | 0.24317279958228497 | 0.0025064445705101461 | 0.10000000000000001 | 0.0010638983890221478 | PASS |
| r_K_ratio | 1.3411611329909079 | 1.3279515421143209 | 0.0011447918097936809 | 0.0098493689920229523 | 0.10000000000000001 | 0.0018189845251703042 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r99 | 25.8125 | 26.1875 | 0.59773879468097657 | 0.014527845036319613 | 0.10000000000000001 | 0.049347462332790734 | PASS |
| radial_profile_L2 | null | null | null | 0.0013084190740619605 | 0.14999999999999999 | 0 | PASS |

### validation-hybrid-volume: exchange/coarse_vs_fine

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 485.70751326384465 | 485.71918987181817 | 0.58857217199244261 | 2.4040410441796409e-05 | 0.10000000000000001 | 0.0025823098557558577 | PASS |
| r_K_ratio | 1.3275528071184484 | 1.3275877148472155 | 0.0027321546374173091 | 2.6294794888724971e-05 | 0.10000000000000001 | 0.0043856798020516022 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 18.9375 | 0.0625 | 0.0032894736842105261 | 0.10000000000000001 | 0.007009868421052631 | PASS |
| r99 | 24.625 | 25.375 | 0.7772815877574013 | 0.030456852791878174 | 0.10000000000000001 | 0.067264449279635416 | PASS |
| radial_profile_L2 | null | null | null | 0.001356710942527452 | 0.14999999999999999 | 0 | PASS |

### validation-hybrid-volume: exchange/middle_vs_fine

| Metric | Reference mean | Candidate mean | Paired SE | Error | Margin | Smoke allowance | Smoke result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 485.70751326384465 | 485.85691882733988 | 0.6742145635088862 | 0.00030760397855751588 | 0.10000000000000001 | 0.0029580584932334952 | PASS |
| r_K_ratio | 1.3275528071184484 | 1.3279515421143209 | 0.0032461253439093909 | 0.00030035339742001758 | 0.10000000000000001 | 0.0052107103165906013 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r99 | 24.625 | 26.1875 | 1.0245679983941198 | 0.063451776649746189 | 0.10000000000000001 | 0.088664138256969297 | PASS |
| radial_profile_L2 | null | null | null | 0.0022897168552985488 | 0.14999999999999999 | 0 | PASS |

### validation-hybrid-volume: step refinement order

| Metric | Coarse error | Middle error | Sampling allowance | Result |
| --- | ---: | ---: | ---: | --- |
| total_mass | 0.0097367487538173453 | 0.0025064445705101461 | 0.0019109650735822945 | PASS |
| r_K_ratio | 0.034573019170252385 | 0.0098493689920229523 | 0.0035043289988511456 | PASS |
| active_fraction | 0 | 0 | 0 | PASS |
| r50 | 0 | 0 | 0 | PASS |
| r90 | 0 | 0 | 0 | PASS |
| r99 | 0.0096852300242130755 | 0.014527845036319613 | 0.041876101214552701 | PASS |
| radial_profile_L2 | 0.0045289386677212573 | 0.0013084190740619605 | 0.0012330365761312973 | PASS |

### validation-hybrid-volume: exchange refinement order

| Metric | Coarse error | Middle error | Sampling allowance | Result |
| --- | ---: | ---: | ---: | --- |
| total_mass | 2.4040410441796409e-05 | 0.00030760397855751588 | 0.0026007467519915807 | PASS |
| r_K_ratio | 2.6294794888724971e-05 | 0.00030035339742001758 | 0.0045504776924638734 | PASS |
| active_fraction | 0 | 0 | 0 | PASS |
| r50 | 0 | 0 | 0 | PASS |
| r90 | 0.0032894736842105261 | 0 | 0.0070098684210526301 | PASS |
| r99 | 0.030456852791878174 | 0.063451776649746189 | 0.090837068674312224 | PASS |
| radial_profile_L2 | 0.001356710942527452 | 0.0022897168552985488 | 0.003676201764210134 | PASS |

Mixed coverage: `{"abm_mass": 5.5, "pde_mass": 476.8351985179054, "to_abm": 8.9375, "to_pde": 256.5}`.

## Statistical trajectories

### validation-statistical-regular

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 255.99952363215377 | 1.8608118993179712e-06 | 0 | null | -1.051936204319238 | 0.30948091139523803 | 2.4045068684167232e-72 | 2.4041233845721485e-72 | PASS |
| 4 | r_K_ratio | 1 | 0.99999627838045557 | 3.7216195444764177e-06 | 0 | null | -1.0519350021930869 | 0.30948144517619536 | 7.8797163836281917e-68 | 7.8772031868028145e-68 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13.0625 | 13.0625 | 0 | 0.0047846889952153108 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 4 | r90 | 17.9375 | 17.6875 | 0.013937282229965157 | 0.0069686411149825784 | 2 | -2.2360679774997898 | 0.040968955955836141 | 1.1148111900685798e-16 | 3.8956480726576652e-17 | PASS |
| 4 | r99 | 19.25 | 19.25 | 0 | 0.00974025974025974 | 0 | 0 | 1 | 1.6302193524528563e-20 | 9.3136886168229159e-19 | PASS |
| 4 | radial_profile_L2 | null | null | 0.018038965598942157 | 0.026445474852675486 | 0.68211917915768328 | -6.9561301233640789 | 4.6059224602554363e-06 | 1.1097248239271743e-18 | 5.1697635186980713e-18 | PASS |
| 8 | total_mass | 256 | 255.99952363105518 | 1.8608161906907839e-06 | 0 | null | -1.0519386328313456 | 0.30947983306326909 | 2.4045067812667786e-72 | 2.4041232965518539e-72 | PASS |
| 8 | r_K_ratio | 1 | 0.99999627837972738 | 3.721620272581494e-06 | 0 | null | -1.0519352085026517 | 0.30948135356836082 | 7.8797163269185298e-68 | 7.8772031296198129e-68 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 8 | r90 | 18.0625 | 17.9375 | 0.006920415224913495 | 0.0034602076124567475 | 2 | -1.4638501094227998 | 0.16387561365565406 | 2.2451366252918472e-19 | 3.9072803388172492e-18 | PASS |
| 8 | r99 | 20.0625 | 20 | 0.0031152647975077881 | 0.0031152647975077881 | 1 | -1 | 0.33317013591547739 | 6.5149803197211003e-25 | 5.0766352092358305e-21 | PASS |
| 8 | radial_profile_L2 | null | null | 0.076347291770587397 | 0.031674679515173401 | 2.4103571982161358 | -0.25602421236708434 | 0.80140998523340368 | 4.1542646456033906e-17 | 2.940075560798739e-14 | PASS |
| 12 | total_mass | 256 | 255.99952362981946 | 1.8608210177253892e-06 | 0 | null | -1.0519413621231255 | 0.30947862117940073 | 2.40450676396229e-72 | 2.4041232782553873e-72 | PASS |
| 12 | r_K_ratio | 1 | 0.9999962783799754 | 3.7216200246062425e-06 | 0 | null | -1.0519351409582454 | 0.30948138356016608 | 7.8797160406440131e-68 | 7.8772028436040692e-68 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 12 | r90 | 18.6875 | 18 | 0.03678929765886288 | 0 | null | -5.7445626465380286 | 3.8761814875216441e-05 | 1.2411344634416193e-17 | 3.0051609001276733e-16 | PASS |
| 12 | r99 | 20.5625 | 20 | 0.027355623100303952 | 0.015197568389057751 | 1.8 | -4.3915503282683988 | 0.00052573078993947925 | 7.4736080236482832e-20 | 4.9956120358567663e-17 | PASS |
| 12 | radial_profile_L2 | null | null | 0.086013650856367627 | 0.033001760468462815 | 2.6063352268302205 | 0.83491828510482491 | 0.41686386867296688 | 5.7822343406166553e-19 | 1.0101692911102647e-15 | PASS |
| 16 | total_mass | 256 | 256.0039221408303 | 1.5320862618309339e-05 | 0 | null | 8.2445046165417484 | 5.931716341289307e-07 | 5.0325713511847685e-72 | 5.0391845497421992e-72 | PASS |
| 16 | r_K_ratio | 1 | 1.0000306417564719 | 3.0641756472001014e-05 | 0 | null | 8.2445130770941919 | 5.9316405906030129e-07 | 1.6479904258231811e-67 | 1.652324456689296e-67 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.0047619047619047623 | 1 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 16 | r90 | 19.0625 | 18 | 0.055737704918032788 | 0.0032786885245901639 | 17 | -17 | 3.2769117377966068e-11 | 2.2623006381054966e-21 | 5.3331248280097439e-21 | PASS |
| 16 | r99 | 21.1875 | 20.3125 | 0.041297935103244837 | 0.014749262536873156 | 2.7999999999999998 | -4.8692584054817667 | 0.00020423834070405623 | 2.6523458183220211e-16 | 6.1694084828562171e-16 | PASS |
| 16 | radial_profile_L2 | null | null | 0.10482530053006356 | 0.029422120612765671 | 3.5628057511456857 | 1.4430400601137259 | 0.16956690828673149 | 2.3450427511053387e-19 | 2.3043438284859061e-15 | PASS |
| 20 | total_mass | 256 | 258.86102431150124 | 0.011175876216801689 | 0 | null | 47.823002382253392 | 8.179980011291618e-18 | 9.7765779169300254e-41 | 2.5488859870111762e-40 | PASS |
| 20 | r_K_ratio | 1 | 1.0223517162604097 | 0.022351716260409549 | 0 | null | 47.822961258902559 | 8.1800848752947955e-18 | 2.0279595428794984e-36 | 1.381084039653937e-35 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0 | null | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 20 | r90 | 19.375 | 18 | 0.070967741935483872 | 0 | null | -11 | 1.4062516106729139e-08 | 1.8615816799604571e-16 | 6.2976132309167012e-17 | PASS |
| 20 | r99 | 21.9375 | 20.8125 | 0.05128205128205128 | 0 | null | -7.2681556777852343 | 2.7487748593004314e-06 | 1.4996759808931379e-17 | 4.7769350818866047e-17 | PASS |
| 20 | radial_profile_L2 | null | null | 0.11212977992124216 | 0.040888641583758337 | 2.7423209864174609 | 2.5808673065357017 | 0.020878654508262644 | 4.9338502375754997e-19 | 9.5064448425540049e-15 | PASS |
| 24 | total_mass | 284.625 | 295.54288516325516 | 0.038358841153290094 | 0.0083443126921387799 | 4.5970042792653443 | 15.085522450392661 | 1.7915079529863427e-10 | 1.6185834028470003e-27 | 2.3745848934894635e-22 | PASS |
| 24 | r_K_ratio | 1.2236328125 | 1.3089045182921968 | 0.069687331788674811 | 0.015163607342378291 | 4.5956961437478698 | 15.081266344658202 | 1.7986632909266017e-10 | 3.915099904104243e-24 | 8.7146401729868153e-18 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.5625 | 13.0625 | 0.036866359447004608 | 0.018433179723502304 | 2 | -3.8729833462074166 | 0.0015017735386984323 | 9.6072681659644597e-15 | 7.6467248765302002e-14 | PASS |
| 24 | r90 | 19.9375 | 18.0625 | 0.094043887147335428 | 0.003134796238244514 | 30 | -21.957751641341996 | 8.1207854767429677e-13 | 3.8511104687127252e-17 | 1.1403537175415299e-20 | PASS |
| 24 | r99 | 22.875 | 21 | 0.081967213114754092 | 0.0054644808743169399 | 15 | -9.3026050941906355 | 1.2823132219716826e-07 | 3.6636365223864341e-16 | 8.7002119456432826e-16 | PASS |
| 24 | radial_profile_L2 | null | null | 0.1549705493726754 | 0.031998880089473554 | 4.842999159325422 | 6.2088751769977391 | 1.673116493934861e-05 | 8.4074367434956288e-19 | 1.0745435957428209e-12 | PASS |
| 28 | total_mass | 360.0625 | 352.85433092867186 | 0.020019216306413928 | 0.0034716195105016492 | 5.7665352570625314 | -8.6085890851100029 | 3.4511311655701212e-07 | 5.2178222392982184e-27 | 6.405122629200521e-24 | PASS |
| 28 | r_K_ratio | 1.81298828125 | 1.7548196381605756 | 0.032084401036127468 | 0.0053864799353622404 | 5.9564690523570647 | -8.8990356437955676 | 2.2654371377277729e-07 | 6.4042098163988693e-24 | 2.8410539777923426e-21 | PASS |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.875 | 13.0625 | 0.058558558558558557 | 0.022522522522522521 | 2.6000000000000001 | -8.0622577482985509 | 7.8258249630702749e-07 | 2.2778649613890186e-15 | 2.5893641111623486e-16 | PASS |
| 28 | r90 | 20.5625 | 18.1875 | 0.11550151975683891 | 0.0060790273556231003 | 19 | -19 | 6.625544366525624e-12 | 3.1383809442133001e-14 | 1.7974071832715349e-18 | PASS |
| 28 | r99 | 23.625 | 21.3125 | 0.097883597883597878 | 0 | null | -10.593069184806563 | 2.3293939237721809e-08 | 1.1628481028638976e-14 | 4.8813247618385495e-16 | PASS |
| 28 | radial_profile_L2 | null | null | 0.19639608022598237 | 0.025528598153441336 | 7.6931791963479812 | 10.104738181600291 | 4.3560148948139394e-08 | 1.1925294546933859e-19 | 1.5887187151941179e-11 | PASS |
| 32 | total_mass | 385.5625 | 378.21658316815626 | 0.019052467062651863 | 0.0024315124007132437 | 7.8356446206332899 | -6.2508629222441829 | 1.552976500195345e-05 | 1.32766745996875e-25 | 1.7383624492422328e-22 | PASS |
| 32 | r_K_ratio | 1.8212005047589639 | 1.9237669667558863 | 0.056318050499605525 | 0.0014984898683687194 | 37.583204056570814 | 13.043287361456954 | 1.3722633656846121e-09 | 2.4042894780627876e-25 | 2.3397825189013335e-18 | PASS |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.022421524663677129 | 2.7999999999999998 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 32 | r90 | 20.8125 | 18.375 | 0.11711711711711711 | 0.006006006006006006 | 19.5 | -19.030051422496395 | 6.4760803247711866e-12 | 1.9103560717606373e-13 | 6.0923273587867951e-19 | PASS |
| 32 | r99 | 24.375 | 21.8125 | 0.10512820512820513 | 0.0051282051282051282 | 20.5 | -16.291747991900039 | 6.0154840401030353e-11 | 1.0891629608238304e-16 | 1.894984565766644e-18 | PASS |
| 32 | radial_profile_L2 | null | null | 0.21017492230657914 | 0.026253620259367155 | 8.0055596230234123 | 13.506167789354492 | 8.4537564388173067e-10 | 1.1795942561568228e-20 | 9.5046041434704363e-12 | PASS |
| 36 | total_mass | 404.5 | 390.49257084145336 | 0.034628996683675278 | 0.001854140914709518 | 18.676572211395534 | -10.952989422141025 | 1.4895456776392415e-08 | 1.4823099630965238e-25 | 2.817115659160833e-22 | PASS |
| 36 | r_K_ratio | 1.6452597661040846 | 1.8638847743910563 | 0.13288175690618612 | 0.0025243036803777237 | 52.640955182659482 | 16.836244400209452 | 3.7632573463113974e-11 | 9.5821375511934502e-23 | 1.5870311004069577e-12 | PASS |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0089686098654708519 | 7 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 36 | r90 | 21.3125 | 18.625 | 0.12609970674486803 | 0 | null | -17.854778168554091 | 1.6220160419503157e-11 | 2.1663433333220131e-12 | 5.0379176673465418e-18 | PASS |
| 36 | r99 | 24.625 | 22 | 0.1065989847715736 | 0.0025380710659898475 | 42 | -16.959029914832215 | 3.3919667718492142e-11 | 1.4370664950489062e-17 | 2.2632691711460884e-18 | PASS |
| 36 | radial_profile_L2 | null | null | 0.22528787583791104 | 0.021200415960106734 | 10.626578094592103 | 16.530925563679663 | 4.8870224731012654e-11 | 1.0067641583813752e-20 | 6.1072726191039248e-11 | PASS |
| 40 | total_mass | 421.125 | 409.96644251661564 | 0.02649701984775148 | 0.0014841199168892847 | 17.853691973414946 | -9.4131666753432288 | 1.100605517036285e-07 | 1.0051708933424404e-26 | 1.0210541251396914e-22 | PASS |
| 40 | r_K_ratio | 1.4862818219349081 | 1.649570562355545 | 0.10986391545047657 | 0.0052088165665525443 | 21.091914842221062 | 12.942794502386514 | 1.5273927728310308e-09 | 3.4817355713089573e-22 | 1.4462708949639824e-12 | PASS |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0044843049327354259 | 14 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 40 | r90 | 21.5 | 18.8125 | 0.125 | 0 | null | -22.456017618285021 | 5.8544703070058476e-13 | 2.3582850886443684e-14 | 3.6683425778814834e-19 | PASS |
| 40 | r99 | 25 | 22 | 0.12 | 0.0074999999999999997 | 16 | -23.2379000772445 | 3.5513304925948273e-13 | 2.2385921293450308e-18 | 7.4091209429095752e-20 | PASS |
| 40 | radial_profile_L2 | null | null | 0.23931975701012922 | 0.019135549966259355 | 12.506552329674786 | 15.167363007239185 | 1.6596426621118105e-10 | 7.9124143814230561e-22 | 4.19301227131262e-11 | PASS |
| 44 | total_mass | 433 | 441.50655922322505 | 0.019645633309988586 | 0.0072170900692840644 | 2.7220989514320184 | 6.1204274535427654 | 1.9590294389756984e-05 | 3.9524299841430892e-26 | 3.1226975632568035e-21 | PASS |
| 44 | r_K_ratio | 1.3851122690599202 | 1.4186749723571346 | 0.024231034586094116 | 0.023883569861476831 | 1.0145482742585199 | 3.4856758958242473 | 0.0033196876379717162 | 1.1600118600015103e-21 | 6.7929423477361276e-16 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 14.0625 | 13.0625 | 0.071111111111111111 | 0.0088888888888888889 | 8 | -7.7459666924148332 | 1.2783477850604916e-06 | 1.25956553381152e-13 | 5.7635238603637946e-15 | PASS |
| 44 | r90 | 21.875 | 18.9375 | 0.13428571428571429 | 0.0057142857142857143 | 23.5 | -47 | 1.0595740345741668e-17 | 1.9641369866418272e-18 | 2.4448502950007348e-23 | PASS |
| 44 | r99 | 25.4375 | 22.3125 | 0.12285012285012285 | 0.0024570024570024569 | 50 | -20.189321327181208 | 2.7524505961136397e-12 | 3.0772991463099793e-16 | 3.5800964663009798e-19 | PASS |
| 44 | radial_profile_L2 | null | null | 0.2444336273714921 | 0.022666905945975712 | 10.783722663961065 | 15.812234710545606 | 9.20033334826029e-11 | 1.4422645201058869e-22 | 1.8414881454564072e-11 | PASS |
| 48 | total_mass | 468.375 | 489.72804166911277 | 0.045589627262583969 | 0.0066720042700827327 | 6.832973334116085 | 11.27344080509701 | 1.0100996484623647e-08 | 5.3377976096866409e-25 | 3.656817895487457e-19 | PASS |
| 48 | r_K_ratio | 1.4233832881828432 | 1.3285012486394163 | 0.066659514925566946 | 0.023969344303723984 | 2.7810320583200276 | -8.28438111864814 | 5.585754708074692e-07 | 6.9070416872612281e-19 | 1.4957882213155598e-16 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 14.4375 | 13.1875 | 0.086580086580086577 | 0.030303030303030304 | 2.8571428571428572 | -11.180339887498949 | 1.129731174656377e-08 | 5.7520127976418164e-14 | 2.7873182537487998e-16 | PASS |
| 48 | r90 | 22.375 | 19 | 0.15083798882681565 | 0 | null | -21.804467033355703 | 8.9935030956328787e-13 | 3.0419385484799993e-12 | 6.5100670518587811e-18 | PASS |
| 48 | r99 | 26.1875 | 22.8125 | 0.12887828162291171 | 0.0071599045346062056 | 18 | -27 | 3.9304962308691221e-14 | 3.8161521502635611e-17 | 3.2641019302848734e-21 | PASS |
| 48 | radial_profile_L2 | null | null | 0.27458671035633075 | 0.017773720507827693 | 15.449028256937005 | 15.456422359839008 | 1.2705593178536584e-10 | 3.043978705295418e-21 | 5.9937116945256635e-08 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1796875, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

### validation-statistical-regular-v14

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 255.99987439227883 | 4.9065516088270256e-07 | 0 | null | -1.0352175095070719 | 0.31696944719393905 | 6.3310518845588904e-81 | 6.3307856304481416e-81 | PASS |
| 4 | r_K_ratio | 1 | 0.9999990186955926 | 9.8130440741306391e-07 | 0 | null | -1.0352112711393178 | 0.31697226570925685 | 2.0746026777977294e-76 | 2.0744281865575098e-76 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13.0625 | 13.0625 | 0 | 0.0047846889952153108 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 4 | r90 | 17.9375 | 17.6875 | 0.013937282229965157 | 0.0069686411149825784 | 2 | -2.2360679774997898 | 0.040968955955836141 | 1.1148111900685798e-16 | 3.8956480726576652e-17 | PASS |
| 4 | r99 | 19.25 | 19.375 | 0.0064935064935064939 | 0.00974025974025974 | 0.66666666666666663 | 1 | 0.33317013591547739 | 1.2617854706008803e-18 | 7.1507823999893397e-17 | PASS |
| 4 | radial_profile_L2 | null | null | 0.01390051111786134 | 0.026445474852675486 | 0.52562909894034393 | -7.7751682619068863 | 1.2210780923911863e-06 | 1.7063425300474121e-19 | 5.5915970512964533e-19 | PASS |
| 8 | total_mass | 256 | 255.99987439079101 | 4.9066097267125297e-07 | 0 | null | -1.0352297764902778 | 0.31696390498294125 | 6.331051438669818e-81 | 6.3307851814239781e-81 | PASS |
| 8 | r_K_ratio | 1 | 0.99999901869484664 | 9.8130515333721968e-07 | 0 | null | -1.0352120739000901 | 0.3169719030182524 | 2.0746022011006527e-76 | 2.0744277097679295e-76 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 8 | r90 | 18.0625 | 18 | 0.0034602076124567475 | 0.0034602076124567475 | 1 | -1 | 0.33317013591547739 | 1.4301148022486204e-22 | 1.9701065876166452e-19 | PASS |
| 8 | r99 | 20.0625 | 20 | 0.0031152647975077881 | 0.0031152647975077881 | 1 | -1 | 0.33317013591547739 | 6.5149803197211003e-25 | 5.0766352092358305e-21 | PASS |
| 8 | radial_profile_L2 | null | null | 0.067252396223438593 | 0.031674679515173401 | 2.1232226261743889 | -1.1211565613575198 | 0.2798500299122274 | 1.1233733444798528e-17 | 3.6201109868020346e-15 | PASS |
| 12 | total_mass | 256 | 255.99987438921451 | 4.9066713076612034e-07 | 0 | null | -1.035242772400685 | 0.31695803351978358 | 6.3310511518897591e-81 | 6.3307848913142876e-81 | PASS |
| 12 | r_K_ratio | 1 | 0.99999901869489838 | 9.8130510160082673e-07 | 0 | null | -1.0352120282386914 | 0.31697192364827376 | 2.0746019330439565e-76 | 2.0744274417429762e-76 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 12 | r90 | 18.6875 | 18 | 0.03678929765886288 | 0 | null | -5.7445626465380286 | 3.8761814875216441e-05 | 1.2411344634416193e-17 | 3.0051609001276733e-16 | PASS |
| 12 | r99 | 20.5625 | 20.125 | 0.021276595744680851 | 0.015197568389057751 | 1.3999999999999999 | -3.415650255319866 | 0.0038327022878305922 | 3.2453744019125243e-19 | 3.5829058514479523e-17 | PASS |
| 12 | radial_profile_L2 | null | null | 0.072893691384749507 | 0.033001760468462815 | 2.2087819058746359 | -0.25736842332843862 | 0.80039152143350745 | 4.3685958447085025e-20 | 2.3969087135186269e-17 | PASS |
| 16 | total_mass | 256 | 256.00565975007765 | 2.2108398740866564e-05 | 0 | null | 29.0549503061726 | 1.3336686824631748e-14 | 7.6725742637134177e-78 | 7.6871276199674024e-78 | PASS |
| 16 | r_K_ratio | 1 | 1.0000442168399475 | 4.4216839947416875e-05 | 0 | null | 29.054978697547448 | 1.3336494557056482e-14 | 2.5117676850929536e-73 | 2.5213053845686074e-73 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.0047619047619047623 | 1 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 16 | r90 | 19.0625 | 18 | 0.055737704918032788 | 0.0032786885245901639 | 17 | -17 | 3.2769117377966068e-11 | 2.2623006381054966e-21 | 5.3331248280097439e-21 | PASS |
| 16 | r99 | 21.1875 | 20.8125 | 0.017699115044247787 | 0.014749262536873156 | 1.2 | -2.4227185592617446 | 0.028528068449741654 | 4.6950640203719188e-18 | 3.0093626661365215e-16 | PASS |
| 16 | radial_profile_L2 | null | null | 0.08775376074242433 | 0.029422120612765671 | 2.9825776971476938 | 0.20757506469588111 | 0.83835262262231791 | 5.9173726605820711e-20 | 1.2294617205940524e-16 | PASS |
| 20 | total_mass | 256 | 259.38930481728698 | 0.013239471942527239 | 0 | null | 48.255029981106148 | 7.1535576417554098e-18 | 9.9612836898674566e-40 | 3.1001936901933769e-39 | PASS |
| 20 | r_K_ratio | 1 | 1.0264789234962126 | 0.026478923496212517 | 0 | null | 48.255010551841188 | 7.1536005859322785e-18 | 1.9075543162085774e-35 | 1.8535898820419411e-34 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0 | null | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 20 | r90 | 19.375 | 18.125 | 0.064516129032258063 | 0 | null | -11.180339887498949 | 1.129731174656377e-08 | 6.9886341424657286e-17 | 1.1645945213682297e-17 | PASS |
| 20 | r99 | 21.9375 | 21 | 0.042735042735042736 | 0 | null | -6.5361700549589266 | 9.4202794798815216e-06 | 3.6683425778815099e-19 | 5.1861867803112176e-17 | PASS |
| 20 | radial_profile_L2 | null | null | 0.092421060738694855 | 0.040888641583758337 | 2.2603113519771729 | 1.0369297345499009 | 0.31619654317416779 | 1.1226962030547471e-19 | 3.5490613354533042e-16 | PASS |
| 24 | total_mass | 284.625 | 298.98754644036325 | 0.050461296233160362 | 0.0083443126921387799 | 6.047387974889797 | 19.730067038543012 | 3.8414531453545918e-12 | 1.6655638915222936e-27 | 4.0384969275713948e-22 | PASS |
| 24 | r_K_ratio | 1.2236328125 | 1.3358250453136638 | 0.091687826337742723 | 0.015163607342378291 | 6.0465708632206123 | 19.727446316178806 | 3.8488472047269843e-12 | 2.9675877127982947e-24 | 2.7669535033326392e-17 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.5625 | 13.0625 | 0.036866359447004608 | 0.018433179723502304 | 2 | -3.8729833462074166 | 0.0015017735386984323 | 9.6072681659644597e-15 | 7.6467248765302002e-14 | PASS |
| 24 | r90 | 19.9375 | 18.25 | 0.084639498432601878 | 0.003134796238244514 | 27 | -14.100290132411525 | 4.6333644836412128e-10 | 7.7245411315955695e-15 | 9.481589102890017e-19 | PASS |
| 24 | r99 | 22.875 | 21.25 | 0.071038251366120214 | 0.0054644808743169399 | 13 | -6.7890285822722154 | 6.1053729318667233e-06 | 1.2337518345797697e-14 | 7.416708381213179e-15 | PASS |
| 24 | radial_profile_L2 | null | null | 0.1311571822925125 | 0.031998880089473554 | 4.0988053933693243 | 4.1375477785895995 | 0.00087694673719733397 | 3.3174571066196259e-19 | 3.9701746118309766e-14 | PASS |
| 28 | total_mass | 360.0625 | 355.79428911997218 | 0.011854083332831953 | 0.0034716195105016492 | 3.4145687040222441 | -4.9981373855468707 | 0.00015894216488354686 | 4.9346299105874915e-27 | 1.1189772562667494e-23 | PASS |
| 28 | r_K_ratio | 1.81298828125 | 1.7783820031380952 | 0.019087976723183724 | 0.0053864799353622404 | 3.543682878659058 | -5.1899130998290026 | 0.00010985368152231689 | 4.8401929701778768e-24 | 5.9660093150193663e-21 | PASS |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.875 | 13.0625 | 0.058558558558558557 | 0.022522522522522521 | 2.6000000000000001 | -8.0622577482985509 | 7.8258249630702749e-07 | 2.2778649613890186e-15 | 2.5893641111623486e-16 | PASS |
| 28 | r90 | 20.5625 | 18.75 | 0.088145896656534953 | 0.0060790273556231003 | 14.5 | -11.066875266096559 | 1.2961577489987425e-08 | 1.5680800107416445e-13 | 1.9477056267575937e-16 | PASS |
| 28 | r99 | 23.625 | 21.9375 | 0.071428571428571425 | 0 | null | -8.5098307000225546 | 3.9910679004152207e-07 | 1.2402986697836266e-16 | 5.9545457575235637e-16 | PASS |
| 28 | radial_profile_L2 | null | null | 0.16812366489946665 | 0.025528598153441336 | 6.5856990614583761 | 8.2179631045461328 | 6.1744432441080883e-07 | 2.0668183461795087e-20 | 1.1633368762049021e-13 | PASS |
| 32 | total_mass | 385.5625 | 379.14312167383497 | 0.016649384538602752 | 0.0024315124007132437 | 6.8473368812426925 | -5.396509588004732 | 7.4173402428140978e-05 | 1.4730519212317234e-25 | 2.2336216332001333e-22 | PASS |
| 32 | r_K_ratio | 1.8212005047589639 | 1.9391831914335063 | 0.064782920038866071 | 0.0014984898683687194 | 43.232137504799965 | 14.726328885627707 | 2.5168461149834619e-10 | 2.6740318924724625e-25 | 4.4405753482401013e-18 | PASS |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.022421524663677129 | 2.7999999999999998 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 32 | r90 | 20.8125 | 18.9375 | 0.090090090090090086 | 0.006006006006006006 | 15 | -21.957751641341996 | 8.1207854767429677e-13 | 2.5303126792687863e-18 | 2.9422487852212645e-20 | PASS |
| 32 | r99 | 24.375 | 22 | 0.097435897435897437 | 0.0051282051282051282 | 19 | -15.343884208657718 | 1.4090734244550085e-10 | 8.3953741872725587e-18 | 3.7020237911275972e-18 | PASS |
| 32 | radial_profile_L2 | null | null | 0.17938424840949052 | 0.026253620259367155 | 6.8327433183424358 | 10.818675251643851 | 1.7576228289281234e-08 | 2.8179852830823343e-21 | 5.7756288571017507e-14 | PASS |
| 36 | total_mass | 404.5 | 389.43826184084367 | 0.037235446623377903 | 0.001854140914709518 | 20.082317545541816 | -11.824525086757328 | 5.2854372155451087e-09 | 1.4574793360891194e-25 | 2.4656946513950525e-22 | PASS |
| 36 | r_K_ratio | 1.6452597661040846 | 1.8931472248902017 | 0.15066767199510719 | 0.0025243036803777237 | 59.686824991104899 | 18.653271937262037 | 8.6419882044172036e-12 | 9.373935059396823e-23 | 6.9902926261365534e-12 | PASS |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0089686098654708519 | 7 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 36 | r90 | 21.3125 | 19 | 0.10850439882697947 | 0 | null | -19.322619811082461 | 5.195040620846464e-12 | 7.7864380968352111e-16 | 1.5335258890741899e-18 | PASS |
| 36 | r99 | 24.625 | 22.25 | 0.096446700507614211 | 0.0025380710659898475 | 38 | -11.783299785974805 | 5.5430657455111653e-09 | 3.3901457391067851e-15 | 5.0500172350014425e-17 | PASS |
| 36 | radial_profile_L2 | null | null | 0.19113782511778182 | 0.021200415960106734 | 9.0157582510385588 | 12.937967455874144 | 1.5352977834313781e-09 | 2.3205167683510447e-21 | 1.873559415950491e-13 | PASS |
| 40 | total_mass | 421.125 | 406.62557094179493 | 0.034430226318088571 | 0.0014841199168892847 | 23.199086493128078 | -12.123072892031722 | 3.7595679941645351e-09 | 1.454436824969024e-26 | 8.6049179868181217e-23 | PASS |
| 40 | r_K_ratio | 1.4862818219349081 | 1.6995080743478872 | 0.1434628677187155 | 0.0052088165665525443 | 27.542315204558335 | 16.515053506398107 | 4.9544479197942529e-11 | 1.9779308949787462e-22 | 1.6280039013798181e-11 | PASS |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0044843049327354259 | 14 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 40 | r90 | 21.5 | 19.0625 | 0.11337209302325581 | 0 | null | -19.030051422496395 | 6.4760803247711866e-12 | 6.3605411260372649e-15 | 2.3244786420838592e-18 | PASS |
| 40 | r99 | 25 | 22.6875 | 0.092499999999999999 | 0.0074999999999999997 | 12.333333333333334 | -15.363413772941893 | 1.3839300528290396e-10 | 5.4380328870719136e-17 | 4.8751207713627678e-19 | PASS |
| 40 | radial_profile_L2 | null | null | 0.20294210478022673 | 0.019135549966259355 | 10.605501547541785 | 12.331063058139474 | 2.9774061335968384e-09 | 2.6639476681283848e-22 | 9.5858061539418275e-14 | PASS |
| 44 | total_mass | 433 | 439.12414900465899 | 0.014143531188588849 | 0.0072170900692840644 | 1.9597276814908708 | 4.3124362386774742 | 0.00061619879125071012 | 5.6783101501054603e-26 | 3.3785584158268306e-21 | PASS |
| 44 | r_K_ratio | 1.3851122690599202 | 1.4826505247964381 | 0.070419025168780944 | 0.023883569861476831 | 2.948429634983663 | 10.150700193299443 | 4.1027820504537489e-08 | 1.7459206768544517e-22 | 6.5239376359162711e-15 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 14.0625 | 13.0625 | 0.071111111111111111 | 0.0088888888888888889 | 8 | -7.7459666924148332 | 1.2783477850604916e-06 | 1.25956553381152e-13 | 5.7635238603637946e-15 | PASS |
| 44 | r90 | 21.875 | 19.1875 | 0.12285714285714286 | 0.0057142857142857143 | 21.5 | -22.456017618285021 | 5.8544703070058476e-13 | 4.6040517988789706e-14 | 9.7976516786557512e-20 | PASS |
| 44 | r99 | 25.4375 | 22.9375 | 0.098280098280098274 | 0.0024570024570024569 | 40 | -15.811388300841896 | 9.2073297610083179e-11 | 1.32240574860815e-17 | 2.0753132400270518e-18 | PASS |
| 44 | radial_profile_L2 | null | null | 0.20578903863170117 | 0.022666905945975712 | 9.0788323347782285 | 13.563001661788283 | 7.9734125142345473e-10 | 3.9888272112335183e-23 | 2.1411534627969286e-14 | PASS |
| 48 | total_mass | 468.375 | 495.63025101823672 | 0.058191088376272695 | 0.0066720042700827327 | 8.7216803258357523 | 13.932439671958022 | 5.478965023639955e-10 | 6.2262438154028136e-25 | 9.8983062404198524e-19 | PASS |
| 48 | r_K_ratio | 1.4233832881828432 | 1.4168010867127527 | 0.00462433521928831 | 0.023969344303723984 | 0.19292706386506586 | -0.56934483893448984 | 0.5775492092836918 | 5.7033749883260893e-20 | 1.7002683649562875e-15 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 14.4375 | 13.1875 | 0.086580086580086577 | 0.030303030303030304 | 2.8571428571428572 | -11.180339887498949 | 1.129731174656377e-08 | 5.7520127976418164e-14 | 2.7873182537487998e-16 | PASS |
| 48 | r90 | 22.375 | 19.5 | 0.12849162011173185 | 0 | null | -15.998991903725805 | 7.7865353236989104e-11 | 1.0061554522861548e-11 | 4.8264550809103502e-17 | PASS |
| 48 | r99 | 26.1875 | 23.1875 | 0.11455847255369929 | 0.0071599045346062056 | 16 | -23.2379000772445 | 3.5513304925948273e-13 | 1.7995348202025614e-17 | 8.6169881800138773e-21 | PASS |
| 48 | radial_profile_L2 | null | null | 0.23341720836842744 | 0.017773720507827693 | 13.1327151378142 | 12.748408916387797 | 1.8827621113569947e-09 | 3.1119166670442979e-21 | 6.3693627289909649e-11 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1875, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

### validation-statistical-spatial-birth

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658184 | 6.5213491597845962e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 255.99987439227883 | 4.9065516088270256e-07 | 0 | null | -1.0352175095070719 | 0.31696944719393905 | 6.3310518845588904e-81 | 6.3307856304481416e-81 | PASS |
| 4 | r_K_ratio | 1 | 0.9999990186955926 | 9.8130440741306391e-07 | 0 | null | -1.0352112711393178 | 0.31697226570925685 | 2.0746026777977294e-76 | 2.0744281865575098e-76 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13.0625 | 13.0625 | 0 | 0.0047846889952153108 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 4 | r90 | 17.9375 | 17.6875 | 0.013937282229965157 | 0.0069686411149825784 | 2 | -2.2360679774997898 | 0.040968955955836141 | 1.1148111900685798e-16 | 3.8956480726576652e-17 | PASS |
| 4 | r99 | 19.25 | 19.375 | 0.0064935064935064939 | 0.00974025974025974 | 0.66666666666666663 | 1 | 0.33317013591547739 | 1.2617854706008803e-18 | 7.1507823999893397e-17 | PASS |
| 4 | radial_profile_L2 | null | null | 0.01390051111786134 | 0.026445474852675489 | 0.52562909894034382 | -7.7751682619068863 | 1.2210780923911863e-06 | 1.706342530047388e-19 | 5.5915970512963349e-19 | PASS |
| 8 | total_mass | 256 | 255.99987439079101 | 4.9066097267125297e-07 | 0 | null | -1.0352297764902778 | 0.31696390498294125 | 6.331051438669818e-81 | 6.3307851814239781e-81 | PASS |
| 8 | r_K_ratio | 1 | 0.99999901869484664 | 9.8130515333721968e-07 | 0 | null | -1.0352120739000901 | 0.3169719030182524 | 2.0746022011006527e-76 | 2.0744277097679295e-76 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 8 | r90 | 18.0625 | 18 | 0.0034602076124567475 | 0.0034602076124567475 | 1 | -1 | 0.33317013591547739 | 1.4301148022486204e-22 | 1.9701065876166452e-19 | PASS |
| 8 | r99 | 20.0625 | 20 | 0.0031152647975077881 | 0.0031152647975077881 | 1 | -1 | 0.33317013591547739 | 6.5149803197211003e-25 | 5.0766352092358305e-21 | PASS |
| 8 | radial_profile_L2 | null | null | 0.067252396223438676 | 0.031674679515173401 | 2.1232226261743916 | -1.121156561357521 | 0.27985002991222763 | 1.1233733444800052e-17 | 3.6201109868025876e-15 | PASS |
| 12 | total_mass | 256 | 255.99987438921454 | 4.9066713065509804e-07 | 0 | null | -1.0352427723769593 | 0.31695803353050223 | 6.3310511325783142e-81 | 6.3307848720036545e-81 | PASS |
| 12 | r_K_ratio | 1 | 0.99999901869489871 | 9.8130510128163762e-07 | 0 | null | -1.0352120281216213 | 0.31697192370116678 | 2.0746019264410209e-76 | 2.0744274351406548e-76 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 12 | r90 | 18.6875 | 18 | 0.03678929765886288 | 0 | null | -5.7445626465380286 | 3.8761814875216441e-05 | 1.2411344634416193e-17 | 3.0051609001276733e-16 | PASS |
| 12 | r99 | 20.5625 | 20.125 | 0.021276595744680851 | 0.015197568389057751 | 1.3999999999999999 | -3.415650255319866 | 0.0038327022878305922 | 3.2453744019125243e-19 | 3.5829058514479523e-17 | PASS |
| 12 | radial_profile_L2 | null | null | 0.0728936913847439 | 0.033001760468462815 | 2.2087819058744658 | -0.25736842332884102 | 0.80039152143320558 | 4.3685958447071682e-20 | 2.3969087135167167e-17 | PASS |
| 16 | total_mass | 256 | 256.00565990929204 | 2.2109020672067548e-05 | 0 | null | 29.055621071987371 | 1.333214515984822e-14 | 7.673154674705337e-78 | 7.6877095417039995e-78 | PASS |
| 16 | r_K_ratio | 1 | 1.0000442180838114 | 4.4218083811317643e-05 | 0 | null | 29.055649465249196 | 1.3331952949265112e-14 | 2.5119576257366105e-73 | 2.5214963152902277e-73 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.0047619047619047623 | 1 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 16 | r90 | 19.0625 | 18 | 0.055737704918032788 | 0.0032786885245901639 | 17 | -17 | 3.2769117377966068e-11 | 2.2623006381054966e-21 | 5.3331248280097439e-21 | PASS |
| 16 | r99 | 21.1875 | 20.8125 | 0.017699115044247787 | 0.014749262536873156 | 1.2 | -2.4227185592617446 | 0.028528068449741654 | 4.6950640203719188e-18 | 3.0093626661365215e-16 | PASS |
| 16 | radial_profile_L2 | null | null | 0.087753521628794323 | 0.029422120612765675 | 2.9825695701457975 | 0.2075579969651104 | 0.8383657119697715 | 5.917506718682807e-20 | 1.2294627722771932e-16 | PASS |
| 20 | total_mass | 256 | 259.38935697918242 | 0.013239675699931355 | 0 | null | 48.255478718939976 | 7.1525658805689808e-18 | 9.9621099804360735e-40 | 3.100505078257499e-39 | PASS |
| 20 | r_K_ratio | 1 | 1.0264793310107851 | 0.026479331010785034 | 0 | null | 48.255459289196338 | 7.1526088194555449e-18 | 1.9076976227082329e-35 | 1.8537942546079051e-34 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0 | null | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 20 | r90 | 19.375 | 18.125 | 0.064516129032258063 | 0 | null | -11.180339887498949 | 1.129731174656377e-08 | 6.9886341424657286e-17 | 1.1645945213682297e-17 | PASS |
| 20 | r99 | 21.9375 | 21 | 0.042735042735042736 | 0 | null | -6.5361700549589266 | 9.4202794798815216e-06 | 3.6683425778815099e-19 | 5.1861867803112176e-17 | PASS |
| 20 | radial_profile_L2 | null | null | 0.09226338129301527 | 0.040888641583758344 | 2.2564550378622465 | 1.0231584195297534 | 0.32245164838009688 | 1.1251404688138574e-19 | 3.5060034693768511e-16 | PASS |
| 24 | total_mass | 284.625 | 298.98832661654501 | 0.050464037300114242 | 0.0083443126921387799 | 6.0477164701242172 | 19.731155400039221 | 3.8383869057078653e-12 | 1.6655364014943582e-27 | 4.0388732636059276e-22 | PASS |
| 24 | r_K_ratio | 1.2236328125 | 1.3358311402728638 | 0.091692807373832841 | 0.015163607342378291 | 6.0468993494427652 | 19.728534634782985 | 3.8457748046898245e-12 | 2.9673416993492329e-24 | 2.7676270594412945e-17 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.5625 | 13.0625 | 0.036866359447004608 | 0.018433179723502304 | 2 | -3.8729833462074166 | 0.0015017735386984323 | 9.6072681659644597e-15 | 7.6467248765302002e-14 | PASS |
| 24 | r90 | 19.9375 | 18.375 | 0.078369905956112859 | 0.003134796238244514 | 25 | -12.198750911856663 | 3.4523439761073071e-09 | 1.4966769468733605e-14 | 2.8116770770607522e-18 | PASS |
| 24 | r99 | 22.875 | 21.375 | 0.065573770491803282 | 0.0054644808743169399 | 12 | -6.7082039324993694 | 7.0065583016648931e-06 | 3.4512116296941117e-15 | 3.5841850263120592e-15 | PASS |
| 24 | radial_profile_L2 | null | null | 0.12927642609054221 | 0.031998880089473547 | 4.0400297050730032 | 3.9655118418184698 | 0.0012436052777004098 | 3.3218560868432211e-19 | 3.3097029976515684e-14 | PASS |
| 28 | total_mass | 360.0625 | 355.79786995269785 | 0.011844138301828507 | 0.0034716195105016492 | 3.4117040378417016 | -4.9934353141941079 | 0.00016039731723934463 | 4.9448840042546594e-27 | 1.120423960280407e-23 | PASS |
| 28 | r_K_ratio | 1.81298828125 | 1.7784099485173064 | 0.019072562735404398 | 0.0053864799353622404 | 3.5408212718278262 | -5.1851935898516457 | 0.00011085077445508441 | 4.8490283073764659e-24 | 5.9749820292080538e-21 | PASS |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.875 | 13.0625 | 0.058558558558558557 | 0.022522522522522521 | 2.6000000000000001 | -8.0622577482985509 | 7.8258249630702749e-07 | 2.2778649613890186e-15 | 2.5893641111623486e-16 | PASS |
| 28 | r90 | 20.5625 | 18.8125 | 0.085106382978723402 | 0.0060790273556231003 | 14 | -12.124355652982143 | 3.7541255191393749e-09 | 1.4192741870443792e-14 | 4.4765474968068944e-17 | PASS |
| 28 | r99 | 23.625 | 21.9375 | 0.071428571428571425 | 0 | null | -8.5098307000225546 | 3.9910679004152207e-07 | 1.2402986697836266e-16 | 5.9545457575235637e-16 | PASS |
| 28 | radial_profile_L2 | null | null | 0.16405669091957084 | 0.025528598153441332 | 6.4263885519094011 | 7.8313246944306272 | 1.1183966370800409e-06 | 1.966819189326734e-20 | 7.1409093913159245e-14 | PASS |
| 32 | total_mass | 385.5625 | 379.16568384801553 | 0.016590866985208501 | 0.0024315124007132437 | 6.8232705621167495 | -5.3686155293505395 | 7.8187578527019736e-05 | 1.5288902634240917e-25 | 2.2786753490476405e-22 | PASS |
| 32 | r_K_ratio | 1.8212005047589639 | 1.9393568638136029 | 0.064878281521384168 | 0.0014984898683687194 | 43.295775894709074 | 14.741142630664925 | 2.4814496209643969e-10 | 2.730946823123324e-25 | 4.471149443760827e-18 | PASS |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.022421524663677129 | 2.7999999999999998 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 32 | r90 | 20.8125 | 18.9375 | 0.090090090090090086 | 0.006006006006006006 | 15 | -21.957751641341996 | 8.1207854767429677e-13 | 2.5303126792687863e-18 | 2.9422487852212645e-20 | PASS |
| 32 | r99 | 24.375 | 22.0625 | 0.094871794871794868 | 0.0051282051282051282 | 18.5 | -13.136324646189436 | 1.2434982186844375e-09 | 1.5609213681963434e-16 | 1.554563354209922e-17 | PASS |
| 32 | radial_profile_L2 | null | null | 0.17487892574276689 | 0.026253620259367158 | 6.6611356458685345 | 10.397644941576296 | 2.9843462695827838e-08 | 2.6100874687192967e-21 | 3.2227504383950966e-14 | PASS |
| 36 | total_mass | 404.5 | 389.45926766739342 | 0.037183516273440319 | 0.001854140914709518 | 20.054309776808811 | -11.805436910281133 | 5.4031100015012418e-09 | 1.4630608635349803e-25 | 2.4759221327288636e-22 | PASS |
| 36 | r_K_ratio | 1.6452597661040846 | 1.893289657624615 | 0.15075424357325423 | 0.0025243036803777237 | 59.721120222229423 | 18.645946032483344 | 8.6910660132474326e-12 | 9.5829111228214293e-23 | 7.1031951951600809e-12 | PASS |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0089686098654708519 | 7 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 36 | r90 | 21.3125 | 19 | 0.10850439882697947 | 0 | null | -19.322619811082461 | 5.195040620846464e-12 | 7.7864380968352111e-16 | 1.5335258890741899e-18 | PASS |
| 36 | r99 | 24.625 | 22.5 | 0.086294416243654817 | 0.0025380710659898475 | 34 | -13.728738502483221 | 6.730952139452569e-10 | 4.5587789692221704e-17 | 1.8183189209815112e-18 | PASS |
| 36 | radial_profile_L2 | null | null | 0.18656080889003238 | 0.021200415960106734 | 8.7998654951434805 | 12.425650837851046 | 2.6806419917659358e-09 | 2.2767547083155321e-21 | 1.0718923015486592e-13 | PASS |
| 40 | total_mass | 421.125 | 406.64803311283043 | 0.034376887829432067 | 0.0014841199168892847 | 23.163147019471324 | -12.095655968352325 | 3.8779148268028664e-09 | 1.473597627386941e-26 | 8.6918298750280041e-23 | PASS |
| 40 | r_K_ratio | 1.4862818219349081 | 1.6996050002855141 | 0.14352808141923748 | 0.0052088165665525443 | 27.554835073455379 | 16.517311941195551 | 4.9447938999950111e-11 | 1.9915408596028513e-22 | 1.6403504379960735e-11 | PASS |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0044843049327354259 | 14 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 40 | r90 | 21.5 | 19.0625 | 0.11337209302325581 | 0 | null | -19.030051422496395 | 6.4760803247711866e-12 | 6.3605411260372649e-15 | 2.3244786420838592e-18 | PASS |
| 40 | r99 | 25 | 22.9375 | 0.082500000000000004 | 0.0074999999999999997 | 11 | -14.379574120909639 | 3.5189622604582387e-10 | 2.4558314713141745e-18 | 7.901250573376388e-19 | PASS |
| 40 | radial_profile_L2 | null | null | 0.19809400512033434 | 0.019135549966259352 | 10.352145899627784 | 11.926715870016155 | 4.7000113633629682e-09 | 2.6669405459628691e-22 | 5.2419008818325543e-14 | PASS |
| 44 | total_mass | 433 | 439.23206584062973 | 0.014392761756650638 | 0.0072170900692840644 | 1.9942610690015123 | 4.3771790308701704 | 0.00054109250360977396 | 5.9224939034279197e-26 | 3.517435842868152e-21 | PASS |
| 44 | r_K_ratio | 1.3851122690599202 | 1.4831263019519321 | 0.070762518736863317 | 0.023883569861476831 | 2.9628116377610789 | 10.19887506362087 | 3.8539975473241908e-08 | 1.7147984527185799e-22 | 6.6649967260819194e-15 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 14.0625 | 13.0625 | 0.071111111111111111 | 0.0088888888888888889 | 8 | -7.7459666924148332 | 1.2783477850604916e-06 | 1.25956553381152e-13 | 5.7635238603637946e-15 | PASS |
| 44 | r90 | 21.875 | 19.1875 | 0.12285714285714286 | 0.0057142857142857143 | 21.5 | -22.456017618285021 | 5.8544703070058476e-13 | 4.6040517988789706e-14 | 9.7976516786557512e-20 | PASS |
| 44 | r99 | 25.4375 | 23 | 0.095823095823095825 | 0.0024570024570024569 | 39 | -15.497028577661013 | 1.2241974350227701e-10 | 5.016695892938304e-18 | 2.6446581536787945e-18 | PASS |
| 44 | radial_profile_L2 | null | null | 0.20022591270254628 | 0.022666905945975712 | 8.8334028993531177 | 13.096159457147007 | 1.2974403191799778e-09 | 4.0060695047060737e-23 | 1.0619373488246274e-14 | PASS |
| 48 | total_mass | 468.375 | 496.1627205009305 | 0.059327932748183647 | 0.0066720042700827327 | 8.8920705602977641 | 14.205603014782996 | 4.1744770622290889e-10 | 5.9167774317339044e-25 | 1.0508057359812269e-18 | PASS |
| 48 | r_K_ratio | 1.4233832881828432 | 1.419135061376732 | 0.0029845979233988911 | 0.023969344303723984 | 0.12451729532439444 | -0.36683445070060317 | 0.71886582056152537 | 5.4701841780063292e-20 | 1.8578973401267354e-15 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 14.4375 | 13.375 | 0.073593073593073599 | 0.030303030303030304 | 2.4285714285714284 | -17 | 3.2769117377966068e-11 | 3.3370933113405186e-17 | 1.4553451008790879e-19 | PASS |
| 48 | r90 | 22.375 | 19.8125 | 0.11452513966480447 | 0 | null | -20.005951495444929 | 3.141619345032701e-12 | 4.6042541344059545e-15 | 1.5128771287691139e-18 | PASS |
| 48 | r99 | 26.1875 | 23.375 | 0.10739856801909307 | 0.0071599045346062056 | 15 | -17.172737481873973 | 2.8355402110316446e-11 | 5.1788953735629637e-16 | 2.4582963915343147e-19 | PASS |
| 48 | radial_profile_L2 | null | null | 0.22627052907292164 | 0.017773720507827696 | 12.730622661320131 | 12.215444655901132 | 3.3882278396906802e-09 | 3.3668681730326525e-21 | 2.4727213384738071e-11 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1875, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

### validation-statistical-active-r20

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 4 | total_mass | 961 | 971.1597207720456 | 0.010572029939693618 | 0.00052029136316337154 | 20.319441544091134 | 7.9957466414937608 | 8.6675086944026473e-07 | 1.6655970241068671e-31 | 2.0256904453958476e-27 | PASS |
| 4 | r_K_ratio | 0.95055031041461968 | 0.95537161884220823 | 0.0050721233529297041 | 0.0013209554508153041 | 3.8397383876944144 | 1.6923800514142024 | 0.1112352122045801 | 1.0681989599361004e-26 | 6.5880253584683872e-22 | PASS |
| 4 | active_fraction | 0.30956904900325755 | 0.29452200988976551 | 0.015047039113492049 | 0.0044779620119824587 | 3.3602426892474924 | -10.156332014753247 | 4.0728409442493956e-08 | 1.4231114756739095e-13 | 1.4616791576237938e-17 | PASS |
| 4 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r90 | 19.0625 | 18.625 | 0.022950819672131147 | 0.016393442622950821 | 1.3999999999999999 | -3.415650255319866 | 0.0038327022878305922 | 4.6454200247466158e-16 | 8.5990288790273657e-17 | PASS |
| 4 | r99 | 42.1875 | 39.4375 | 0.065185185185185179 | 0.014814814814814815 | 4.4000000000000004 | -7.2011903777877482 | 3.0674345362040858e-06 | 3.2295579237308931e-15 | 7.4750622844780019e-16 | PASS |
| 4 | radial_profile_L2 | null | null | 0.034915755052503467 | 0.010491987617756247 | 3.327849433734877 | 4.3630884530263003 | 0.00055660004916722506 | 3.820578568221231e-26 | 7.6802967655642231e-25 | PASS |
| 8 | total_mass | 930.6875 | 970.37748642098722 | 0.04264587890241054 | 0.0012087838291585521 | 35.279987929766406 | 18.224151302852587 | 1.2083610738640033e-11 | 2.4603784844513179e-28 | 5.7265887319147915e-23 | PASS |
| 8 | r_K_ratio | 0.89478658768318697 | 0.91508983096094598 | 0.02269059858209195 | 0.00095859770278864428 | 23.670616480806306 | 5.4168313475588086 | 7.1383887738472675e-05 | 4.92679901207683e-25 | 3.0796918725172109e-19 | PASS |
| 8 | active_fraction | 0.31580496549020443 | 0.27472069115986003 | 0.041084274330344436 | 0.0041874013475623299 | 9.8114011340855569 | -17.071898106328231 | 3.0849170249370308e-11 | 0.0010591161100719955 | 1.3293449287692704e-16 | PASS |
| 8 | r50 | 13.125 | 14 | 0.066666666666666666 | 0 | null | 10.2469507659596 | 3.6215788445874798e-08 | 4.2856040100792459e-20 | 2.7661345311125624e-13 | PASS |
| 8 | r90 | 22.9375 | 22.5 | 0.019073569482288829 | 0.073569482288828342 | 0.25925925925925924 | -0.706866391316882 | 0.49048560198715263 | 8.1642090029138098e-09 | 4.7890684845476589e-07 | PASS |
| 8 | r99 | 63.5 | 57.6875 | 0.091535433070866146 | 0.017716535433070866 | 5.166666666666667 | -9.7994319991289007 | 6.5189937337095449e-08 | 2.637761303481299e-15 | 2.0846796240234074e-15 | PASS |
| 8 | radial_profile_L2 | null | null | 0.055546805567441175 | 0.010917333453772486 | 5.0879462281375103 | 5.9890456730199961 | 2.4812894363272168e-05 | 5.9618048383494215e-25 | 7.2202456818107843e-23 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5546875, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 2.652994138115695e-12, "maximum_r99_to_half_width": 0.4765625, "outer_face_mass_present": true}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

### validation-statistical-active-r200

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: FAIL.

Duration is shorter than the sampling interval: endpoint statistical diagnostic only, without a resolved time-series claim.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 2 | total_mass | 982.375 | 989.7017335674069 | 0.0074581840614906776 | 0.00019086397760529329 | 39.075912359503491 | 4.9986902801061843 | 0.00015877195788384227 | 7.875421813558414e-31 | 1.2674639998173777e-26 | PASS |
| 2 | r_K_ratio | 0.97985804784240171 | 0.98185735548383302 | 0.0020404053891618511 | 0.00063597803570704424 | 3.2082953728007979 | 0.8014649647766281 | 0.43537548831832684 | 1.2799002618489483e-27 | 4.4723440163371157e-23 | PASS |
| 2 | active_fraction | 0.31117281779932982 | 0.29550364174646931 | 0.015669176052860528 | 0.0060028937985193084 | 2.61027040936915 | -18.35507386875776 | 1.0900552845839682e-11 | 5.396353190396315e-17 | 3.3688509997687006e-21 | PASS |
| 2 | r50 | 13 | 14 | 0.076923076923076927 | 0 | null | null | 0 | 0 | 0 | PASS |
| 2 | r90 | 33.25 | 26.6875 | 0.19736842105263158 | 0.080827067669172928 | 2.441860465116279 | -5.5801616220420946 | 5.2560085378087085e-05 | 0.053937023657414949 | 9.7952689369646178e-09 | FAIL |
| 2 | r99 | 71.5625 | 69.25 | 0.032314410480349345 | 0.0096069868995633193 | 3.3636363636363638 | -3.9698595833205057 | 0.0012326459082043812 | 2.1567240381686956e-17 | 1.0476149226677879e-15 | PASS |
| 2 | radial_profile_L2 | null | null | 0.081326884302303698 | 0.013394453118475617 | 6.071683821874414 | 11.035358301157789 | 1.3468597062346127e-08 | 2.4447976064259178e-26 | 2.9475055435703118e-23 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": true, "maximum_boundary_mass": 4, "maximum_r99_to_half_width": 1.203125, "outer_face_mass_present": true}, "baseline": {"annotated": true, "boundary_influenced": true, "maximum_boundary_mass": 7, "maximum_r99_to_half_width": 1.21875, "outer_face_mass_present": true}, "pde": {"annotated": true, "boundary_influenced": true, "maximum_boundary_mass": 7.5663361569535095, "maximum_r99_to_half_width": 1.09375, "outer_face_mass_present": true}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

### validation-statistical-active-r20-v14

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 4 | total_mass | 961 | 957.90656585838929 | 0.003218974132789498 | 0.00052029136316337154 | 6.1868682832214148 | -2.8447421948193705 | 0.012296571838415133 | 6.5482316044096946e-33 | 2.1587448948415302e-28 | PASS |
| 4 | r_K_ratio | 0.95055031041461968 | 0.94601696199436702 | 0.0047691830412166618 | 0.0013209554508153041 | 3.6104041497183603 | -1.6136684927102185 | 0.1274345568774421 | 1.0685895014976849e-26 | 3.7935550722800251e-22 | PASS |
| 4 | active_fraction | 0.30956904900325755 | 0.30629895168952725 | 0.0032700973137303226 | 0.0044779620119824587 | 0.73026463935601882 | -2.63156267037931 | 0.018873585616829489 | 1.4630200768729918e-16 | 2.0861532987049903e-17 | PASS |
| 4 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r90 | 19.0625 | 19 | 0.0032786885245901639 | 0.016393442622950821 | 0.20000000000000001 | -1 | 0.33317013591547739 | 6.3106316067013764e-23 | 8.9016836752560409e-20 | PASS |
| 4 | r99 | 42.1875 | 40.5625 | 0.038518518518518521 | 0.014814814814814815 | 2.6000000000000001 | -4.333333333333333 | 0.00059085778625744573 | 7.2263681135141884e-16 | 1.4960608813144921e-15 | PASS |
| 4 | radial_profile_L2 | null | null | 0.04402952038479812 | 0.010491987617756245 | 4.1964899301143106 | 7.4422381775128921 | 2.0725615554004122e-06 | 1.8320674186031458e-26 | 8.1220604059719425e-25 | PASS |
| 8 | total_mass | 930.6875 | 930.17624675850664 | 0.00054932857859741886 | 0.0012087838291585521 | 0.454447325771898 | -0.25820631899454072 | 0.79975686500682674 | 4.7640347720947908e-29 | 3.9086266215654425e-24 | PASS |
| 8 | r_K_ratio | 0.89478658768318697 | 0.89019784338966224 | 0.005128311439497709 | 0.00095859770278864428 | 5.3498056844690991 | -1.1556781954739936 | 0.26588985593582204 | 3.4656152627587374e-24 | 1.9236040346766121e-19 | PASS |
| 8 | active_fraction | 0.31580496549020443 | 0.30992378821835553 | 0.0058811772718488919 | 0.0041874013475623299 | 1.4044933321886099 | -3.1850839046096624 | 0.0061487594061861101 | 1.1819392671347893e-13 | 3.6532628616503073e-15 | PASS |
| 8 | r50 | 13.125 | 13.9375 | 0.061904761904761907 | 0 | null | 8.0622577482985509 | 7.8258249630702749e-07 | 3.6873899291617687e-18 | 8.5245610309508474e-13 | PASS |
| 8 | r90 | 22.9375 | 21.125 | 0.07901907356948229 | 0.073569482288828342 | 1.0740740740740742 | -3.092582041694214 | 0.0074289371221557829 | 6.4679713206351386e-08 | 2.9962515735592815e-08 | PASS |
| 8 | r99 | 63.5 | 60.9375 | 0.040354330708661415 | 0.017716535433070866 | 2.2777777777777777 | -4.2823103365101671 | 0.00065469959890511616 | 2.3014461693173097e-16 | 1.4852885343350244e-14 | PASS |
| 8 | radial_profile_L2 | null | null | 0.055529674363131512 | 0.010917333453772484 | 5.086377053358504 | 5.248900758347566 | 9.81457975014368e-05 | 4.6965365295240574e-25 | 5.6802490780922257e-23 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5546875, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 1.1642856634157104e-11, "maximum_r99_to_half_width": 0.5, "outer_face_mass_present": true}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

### validation-statistical-vascular-nonlinear

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: FAIL.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | vascular_length | 0 | 0 | null | null | null | null | 1 | null | null | FAIL |
| 0 | perfused_volume | 0 | 0 | null | null | null | null | 1 | null | null | FAIL |
| 0 | lesion_perfused_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | vascular_length | 6.5625 | 11.139162154530835 | 0.69739613783327015 | 0.038095238095238099 | 18.30664861812334 | 3.0731577191543407 | 0.0077294766562070235 | 4.2300259140499634e-14 | 0.44822464543172308 | FAIL |
| 4 | perfused_volume | 6.5625 | 11.139162154530835 | 0.69739613783327015 | 0.038095238095238099 | 18.30664861812334 | 3.0731577191543407 | 0.0077294766562070235 | 4.2300259140499634e-14 | 0.44822464543172308 | FAIL |
| 4 | lesion_perfused_fraction | 0.013263302548422297 | 0.030276155441830516 | 0.017012852893408221 | 0.00027370739056767172 | 62.157082635304086 | 5.6516307594264923 | 4.6022337053035405e-05 | 4.4587727375757038e-19 | 1.3327536049649137e-17 | PASS |
| 8 | vascular_length | 26.6875 | 44.809737788726686 | 0.6790534066033419 | 0.01405152224824356 | 48.325967436604493 | 5.6381327528177705 | 4.7189021802432674e-05 | 4.5665850820706246e-18 | 0.37057439858497371 | FAIL |
| 8 | perfused_volume | 26.6875 | 44.809737788726686 | 0.6790534066033419 | 0.01405152224824356 | 48.325967436604493 | 5.6381327528177705 | 4.7189021802432674e-05 | 4.5665850820706246e-18 | 0.37057439858497371 | FAIL |
| 8 | lesion_perfused_fraction | 0.037530736013472225 | 0.12027296915844937 | 0.082742233144977154 | 0.00049550457743789461 | 166.98580984420445 | 20.534579138960048 | 2.1522335257225706e-12 | 2.4437529466751145e-19 | 2.1282258366215037e-11 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.4583333333333333, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.4583333333333333, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.375, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

## Paired mechanism and correction interventions

### ablation-regular/ablation_report.json

Equivalence by intervention: `{"enlarged_boundary": true, "full": true, "no_activation": true, "no_division": true, "no_exchange": true, "no_migration": false}`.

| Intervention | Hours | Metric | ABM effect | PDE effect | Interaction/error change | SE | Paired t | Paired p |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| no_migration | 0 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_migration | 4 | total_mass | 0 | 0.00012560772118597185 | 0.00012560772118597185 | 0.00012133461811883484 | 1.0352175095070719 | 0.31696944719393905 |
| no_migration | 4 | r_K_ratio | 0 | 9.8130440741306391e-07 | 9.8130440741306391e-07 | 9.4792670324490787e-07 | 1.0352112711393178 | 0.31697226570925685 |
| no_migration | 4 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 4 | r50 | -0.0625 | -0.0625 | 0 | 0 | null | 1 |
| no_migration | 4 | r90 | -0.5625 | -0.3125 | 0.25 | 0.11180339887498948 | 2.2360679774997898 | 0.040968955955836141 |
| no_migration | 4 | r99 | -0.3125 | -0.4375 | -0.125 | 0.125 | -1 | 0.33317013591547739 |
| no_migration | 4 | radial_profile_L2 | null | null | -0.01390051111786134 | null | null | null |
| no_migration | 8 | total_mass | 0 | 0.00012560920900384076 | 0.00012560920900384076 | 0.00012133461754712231 | 1.0352297764902778 | 0.31696390498294125 |
| no_migration | 8 | r_K_ratio | 0 | 9.8130515333721968e-07 | 9.8130515333721968e-07 | 9.4792668872207049e-07 | 1.0352120739000901 | 0.3169719030182524 |
| no_migration | 8 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 8 | r50 | -0.0625 | -0.0625 | 0 | 0 | null | 1 |
| no_migration | 8 | r90 | -0.6875 | -0.625 | 0.0625 | 0.0625 | 1 | 0.33317013591547739 |
| no_migration | 8 | r99 | -1.125 | -1.0625 | 0.0625 | 0.0625 | 1 | 0.33317013591547739 |
| no_migration | 8 | radial_profile_L2 | null | null | -0.067252396223438593 | null | null | null |
| no_migration | 12 | total_mass | 0 | 0.00012561063768146141 | 0.00012561063768146141 | 0.00012133461747717505 | 1.0352415517780056 | 0.31695858498538909 |
| no_migration | 12 | r_K_ratio | 0 | 9.8130394693418665e-07 | 9.8130394693418665e-07 | 9.4792668289082791e-07 | 1.035210807592809 | 0.31697247514183602 |
| no_migration | 12 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 12 | r50 | -0.0625 | -0.0625 | 0 | 0 | null | 1 |
| no_migration | 12 | r90 | -1.3125 | -0.625 | 0.6875 | 0.11967838846954226 | 5.7445626465380286 | 3.8761814875216441e-05 |
| no_migration | 12 | r99 | -1.625 | -1.1875 | 0.4375 | 0.12808688457449499 | 3.415650255319866 | 0.0038327022878305922 |
| no_migration | 12 | radial_profile_L2 | null | null | -0.07289369138471738 | null | null | null |
| no_migration | 16 | total_mass | 0 | -0.010864838734402582 | -0.010864838734402582 | 0.0003033934222984094 | -35.811055665261676 | 6.0455933439191715e-16 |
| no_migration | 16 | r_K_ratio | 0 | -8.4881595078210859e-05 | -8.4881595078210859e-05 | 2.3702610893173368e-06 | -35.811073919564599 | 6.0455476219770827e-16 |
| no_migration | 16 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 16 | r50 | -0.125 | -0.0625 | 0.0625 | 0.0625 | 1 | 0.33317013591547739 |
| no_migration | 16 | r90 | -1.6875 | -0.625 | 1.0625 | 0.0625 | 17 | 3.2769117377966068e-11 |
| no_migration | 16 | r99 | -2.25 | -1.875 | 0.375 | 0.15478479684172258 | 2.4227185592617446 | 0.028528068449741654 |
| no_migration | 16 | radial_profile_L2 | null | null | -0.087752695635057346 | null | null | null |
| no_migration | 20 | total_mass | 0 | -5.8605528237900657 | -5.8605528237900657 | 0.11947625677421662 | -49.052029097841576 | 5.6032345700991519e-18 |
| no_migration | 20 | r_K_ratio | 0 | -0.045785548547017933 | -0.045785548547017933 | 0.00093340805550754197 | -49.0520177931419 | 5.6032538272925652e-18 |
| no_migration | 20 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 20 | r50 | -0.125 | -0.0625 | 0.0625 | 0.0625 | 1 | 0.33317013591547739 |
| no_migration | 20 | r90 | -2 | -0.75 | 1.25 | 0.11180339887498948 | 11.180339887498949 | 1.129731174656377e-08 |
| no_migration | 20 | r99 | -3 | -2.0625 | 0.9375 | 0.14343262065048754 | 6.5361700549589266 | 9.4202794798815216e-06 |
| no_migration | 20 | radial_profile_L2 | null | null | -0.091940361112018074 | null | null | null |
| no_migration | 24 | total_mass | -1.625 | -52.208642137442851 | -50.583642137442851 | 0.85454611061752672 | -59.193578332349126 | 3.3894646994228432e-19 |
| no_migration | 24 | r_K_ratio | -0.0126953125 | -0.40786485544709794 | -0.39516954294709794 | 0.0066760500641724679 | -59.192117966251551 | 3.3907142359259676e-19 |
| no_migration | 24 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 24 | r50 | -0.5 | 0 | 0.5 | 0.12909944487358058 | 3.8729833462074166 | 0.0015017735386984323 |
| no_migration | 24 | r90 | -2.25 | -0.875 | 1.375 | 0.15478479684172258 | 8.883301383959731 | 2.3171165556369958e-07 |
| no_migration | 24 | r99 | -3.8125 | -2.3125 | 1.5 | 0.25819888974716115 | 5.8094750193111251 | 3.4405086533123603e-05 |
| no_migration | 24 | radial_profile_L2 | null | null | -0.11489648391695023 | null | null | null |
| no_migration | 28 | total_mass | -5.8125 | -110.22584557842289 | -104.41334557842289 | 1.2846273607504619 | -81.279092107633474 | 2.9560412752486501e-21 |
| no_migration | 28 | r_K_ratio | -0.04541015625 | -0.85987853796974123 | -0.81446838171974123 | 0.010030108088671483 | -81.202353406305107 | 2.9981321885285906e-21 |
| no_migration | 28 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 28 | r50 | -0.75 | 0 | 0.75 | 0.11180339887498948 | 6.7082039324993694 | 7.0065583016648931e-06 |
| no_migration | 28 | r90 | -2.5625 | -1.375 | 1.1875 | 0.20854156260403664 | 5.6943085357748924 | 4.2526821576255181e-05 |
| no_migration | 28 | r99 | -3.6875 | -3 | 0.6875 | 0.19830006723817989 | 3.466968062972152 | 0.0034496088572533428 |
| no_migration | 28 | radial_profile_L2 | null | null | -0.12802134838502877 | null | null | null |
| no_migration | 32 | total_mass | -8.4375 | -133.73260416234154 | -125.29510416234154 | 1.4019127584014897 | -89.374394670041696 | 7.1350632067270752e-22 |
| no_migration | 32 | r_K_ratio | -0.03243103546162622 | -1.0219135233749639 | -0.98948248791333759 | 0.0083088809038005529 | -119.08733551118056 | 9.6854882537471302e-24 |
| no_migration | 32 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 32 | r50 | -0.8125 | 0 | 0.8125 | 0.10077822185373186 | 8.0622577482985509 | 7.8258249630702749e-07 |
| no_migration | 32 | r90 | -2.8125 | -1.5625 | 1.25 | 0.17078251276599329 | 7.3192505471139997 | 2.5290998830957293e-06 |
| no_migration | 32 | r99 | -4.3125 | -3.0625 | 1.25 | 0.19364916731037085 | 6.4549722436790278 | 1.0848127703079008e-05 |
| no_migration | 32 | radial_profile_L2 | null | null | -0.13647764960964145 | null | null | null |
| no_migration | 36 | total_mass | -10.3125 | -144.04916942694962 | -133.73666942694962 | 1.1943574759727429 | -111.97373660513823 | 2.4374630133662884e-23 |
| no_migration | 36 | r_K_ratio | -0.019031808192100444 | -0.97604494040665457 | -0.95701313221455409 | 0.008058910784688263 | -118.75216859738119 | 1.0103309231461033e-23 |
| no_migration | 36 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 36 | r50 | -0.875 | 0 | 0.875 | 0.085391256382996647 | 10.2469507659596 | 3.6215788445874798e-08 |
| no_migration | 36 | r90 | -3.3125 | -1.625 | 1.6875 | 0.17603858478564674 | 9.5859666337058052 | 8.6904680420919338e-08 |
| no_migration | 36 | r99 | -4.3125 | -3.3125 | 1 | 0.25819888974716115 | 3.8729833462074166 | 0.0015017735386984323 |
| no_migration | 36 | radial_profile_L2 | null | null | -0.14518622028853712 | null | null | null |
| no_migration | 40 | total_mass | -12.4375 | -161.24648992454425 | -148.80898992454425 | 0.97557404208921705 | -152.53479849245062 | 2.3704796524543965e-25 |
| no_migration | 40 | r_K_ratio | -0.014536219382929247 | -0.78248400390061623 | -0.76794778451768697 | 0.0087308617575123457 | -87.957844923717545 | 9.0633606016350317e-22 |
| no_migration | 40 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 40 | r50 | -0.875 | 0 | 0.875 | 0.085391256382996647 | 10.2469507659596 | 3.6215788445874798e-08 |
| no_migration | 40 | r90 | -3.4375 | -1.6875 | 1.75 | 0.19364916731037085 | 9.0369611411506394 | 1.8612315534939703e-07 |
| no_migration | 40 | r99 | -4.8125 | -3.75 | 1.0625 | 0.19297560294849017 | 5.5058773428660137 | 6.0385135056160928e-05 |
| no_migration | 40 | radial_profile_L2 | null | null | -0.15302764833551841 | null | null | null |
| no_migration | 44 | total_mass | -14.125 | -193.74725519981106 | -179.62225519981106 | 1.1204059711495673 | -160.31890209895408 | 1.1240520425816222e-25 |
| no_migration | 44 | r_K_ratio | -0.020820605221394506 | -0.56564354194606381 | -0.5448229367246693 | 0.0097409684950922198 | -55.931084983917849 | 7.9046520268331601e-19 |
| no_migration | 44 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 44 | r50 | -1 | 0 | 1 | 0.12909944487358058 | 7.7459666924148332 | 1.2783477850604916e-06 |
| no_migration | 44 | r90 | -3.8125 | -1.8125 | 2 | 0.091287092917527679 | 21.908902300206645 | 8.3886833819515926e-13 |
| no_migration | 44 | r99 | -5.1875 | -4 | 1.1875 | 0.20854156260403664 | 5.6943085357748924 | 4.2526821576255181e-05 |
| no_migration | 44 | radial_profile_L2 | null | null | -0.15794839798994056 | null | null | null |
| no_migration | 48 | total_mass | -25.4375 | -250.25356870911307 | -224.81606870911307 | 1.2413678535918904 | -181.10350454026107 | 1.8073466102290805e-26 |
| no_migration | 48 | r_K_ratio | -0.064987486008190057 | -0.49979575617272415 | -0.43480827016453411 | 0.011530124694820473 | -37.710630342086183 | 2.8068669962861586e-16 |
| no_migration | 48 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 48 | r50 | -1.3125 | -0.125 | 1.1875 | 0.13597640726733934 | 8.7331326357615122 | 2.8777724060903971e-07 |
| no_migration | 48 | r90 | -4.0625 | -2.125 | 1.9375 | 0.19297560294849017 | 10.040129272285084 | 4.7402960295905572e-08 |
| no_migration | 48 | r99 | -5.5 | -4.25 | 1.25 | 0.17078251276599329 | 7.3192505471139997 | 2.5290998830957293e-06 |
| no_migration | 48 | radial_profile_L2 | null | null | -0.17335345821972351 | null | null | null |
| no_activation | 0 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 0 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 0 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 0 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 0 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 0 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 0 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_activation | 4 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 4 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 4 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 4 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 4 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 4 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 4 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_activation | 8 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 8 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 8 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 8 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 8 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 8 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 8 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_activation | 12 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 12 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 12 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 12 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 12 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 12 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 12 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_activation | 16 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 16 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 16 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 16 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 16 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 16 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 16 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_activation | 20 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 20 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 20 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 20 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 20 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 20 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 20 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_activation | 24 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 24 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 24 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 24 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 24 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 24 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 24 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_activation | 28 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 28 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 28 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 28 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 28 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 28 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 28 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_activation | 32 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 32 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 32 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 32 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 32 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 32 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 32 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_activation | 36 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 36 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 36 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 36 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 36 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 36 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 36 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_activation | 40 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 40 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 40 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 40 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 40 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 40 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 40 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_activation | 44 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 44 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 44 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 44 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 44 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 44 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 44 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_activation | 48 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 48 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 48 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 48 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 48 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 48 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 48 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_division | 0 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 0 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 0 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 0 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 0 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 0 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 0 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_division | 4 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 4 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 4 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 4 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 4 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 4 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 4 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_division | 8 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 8 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 8 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 8 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 8 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 8 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 8 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_division | 12 | total_mass | 0 | -1.6704149174984195e-10 | -1.6704149174984195e-10 | 4.2428366466747232e-12 | -39.370238748352214 | 1.4801308555668398e-16 |
| no_division | 12 | r_K_ratio | 0 | -1.3049075708870816e-12 | -1.3049075708870816e-12 | 3.3175575872480892e-14 | -39.333381156753369 | 1.5008819501016545e-16 |
| no_division | 12 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 12 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 12 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 12 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 12 | radial_profile_L2 | null | null | 2.7935986857130501e-14 | null | null | null |
| no_division | 16 | total_mass | 0 | -0.0057853629485347113 | -0.0057853629485347113 | 0.00014103344494995876 | -41.021212738492373 | 8.0354819255330081e-17 |
| no_division | 16 | r_K_ratio | 0 | -4.5198147455322024e-05 | -4.5198147455322024e-05 | 1.1018237810373101e-06 | -41.021212496221771 | 8.0354826314745964e-17 |
| no_division | 16 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 16 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 16 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 16 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 16 | radial_profile_L2 | null | null | 9.9461584562865468e-07 | null | null | null |
| no_division | 20 | total_mass | 0 | -3.389430432232869 | -3.389430432232869 | 0.070228403626779887 | -48.262957111279007 | 7.1360593534784629e-18 |
| no_division | 20 | r_K_ratio | 0 | -0.026479904805026606 | -0.026479904805026606 | 0.00054865920192191241 | -48.262937561731341 | 7.1361024508245756e-18 |
| no_division | 20 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 20 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 20 | r90 | 0 | -0.0625 | -0.0625 | 0.0625 | -1 | 0.33317013591547739 |
| no_division | 20 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 20 | radial_profile_L2 | null | null | 0.0005885337022462217 | null | null | null |
| no_division | 24 | total_mass | -28.625 | -42.987672057546639 | -14.362672057546641 | 0.72798856400721934 | -19.72925505654538 | 3.8437424392118977e-12 |
| no_division | 24 | r_K_ratio | -0.2236328125 | -0.33582602662384331 | -0.11219321412384331 | 0.0056873976465268549 | -19.726634411145277 | 3.8511409710816958e-12 |
| no_division | 24 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 24 | r50 | -0.25 | 0 | 0.25 | 0.11180339887498948 | 2.2360679774997898 | 0.040968955955836141 |
| no_division | 24 | r90 | -0.0625 | -0.125 | -0.0625 | 0.11063265039459795 | -0.56493268286603204 | 0.58047096825964761 |
| no_division | 24 | r99 | -0.5625 | -0.1875 | 0.375 | 0.20155644370746373 | 1.8605210188381269 | 0.082530706346961857 |
| no_division | 24 | radial_profile_L2 | null | null | -0.014419662601844671 | null | null | null |
| no_division | 28 | total_mass | -104.0625 | -99.794414739433506 | 4.2680852605665009 | 0.85396222351109519 | 4.9979790007784199 | 0.00015899095787275111 |
| no_division | 28 | r_K_ratio | -0.81298828125 | -0.77838298445021348 | 0.034605296799786502 | 0.0066680036310712461 | 5.1897537425646236 | 0.00010988719720218108 |
| no_division | 28 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 28 | r50 | -0.25 | 0 | 0.25 | 0.11180339887498948 | 2.2360679774997898 | 0.040968955955836141 |
| no_division | 28 | r90 | -0.5625 | -0.375 | 0.1875 | 0.20854156260403664 | 0.89910134775393047 | 0.3828053154503297 |
| no_division | 28 | r99 | -1.0625 | -0.3125 | 0.75 | 0.19364916731037085 | 3.8729833462074166 | 0.0015017735386984323 |
| no_division | 28 | radial_profile_L2 | null | null | -0.028337468171891117 | null | null | null |
| no_division | 32 | total_mass | -129.5625 | -123.1432472957172 | 6.4192527042827905 | 1.1895367481275387 | 5.3964307654953894 | 7.4184439518388149e-05 |
| no_division | 32 | r_K_ratio | -0.82120050475896389 | -0.93918417274748689 | -0.11798366798852308 | 0.0080115482042361266 | -14.72670013096082 | 2.5159525249923727e-10 |
| no_division | 32 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 32 | r50 | -0.25 | 0 | 0.25 | 0.11180339887498948 | 2.2360679774997898 | 0.040968955955836141 |
| no_division | 32 | r90 | -0.6875 | -0.1875 | 0.5 | 0.12909944487358058 | 3.8729833462074166 | 0.0015017735386984323 |
| no_division | 32 | r99 | -1.1875 | 0 | 1.1875 | 0.13597640726733934 | 8.7331326357615122 | 2.8777724060903971e-07 |
| no_division | 32 | radial_profile_L2 | null | null | -0.028476163469525861 | null | null | null |
| no_division | 36 | total_mass | -148.5 | -133.43838746528112 | 15.061612534718885 | 1.2737737368967728 | 11.824401852885341 | 5.2861881415197011e-09 |
| no_division | 36 | r_K_ratio | -0.64525976610408464 | -0.8931482062059326 | -0.24788844010184791 | 0.013289294943228636 | -18.653242415103129 | 8.6421853869033513e-12 |
| no_division | 36 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 36 | r50 | 0 | 0 | 0 | 0.091287092917527679 | 0 | 1 |
| no_division | 36 | r90 | -1 | 0 | 1 | 0.12909944487358058 | 7.7459666924148332 | 1.2783477850604916e-06 |
| no_division | 36 | r99 | -0.8125 | -0.25 | 0.5625 | 0.24098322901535424 | 2.3341873303729379 | 0.033903591629393631 |
| no_division | 36 | radial_profile_L2 | null | null | -0.030440681039780981 | null | null | null |
| no_division | 40 | total_mass | -165.125 | -150.62569656877577 | 14.49930343122424 | 1.1960309053924181 | 12.122850141959345 | 3.7605139263964917e-09 |
| no_division | 40 | r_K_ratio | -0.4862818219349081 | -0.69950905566606869 | -0.21322723373116059 | 0.012911194671055348 | -16.514911219577453 | 4.9550568169313348e-11 |
| no_division | 40 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 40 | r50 | 0.0625 | 0 | -0.0625 | 0.0625 | -1 | 0.33317013591547739 |
| no_division | 40 | r90 | -0.8125 | -0.0625 | 0.75 | 0.11180339887498948 | 6.7082039324993694 | 7.0065583016648931e-06 |
| no_division | 40 | r99 | -1.1875 | -0.4375 | 0.75 | 0.26614532371118854 | 2.8180093098831724 | 0.012979205467040181 |
| no_division | 40 | radial_profile_L2 | null | null | -0.030944061815823959 | null | null | null |
| no_division | 44 | total_mass | -177 | -183.12427463431868 | -6.1242746343186827 | 1.420105691444101 | -4.3125484752412522 | 0.00061605976188320178 |
| no_division | 44 | r_K_ratio | -0.38511226905992024 | -0.48265150611580576 | -0.097539237055885492 | 0.0096093019281356211 | -10.150501855945937 | 4.1038407289767279e-08 |
| no_division | 44 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 44 | r50 | -0.0625 | 0 | 0.0625 | 0.11063265039459795 | 0.56493268286603204 | 0.58047096825964761 |
| no_division | 44 | r90 | -1.125 | -0.1875 | 0.9375 | 0.17001838135919303 | 5.5141096657035584 | 5.9461335275977702e-05 |
| no_division | 44 | r99 | -1.1875 | -0.25 | 0.9375 | 0.19297560294849017 | 4.8581270672347179 | 0.00020872892555786631 |
| no_division | 44 | radial_profile_L2 | null | null | -0.021768696833436713 | null | null | null |
| no_division | 48 | total_mass | -212.375 | -239.63037665047639 | -27.255376650476407 | 1.9562335049486281 | -13.932578386746396 | 5.4782020989102741e-10 |
| no_division | 48 | r_K_ratio | -0.4233832881828431 | -0.41680206803374842 | 0.0065812201490947006 | 0.011561299327221252 | 0.56924571908618604 | 0.57761476320715244 |
| no_division | 48 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 48 | r50 | -0.3125 | -0.125 | 0.1875 | 0.13597640726733934 | 1.3789156793307651 | 0.18813706218300907 |
| no_division | 48 | r90 | -1.3125 | -0.4375 | 0.875 | 0.20155644370746373 | 4.3412157106222962 | 0.00058157801081101909 |
| no_division | 48 | r99 | -1.5625 | -0.1875 | 1.375 | 0.25617376914898998 | 5.3674504012169324 | 7.8360069669564539e-05 |
| no_division | 48 | radial_profile_L2 | null | null | -0.03789874174959118 | null | null | null |
| no_exchange | 0 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 0 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 0 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 0 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 0 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 0 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 0 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_exchange | 4 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 4 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 4 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 4 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 4 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 4 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 4 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_exchange | 8 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 8 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 8 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 8 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 8 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 8 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 8 | radial_profile_L2 | null | null | -0.0010035764907716516 | null | null | null |
| no_exchange | 12 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 12 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 12 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 12 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 12 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 12 | r99 | 0.0625 | 0 | -0.0625 | 0.0625 | -1 | 0.33317013591547739 |
| no_exchange | 12 | radial_profile_L2 | null | null | -0.0009517645671109215 | null | null | null |
| no_exchange | 16 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 16 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 16 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 16 | r50 | -0.0625 | 0 | 0.0625 | 0.0625 | 1 | 0.33317013591547739 |
| no_exchange | 16 | r90 | 0.0625 | 0 | -0.0625 | 0.0625 | -1 | 0.33317013591547739 |
| no_exchange | 16 | r99 | -0.0625 | 0 | 0.0625 | 0.0625 | 1 | 0.33317013591547739 |
| no_exchange | 16 | radial_profile_L2 | null | null | -0.00082686899493643329 | null | null | null |
| no_exchange | 20 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 20 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 20 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 20 | r50 | 0.125 | 0 | -0.125 | 0.085391256382996647 | -1.4638501094227998 | 0.16387561365565406 |
| no_exchange | 20 | r90 | 0.0625 | 0 | -0.0625 | 0.11063265039459795 | -0.56493268286603204 | 0.58047096825964761 |
| no_exchange | 20 | r99 | 0 | 0 | 0 | 0.091287092917527679 | 0 | 1 |
| no_exchange | 20 | radial_profile_L2 | null | null | 0.0012586877602349944 | null | null | null |
| no_exchange | 24 | total_mass | 0.1875 | 0 | -0.1875 | 0.33189795118379384 | -0.56493268286603204 | 0.58047096825964761 |
| no_exchange | 24 | r_K_ratio | 0.00146484375 | 0 | -0.00146484375 | 0.0025929527436233894 | -0.56493268286603204 | 0.58047096825964761 |
| no_exchange | 24 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 24 | r50 | -0.125 | 0 | 0.125 | 0.125 | 1 | 0.33317013591547739 |
| no_exchange | 24 | r90 | -0.0625 | 0 | 0.0625 | 0.0625 | 1 | 0.33317013591547739 |
| no_exchange | 24 | r99 | 0 | 0 | 0 | 0.15811388300841897 | 0 | 1 |
| no_exchange | 24 | radial_profile_L2 | null | null | -0.0077316661724574354 | null | null | null |
| no_exchange | 28 | total_mass | 0.5625 | 0 | -0.5625 | 0.63224962106222993 | -0.8896802485305646 | 0.38768392744548597 |
| no_exchange | 28 | r_K_ratio | 0.00439453125 | 0 | -0.00439453125 | 0.0049394501645486713 | -0.8896802485305646 | 0.38768392744548597 |
| no_exchange | 28 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 28 | r50 | -0.125 | 0 | 0.125 | 0.085391256382996647 | 1.4638501094227998 | 0.16387561365565406 |
| no_exchange | 28 | r90 | -0.125 | 0 | 0.125 | 0.125 | 1 | 0.33317013591547739 |
| no_exchange | 28 | r99 | 0 | 0 | 0 | 0.12909944487358058 | 0 | 1 |
| no_exchange | 28 | radial_profile_L2 | null | null | -0.0026976319879465915 | null | null | null |
| no_exchange | 32 | total_mass | 0.875 | 0 | -0.875 | 0.83603727987054099 | -1.0466040463356996 | 0.31185515076135295 |
| no_exchange | 32 | r_K_ratio | 0.0065694230851235796 | 0 | -0.0065694230851235796 | 0.0058090117776161141 | -1.1309020082275545 | 0.27585445982534262 |
| no_exchange | 32 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 32 | r50 | -0.25 | 0 | 0.25 | 0.11180339887498948 | 2.2360679774997898 | 0.040968955955836141 |
| no_exchange | 32 | r90 | 0.0625 | 0 | -0.0625 | 0.0625 | -1 | 0.33317013591547739 |
| no_exchange | 32 | r99 | 0.125 | 0 | -0.125 | 0.125 | -1 | 0.33317013591547739 |
| no_exchange | 32 | radial_profile_L2 | null | null | -0.0014806734014796707 | null | null | null |
| no_exchange | 36 | total_mass | 0.25 | 0 | -0.25 | 0.68007352543677213 | -0.36760731104690392 | 0.71830125679313472 |
| no_exchange | 36 | r_K_ratio | 0.0069347190839245537 | 0 | -0.0069347190839245537 | 0.0054164797765932351 | -1.2803000047913471 | 0.21988483794049785 |
| no_exchange | 36 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 36 | r50 | -0.125 | 0 | 0.125 | 0.125 | 1 | 0.33317013591547739 |
| no_exchange | 36 | r90 | -0.125 | 0 | 0.125 | 0.125 | 1 | 0.33317013591547739 |
| no_exchange | 36 | r99 | 0.1875 | 0 | -0.1875 | 0.13597640726733934 | -1.3789156793307651 | 0.18813706218300907 |
| no_exchange | 36 | radial_profile_L2 | null | null | -0.0023135721269865739 | null | null | null |
| no_exchange | 40 | total_mass | 0.25 | 0 | -0.25 | 0.68007352543677213 | -0.36760731104690392 | 0.71830125679313472 |
| no_exchange | 40 | r_K_ratio | 0.0049084555455160134 | 0 | -0.0049084555455160134 | 0.0047817147511091465 | -1.0265053021779411 | 0.32092331968445276 |
| no_exchange | 40 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 40 | r50 | 0.0625 | 0 | -0.0625 | 0.0625 | -1 | 0.33317013591547739 |
| no_exchange | 40 | r90 | 0.25 | 0 | -0.25 | 0.11180339887498948 | -2.2360679774997898 | 0.040968955955836141 |
| no_exchange | 40 | r99 | 0.25 | 0 | -0.25 | 0.14433756729740643 | -1.7320508075688774 | 0.10377091535383504 |
| no_exchange | 40 | radial_profile_L2 | null | null | -0.003476029536625902 | null | null | null |
| no_exchange | 44 | total_mass | -0.0625 | 0 | 0.0625 | 0.74982636879035758 | 0.083352630157334795 | 0.93467335013319341 |
| no_exchange | 44 | r_K_ratio | 0.0006877176901675941 | 0 | -0.0006877176901675941 | 0.0047267990552590745 | -0.14549332055959388 | 0.88625848839631893 |
| no_exchange | 44 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 44 | r50 | -0.0625 | 0 | 0.0625 | 0.11063265039459795 | 0.56493268286603204 | 0.58047096825964761 |
| no_exchange | 44 | r90 | 0 | 0 | 0 | 0.091287092917527679 | 0 | 1 |
| no_exchange | 44 | r99 | 0.1875 | 0 | -0.1875 | 0.13597640726733934 | -1.3789156793307651 | 0.18813706218300907 |
| no_exchange | 44 | radial_profile_L2 | null | null | -0.00043652556331930104 | null | null | null |
| no_exchange | 48 | total_mass | 0.8125 | 0 | -0.8125 | 1.2049161450767711 | -0.6743207843299599 | 0.51036689320579987 |
| no_exchange | 48 | r_K_ratio | 0.0027878154931390148 | 0 | -0.0027878154931390148 | 0.006315749700019707 | -0.44140689950558298 | 0.66521490685887585 |
| no_exchange | 48 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 48 | r50 | -0.0625 | 0 | 0.0625 | 0.0625 | 1 | 0.33317013591547739 |
| no_exchange | 48 | r90 | 0.125 | 0 | -0.125 | 0.125 | -1 | 0.33317013591547739 |
| no_exchange | 48 | r99 | 0 | 0 | 0 | 0.15811388300841897 | 0 | 1 |
| no_exchange | 48 | radial_profile_L2 | null | null | 0.0010640748653448218 | null | null | null |
| enlarged_boundary | 0 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 0 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 0 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 0 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 0 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 0 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 0 | radial_profile_L2 | null | null | 0 | null | null | null |
| enlarged_boundary | 4 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 4 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 4 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 4 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 4 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 4 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 4 | radial_profile_L2 | null | null | 0 | null | null | null |
| enlarged_boundary | 8 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 8 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 8 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 8 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 8 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 8 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 8 | radial_profile_L2 | null | null | 0 | null | null | null |
| enlarged_boundary | 12 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 12 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 12 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 12 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 12 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 12 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 12 | radial_profile_L2 | null | null | 0 | null | null | null |
| enlarged_boundary | 16 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 16 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 16 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 16 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 16 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 16 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 16 | radial_profile_L2 | null | null | 0 | null | null | null |
| enlarged_boundary | 20 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 20 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 20 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 20 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 20 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 20 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 20 | radial_profile_L2 | null | null | 0 | null | null | null |
| enlarged_boundary | 24 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 24 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 24 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 24 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 24 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 24 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 24 | radial_profile_L2 | null | null | 0 | null | null | null |
| enlarged_boundary | 28 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 28 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 28 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 28 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 28 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 28 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 28 | radial_profile_L2 | null | null | 0 | null | null | null |
| enlarged_boundary | 32 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 32 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 32 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 32 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 32 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 32 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 32 | radial_profile_L2 | null | null | 0 | null | null | null |
| enlarged_boundary | 36 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 36 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 36 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 36 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 36 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 36 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 36 | radial_profile_L2 | null | null | 0 | null | null | null |
| enlarged_boundary | 40 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 40 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 40 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 40 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 40 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 40 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 40 | radial_profile_L2 | null | null | 0 | null | null | null |
| enlarged_boundary | 44 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 44 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 44 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 44 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 44 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 44 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 44 | radial_profile_L2 | null | null | 0 | null | null | null |
| enlarged_boundary | 48 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 48 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 48 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 48 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 48 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 48 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 48 | radial_profile_L2 | null | null | 0 | null | null | null |

#### full: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 255.99987439227883 | 4.9065516088270256e-07 | 0 | null | -1.0352175095070719 | 0.31696944719393905 | 6.3310518845588904e-81 | 6.3307856304481416e-81 | PASS |
| 4 | r_K_ratio | 1 | 0.9999990186955926 | 9.8130440741306391e-07 | 0 | null | -1.0352112711393178 | 0.31697226570925685 | 2.0746026777977294e-76 | 2.0744281865575098e-76 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13.0625 | 13.0625 | 0 | 0.0047846889952153108 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 4 | r90 | 17.9375 | 17.6875 | 0.013937282229965157 | 0.0069686411149825784 | 2 | -2.2360679774997898 | 0.040968955955836141 | 1.1148111900685798e-16 | 3.8956480726576652e-17 | PASS |
| 4 | r99 | 19.25 | 19.375 | 0.0064935064935064939 | 0.00974025974025974 | 0.66666666666666663 | 1 | 0.33317013591547739 | 1.2617854706008803e-18 | 7.1507823999893397e-17 | PASS |
| 4 | radial_profile_L2 | null | null | 0.01390051111786134 | 0.026445474852675486 | 0.52562909894034393 | -7.7751682619068863 | 1.2210780923911863e-06 | 1.7063425300474121e-19 | 5.5915970512964533e-19 | PASS |
| 8 | total_mass | 256 | 255.99987439079101 | 4.9066097267125297e-07 | 0 | null | -1.0352297764902778 | 0.31696390498294125 | 6.331051438669818e-81 | 6.3307851814239781e-81 | PASS |
| 8 | r_K_ratio | 1 | 0.99999901869484664 | 9.8130515333721968e-07 | 0 | null | -1.0352120739000901 | 0.3169719030182524 | 2.0746022011006527e-76 | 2.0744277097679295e-76 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 8 | r90 | 18.0625 | 18 | 0.0034602076124567475 | 0.0034602076124567475 | 1 | -1 | 0.33317013591547739 | 1.4301148022486204e-22 | 1.9701065876166452e-19 | PASS |
| 8 | r99 | 20.0625 | 20 | 0.0031152647975077881 | 0.0031152647975077881 | 1 | -1 | 0.33317013591547739 | 6.5149803197211003e-25 | 5.0766352092358305e-21 | PASS |
| 8 | radial_profile_L2 | null | null | 0.067252396223438593 | 0.031674679515173401 | 2.1232226261743889 | -1.1211565613575198 | 0.2798500299122274 | 1.1233733444798528e-17 | 3.6201109868020346e-15 | PASS |
| 12 | total_mass | 256 | 255.99987438921451 | 4.9066713076612034e-07 | 0 | null | -1.035242772400685 | 0.31695803351978358 | 6.3310511518897591e-81 | 6.3307848913142876e-81 | PASS |
| 12 | r_K_ratio | 1 | 0.99999901869489838 | 9.8130510160082673e-07 | 0 | null | -1.0352120282386914 | 0.31697192364827376 | 2.0746019330439565e-76 | 2.0744274417429762e-76 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 12 | r90 | 18.6875 | 18 | 0.03678929765886288 | 0 | null | -5.7445626465380286 | 3.8761814875216441e-05 | 1.2411344634416193e-17 | 3.0051609001276733e-16 | PASS |
| 12 | r99 | 20.5625 | 20.125 | 0.021276595744680851 | 0.015197568389057751 | 1.3999999999999999 | -3.415650255319866 | 0.0038327022878305922 | 3.2453744019125243e-19 | 3.5829058514479523e-17 | PASS |
| 12 | radial_profile_L2 | null | null | 0.072893691384749507 | 0.033001760468462815 | 2.2087819058746359 | -0.25736842332843862 | 0.80039152143350745 | 4.3685958447085025e-20 | 2.3969087135186269e-17 | PASS |
| 16 | total_mass | 256 | 256.00565975007765 | 2.2108398740866564e-05 | 0 | null | 29.0549503061726 | 1.3336686824631748e-14 | 7.6725742637134177e-78 | 7.6871276199674024e-78 | PASS |
| 16 | r_K_ratio | 1 | 1.0000442168399475 | 4.4216839947416875e-05 | 0 | null | 29.054978697547448 | 1.3336494557056482e-14 | 2.5117676850929536e-73 | 2.5213053845686074e-73 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.0047619047619047623 | 1 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 16 | r90 | 19.0625 | 18 | 0.055737704918032788 | 0.0032786885245901639 | 17 | -17 | 3.2769117377966068e-11 | 2.2623006381054966e-21 | 5.3331248280097439e-21 | PASS |
| 16 | r99 | 21.1875 | 20.8125 | 0.017699115044247787 | 0.014749262536873156 | 1.2 | -2.4227185592617446 | 0.028528068449741654 | 4.6950640203719188e-18 | 3.0093626661365215e-16 | PASS |
| 16 | radial_profile_L2 | null | null | 0.08775376074242433 | 0.029422120612765671 | 2.9825776971476938 | 0.20757506469588111 | 0.83835262262231791 | 5.9173726605820711e-20 | 1.2294617205940524e-16 | PASS |
| 20 | total_mass | 256 | 259.38930481728698 | 0.013239471942527239 | 0 | null | 48.255029981106148 | 7.1535576417554098e-18 | 9.9612836898674566e-40 | 3.1001936901933769e-39 | PASS |
| 20 | r_K_ratio | 1 | 1.0264789234962126 | 0.026478923496212517 | 0 | null | 48.255010551841188 | 7.1536005859322785e-18 | 1.9075543162085774e-35 | 1.8535898820419411e-34 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0 | null | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 20 | r90 | 19.375 | 18.125 | 0.064516129032258063 | 0 | null | -11.180339887498949 | 1.129731174656377e-08 | 6.9886341424657286e-17 | 1.1645945213682297e-17 | PASS |
| 20 | r99 | 21.9375 | 21 | 0.042735042735042736 | 0 | null | -6.5361700549589266 | 9.4202794798815216e-06 | 3.6683425778815099e-19 | 5.1861867803112176e-17 | PASS |
| 20 | radial_profile_L2 | null | null | 0.092421060738694855 | 0.040888641583758337 | 2.2603113519771729 | 1.0369297345499009 | 0.31619654317416779 | 1.1226962030547471e-19 | 3.5490613354533042e-16 | PASS |
| 24 | total_mass | 284.625 | 298.98754644036325 | 0.050461296233160362 | 0.0083443126921387799 | 6.047387974889797 | 19.730067038543012 | 3.8414531453545918e-12 | 1.6655638915222936e-27 | 4.0384969275713948e-22 | PASS |
| 24 | r_K_ratio | 1.2236328125 | 1.3358250453136638 | 0.091687826337742723 | 0.015163607342378291 | 6.0465708632206123 | 19.727446316178806 | 3.8488472047269843e-12 | 2.9675877127982947e-24 | 2.7669535033326392e-17 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.5625 | 13.0625 | 0.036866359447004608 | 0.018433179723502304 | 2 | -3.8729833462074166 | 0.0015017735386984323 | 9.6072681659644597e-15 | 7.6467248765302002e-14 | PASS |
| 24 | r90 | 19.9375 | 18.25 | 0.084639498432601878 | 0.003134796238244514 | 27 | -14.100290132411525 | 4.6333644836412128e-10 | 7.7245411315955695e-15 | 9.481589102890017e-19 | PASS |
| 24 | r99 | 22.875 | 21.25 | 0.071038251366120214 | 0.0054644808743169399 | 13 | -6.7890285822722154 | 6.1053729318667233e-06 | 1.2337518345797697e-14 | 7.416708381213179e-15 | PASS |
| 24 | radial_profile_L2 | null | null | 0.1311571822925125 | 0.031998880089473554 | 4.0988053933693243 | 4.1375477785895995 | 0.00087694673719733397 | 3.3174571066196259e-19 | 3.9701746118309766e-14 | PASS |
| 28 | total_mass | 360.0625 | 355.79428911997218 | 0.011854083332831953 | 0.0034716195105016492 | 3.4145687040222441 | -4.9981373855468707 | 0.00015894216488354686 | 4.9346299105874915e-27 | 1.1189772562667494e-23 | PASS |
| 28 | r_K_ratio | 1.81298828125 | 1.7783820031380952 | 0.019087976723183724 | 0.0053864799353622404 | 3.543682878659058 | -5.1899130998290026 | 0.00010985368152231689 | 4.8401929701778768e-24 | 5.9660093150193663e-21 | PASS |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.875 | 13.0625 | 0.058558558558558557 | 0.022522522522522521 | 2.6000000000000001 | -8.0622577482985509 | 7.8258249630702749e-07 | 2.2778649613890186e-15 | 2.5893641111623486e-16 | PASS |
| 28 | r90 | 20.5625 | 18.75 | 0.088145896656534953 | 0.0060790273556231003 | 14.5 | -11.066875266096559 | 1.2961577489987425e-08 | 1.5680800107416445e-13 | 1.9477056267575937e-16 | PASS |
| 28 | r99 | 23.625 | 21.9375 | 0.071428571428571425 | 0 | null | -8.5098307000225546 | 3.9910679004152207e-07 | 1.2402986697836266e-16 | 5.9545457575235637e-16 | PASS |
| 28 | radial_profile_L2 | null | null | 0.16812366489946665 | 0.025528598153441336 | 6.5856990614583761 | 8.2179631045461328 | 6.1744432441080883e-07 | 2.0668183461795087e-20 | 1.1633368762049021e-13 | PASS |
| 32 | total_mass | 385.5625 | 379.14312167383497 | 0.016649384538602752 | 0.0024315124007132437 | 6.8473368812426925 | -5.396509588004732 | 7.4173402428140978e-05 | 1.4730519212317234e-25 | 2.2336216332001333e-22 | PASS |
| 32 | r_K_ratio | 1.8212005047589639 | 1.9391831914335063 | 0.064782920038866071 | 0.0014984898683687194 | 43.232137504799965 | 14.726328885627707 | 2.5168461149834619e-10 | 2.6740318924724625e-25 | 4.4405753482401013e-18 | PASS |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.022421524663677129 | 2.7999999999999998 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 32 | r90 | 20.8125 | 18.9375 | 0.090090090090090086 | 0.006006006006006006 | 15 | -21.957751641341996 | 8.1207854767429677e-13 | 2.5303126792687863e-18 | 2.9422487852212645e-20 | PASS |
| 32 | r99 | 24.375 | 22 | 0.097435897435897437 | 0.0051282051282051282 | 19 | -15.343884208657718 | 1.4090734244550085e-10 | 8.3953741872725587e-18 | 3.7020237911275972e-18 | PASS |
| 32 | radial_profile_L2 | null | null | 0.17938424840949052 | 0.026253620259367155 | 6.8327433183424358 | 10.818675251643851 | 1.7576228289281234e-08 | 2.8179852830823343e-21 | 5.7756288571017507e-14 | PASS |
| 36 | total_mass | 404.5 | 389.43826184084367 | 0.037235446623377903 | 0.001854140914709518 | 20.082317545541816 | -11.824525086757328 | 5.2854372155451087e-09 | 1.4574793360891194e-25 | 2.4656946513950525e-22 | PASS |
| 36 | r_K_ratio | 1.6452597661040846 | 1.8931472248902017 | 0.15066767199510719 | 0.0025243036803777237 | 59.686824991104899 | 18.653271937262037 | 8.6419882044172036e-12 | 9.373935059396823e-23 | 6.9902926261365534e-12 | PASS |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0089686098654708519 | 7 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 36 | r90 | 21.3125 | 19 | 0.10850439882697947 | 0 | null | -19.322619811082461 | 5.195040620846464e-12 | 7.7864380968352111e-16 | 1.5335258890741899e-18 | PASS |
| 36 | r99 | 24.625 | 22.25 | 0.096446700507614211 | 0.0025380710659898475 | 38 | -11.783299785974805 | 5.5430657455111653e-09 | 3.3901457391067851e-15 | 5.0500172350014425e-17 | PASS |
| 36 | radial_profile_L2 | null | null | 0.19113782511778182 | 0.021200415960106734 | 9.0157582510385588 | 12.937967455874144 | 1.5352977834313781e-09 | 2.3205167683510447e-21 | 1.873559415950491e-13 | PASS |
| 40 | total_mass | 421.125 | 406.62557094179493 | 0.034430226318088571 | 0.0014841199168892847 | 23.199086493128078 | -12.123072892031722 | 3.7595679941645351e-09 | 1.454436824969024e-26 | 8.6049179868181217e-23 | PASS |
| 40 | r_K_ratio | 1.4862818219349081 | 1.6995080743478872 | 0.1434628677187155 | 0.0052088165665525443 | 27.542315204558335 | 16.515053506398107 | 4.9544479197942529e-11 | 1.9779308949787462e-22 | 1.6280039013798181e-11 | PASS |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0044843049327354259 | 14 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 40 | r90 | 21.5 | 19.0625 | 0.11337209302325581 | 0 | null | -19.030051422496395 | 6.4760803247711866e-12 | 6.3605411260372649e-15 | 2.3244786420838592e-18 | PASS |
| 40 | r99 | 25 | 22.6875 | 0.092499999999999999 | 0.0074999999999999997 | 12.333333333333334 | -15.363413772941893 | 1.3839300528290396e-10 | 5.4380328870719136e-17 | 4.8751207713627678e-19 | PASS |
| 40 | radial_profile_L2 | null | null | 0.20294210478022673 | 0.019135549966259355 | 10.605501547541785 | 12.331063058139474 | 2.9774061335968384e-09 | 2.6639476681283848e-22 | 9.5858061539418275e-14 | PASS |
| 44 | total_mass | 433 | 439.12414900465899 | 0.014143531188588849 | 0.0072170900692840644 | 1.9597276814908708 | 4.3124362386774742 | 0.00061619879125071012 | 5.6783101501054603e-26 | 3.3785584158268306e-21 | PASS |
| 44 | r_K_ratio | 1.3851122690599202 | 1.4826505247964381 | 0.070419025168780944 | 0.023883569861476831 | 2.948429634983663 | 10.150700193299443 | 4.1027820504537489e-08 | 1.7459206768544517e-22 | 6.5239376359162711e-15 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 14.0625 | 13.0625 | 0.071111111111111111 | 0.0088888888888888889 | 8 | -7.7459666924148332 | 1.2783477850604916e-06 | 1.25956553381152e-13 | 5.7635238603637946e-15 | PASS |
| 44 | r90 | 21.875 | 19.1875 | 0.12285714285714286 | 0.0057142857142857143 | 21.5 | -22.456017618285021 | 5.8544703070058476e-13 | 4.6040517988789706e-14 | 9.7976516786557512e-20 | PASS |
| 44 | r99 | 25.4375 | 22.9375 | 0.098280098280098274 | 0.0024570024570024569 | 40 | -15.811388300841896 | 9.2073297610083179e-11 | 1.32240574860815e-17 | 2.0753132400270518e-18 | PASS |
| 44 | radial_profile_L2 | null | null | 0.20578903863170117 | 0.022666905945975712 | 9.0788323347782285 | 13.563001661788283 | 7.9734125142345473e-10 | 3.9888272112335183e-23 | 2.1411534627969286e-14 | PASS |
| 48 | total_mass | 468.375 | 495.63025101823672 | 0.058191088376272695 | 0.0066720042700827327 | 8.7216803258357523 | 13.932439671958022 | 5.478965023639955e-10 | 6.2262438154028136e-25 | 9.8983062404198524e-19 | PASS |
| 48 | r_K_ratio | 1.4233832881828432 | 1.4168010867127527 | 0.00462433521928831 | 0.023969344303723984 | 0.19292706386506586 | -0.56934483893448984 | 0.5775492092836918 | 5.7033749883260893e-20 | 1.7002683649562875e-15 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 14.4375 | 13.1875 | 0.086580086580086577 | 0.030303030303030304 | 2.8571428571428572 | -11.180339887498949 | 1.129731174656377e-08 | 5.7520127976418164e-14 | 2.7873182537487998e-16 | PASS |
| 48 | r90 | 22.375 | 19.5 | 0.12849162011173185 | 0 | null | -15.998991903725805 | 7.7865353236989104e-11 | 1.0061554522861548e-11 | 4.8264550809103502e-17 | PASS |
| 48 | r99 | 26.1875 | 23.1875 | 0.11455847255369929 | 0.0071599045346062056 | 16 | -23.2379000772445 | 3.5513304925948273e-13 | 1.7995348202025614e-17 | 8.6169881800138773e-21 | PASS |
| 48 | radial_profile_L2 | null | null | 0.23341720836842744 | 0.017773720507827693 | 13.1327151378142 | 12.748408916387797 | 1.8827621113569947e-09 | 3.1119166670442979e-21 | 6.3693627289909649e-11 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1875, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

#### no_migration: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: FAIL.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 4 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 4 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 8 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 8 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 8 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 12 | total_mass | 256 | 255.99999999985221 | 5.7732291169898531e-13 | 0 | null | -39.46596347006362 | 1.427652693410483e-16 | 1.3907099505504748e-193 | 1.3907099504816198e-193 | PASS |
| 12 | r_K_ratio | 1 | 0.99999999999884537 | 1.1546666400796823e-12 | 0 | null | -39.470000494738905 | 1.4254835719827034e-16 | 4.5513224768269243e-189 | 4.5513224763765052e-189 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 12 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 12 | radial_profile_L2 | null | null | 3.2132089803357784e-14 | 0.036932043236795692 | 8.7003282210349854e-13 | -13.759882903644412 | 6.5213491598757298e-10 | 5.2388483704395915e-194 | 5.2388483704541834e-194 | PASS |
| 16 | total_mass | 256 | 255.99479491134326 | 2.0332377565393522e-05 | 0 | null | -41.039370124078054 | 7.9827591650811398e-17 | 1.2316461238656817e-80 | 1.2295015108509992e-80 | PASS |
| 16 | r_K_ratio | 1 | 0.99995933524486924 | 4.0664755130793984e-05 | 0 | null | -41.03937012388127 | 7.9827591656505162e-17 | 4.039376652521207e-76 | 4.0253217102264507e-76 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 16 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 16 | radial_profile_L2 | null | null | 1.0651073669858551e-06 | 0.036932043236795692 | 2.8839654501559774e-05 | -13.759494193661039 | 6.5239220566091882e-10 | 1.0285338594359742e-80 | 1.0286277636359151e-80 | PASS |
| 20 | total_mass | 256 | 253.52875199349691 | 0.0096533125254027047 | 0 | null | -50.135799233660236 | 4.0447037778220632e-18 | 1.3043816924938161e-41 | 5.7012250344599849e-42 | PASS |
| 20 | r_K_ratio | 1 | 0.98069337494919462 | 0.019306625050805416 | 0 | null | -50.135799233659661 | 4.0447037778227542e-18 | 6.5803273139127886e-37 | 1.2555551427673749e-37 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 20 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 20 | radial_profile_L2 | null | null | 0.00048069962667678665 | 0.036932043236795692 | 0.013015787499075105 | -13.595089949751127 | 7.7150843682299059e-10 | 3.3737108636554468e-41 | 3.5156200969121823e-41 | PASS |
| 24 | total_mass | 283 | 246.77890430292041 | 0.1279897374455109 | 0.0079505300353356883 | 16.098264754257592 | -40.484230232592964 | 9.7752489848969245e-17 | 3.2818362105478331e-23 | 9.7815777769852664e-24 | PASS |
| 24 | r_K_ratio | 1.2109375 | 0.9279601898665657 | 0.23368448836825537 | 0.014516129032258065 | 16.098264754257592 | -40.484230232592964 | 9.7752489848969245e-17 | 4.0061669251704097e-15 | 4.040582789909796e-21 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 24 | r90 | 17.6875 | 17.375 | 0.017667844522968199 | 0.0070671378091872791 | 2.5 | -2.6111648393354674 | 0.019657034541711152 | 1.4826735422532143e-16 | 4.5581962869429108e-16 | PASS |
| 24 | r99 | 19.0625 | 18.9375 | 0.0065573770491803279 | 0.0065573770491803279 | 1 | -1.4638501094227998 | 0.16387561365565406 | 4.2575407832340205e-21 | 1.8076952773475326e-19 | PASS |
| 24 | radial_profile_L2 | null | null | 0.016260698375562262 | 0.031004829841119182 | 0.52445694618833294 | -13.730419384384332 | 6.7194590618668282e-10 | 4.2748402464487196e-19 | 1.7125541349342287e-18 | PASS |
| 28 | total_mass | 354.25 | 245.5684435415493 | 0.30679338449809657 | 0.0061750176429075515 | 49.682997238148893 | -81.007231745106608 | 3.108061399370375e-21 | 1.1052934409010355e-11 | 1.5069501925999661e-24 | PASS |
| 28 | r_K_ratio | 1.767578125 | 0.91850346516835391 | 0.48036047053458875 | 0.0096685082872928173 | 49.682997238148893 | -81.007231745106537 | 3.108061399370419e-21 | 0.99999999999999922 | 3.7142410684516703e-23 | FAIL |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.014285714285714285 | 0.33333333333333331 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 28 | r90 | 18 | 17.375 | 0.034722222222222224 | 0 | null | -5 | 0.0001583695146220272 | 2.5609883229122797e-15 | 4.0487734815626444e-17 | PASS |
| 28 | r99 | 19.9375 | 18.9375 | 0.050156739811912224 | 0.006269592476489028 | 8 | null | 0 | 2.8880752909119404e-29 | 1.8281234277735476e-31 | PASS |
| 28 | radial_profile_L2 | null | null | 0.040102316514437887 | 0.035579270731715518 | 1.12712587103958 | -3.3501685779973447 | 0.0043838252798495438 | 1.1544744828317028e-18 | 3.5639929123899286e-17 | PASS |
| 32 | total_mass | 377.125 | 245.41051751149342 | 0.34925948289958653 | 0.0041431885979449782 | 84.297268792644203 | -122.67218569008094 | 6.210122896748211e-24 | 0.34775891277754345 | 8.2675696869125977e-27 | FAIL |
| 32 | r_K_ratio | 1.7887694692973377 | 0.91726966805854238 | 0.48720632602318326 | 0.0027342015763558028 | 178.18961492683408 | -83.286861360090114 | 2.0515236723121589e-21 | 0.99999999999999967 | 2.7215099912080891e-23 | FAIL |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.0095238095238095247 | 0.5 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 32 | r90 | 18 | 17.375 | 0.034722222222222224 | 0 | null | -5 | 0.0001583695146220272 | 2.5609883229122797e-15 | 4.0487734815626444e-17 | PASS |
| 32 | r99 | 20.0625 | 18.9375 | 0.056074766355140186 | 0 | null | -13.174650984805199 | 1.1942408086569779e-09 | 3.1450495481497661e-20 | 8.9736644299638471e-21 | PASS |
| 32 | radial_profile_L2 | null | null | 0.042906598799849066 | 0.026598751298356628 | 1.6131057551750558 | -2.2952229889878817 | 0.036560204498734393 | 1.2833061626887599e-19 | 5.0769430755606313e-18 | PASS |
| 36 | total_mass | 394.1875 | 245.38909241389402 | 0.37748129401897823 | 0.0019026478515934676 | 198.39787678147465 | -132.13782890416067 | 2.0384777019909154e-24 | 0.99999999987985011 | 4.846584309414167e-27 | FAIL |
| 36 | r_K_ratio | 1.6262279579119843 | 0.91710228448354703 | 0.43605551729594411 | 0.0062477019091505366 | 69.794545840493214 | -42.145769418517531 | 5.3734857948289266e-17 | 0.99999999910111659 | 3.5945722125298909e-19 | FAIL |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.0625 | 13.0625 | 0 | 0 | null | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 36 | r90 | 18 | 17.375 | 0.034722222222222224 | 0 | null | -5 | 0.0001583695146220272 | 2.5609883229122797e-15 | 4.0487734815626444e-17 | PASS |
| 36 | r99 | 20.3125 | 18.9375 | 0.067692307692307691 | 0 | null | -11 | 1.4062516106729139e-08 | 2.4748104659445616e-18 | 4.5290686638199662e-18 | PASS |
| 36 | radial_profile_L2 | null | null | 0.045951604829244722 | 0.022848576447345625 | 2.0111364458586669 | -2.8877373493029102 | 0.011271435637783987 | 7.1976174195233755e-19 | 3.6869875041834378e-17 | PASS |
| 40 | total_mass | 408.6875 | 245.37908101725068 | 0.39959240001896146 | 0.0013763572411683743 | 290.32607819155436 | -135.26068491136041 | 1.4363485464011021e-24 | 0.99999999999996203 | 5.1680317959373897e-27 | FAIL |
| 40 | r_K_ratio | 1.4717456025519788 | 0.91702407044727097 | 0.3769140068384314 | 0.0092852152085018697 | 40.592920936642983 | -30.769246779698943 | 5.720061734805753e-15 | 0.99793893807091993 | 1.436298348042337e-17 | FAIL |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 13.0625 | 13.0625 | 0 | 0 | null | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 40 | r90 | 18.0625 | 17.375 | 0.038062283737024222 | 0.0034602076124567475 | 11 | -5.7445626465380286 | 3.8761814875216441e-05 | 1.3017439639714619e-15 | 2.6565258179719937e-17 | PASS |
| 40 | r99 | 20.1875 | 18.9375 | 0.061919504643962849 | 0.0030959752321981426 | 20 | -11.180339887498949 | 1.129731174656377e-08 | 5.8869867743851206e-19 | 8.9512236720867488e-19 | PASS |
| 40 | radial_profile_L2 | null | null | 0.049914456444708329 | 0.024227567023808398 | 2.0602339638832681 | -1.6439962716571783 | 0.12096705294937792 | 7.4091061495801653e-19 | 5.3490445234045498e-17 | PASS |
| 44 | total_mass | 418.875 | 245.3768938048479 | 0.4142001938410077 | 0.0032826022082960309 | 126.18044086920153 | -115.47774215397257 | 1.536120499926264e-23 | 0.99999999999998535 | 7.1127855299788061e-26 | FAIL |
| 44 | r_K_ratio | 1.3642916638385258 | 0.91700698285037419 | 0.32785121601467987 | 0.02354061237077195 | 13.927047047499084 | -28.5953319470241 | 1.6872866009022303e-14 | 0.0047490299904603987 | 1.5175753873776502e-17 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 13.0625 | 13.0625 | 0 | 0 | null | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 44 | r90 | 18.0625 | 17.375 | 0.038062283737024222 | 0.0034602076124567475 | 11 | -5.7445626465380286 | 3.8761814875216441e-05 | 1.3017439639714619e-15 | 2.6565258179719937e-17 | PASS |
| 44 | r99 | 20.25 | 18.9375 | 0.064814814814814811 | 0.0030864197530864196 | 21 | -10.966892325208963 | 1.464381425614242e-08 | 1.3801419874098601e-18 | 2.4873321813140284e-18 | PASS |
| 44 | radial_profile_L2 | null | null | 0.047840640641760614 | 0.021099978364348116 | 2.2673312652583211 | -1.6676185722967594 | 0.11612740088911967 | 1.9332203664977814e-19 | 1.1712153677336303e-17 | PASS |
| 48 | total_mass | 442.9375 | 245.37668230912368 | 0.44602414040553423 | 0.0046564131508395655 | 95.787063122849119 | -100.85352470431691 | 1.1680519398189367e-22 | 0.99999999999999911 | 8.9094608732779362e-25 | FAIL |
| 48 | r_K_ratio | 1.3583958021746529 | 0.91700533054002864 | 0.32493509691947176 | 0.023109012801366548 | 14.060968320561784 | -25.478401460891579 | 9.219913136519297e-14 | 0.0042755581864384884 | 7.8882396414004588e-17 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.0047619047619047623 | 1 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 48 | r90 | 18.3125 | 17.375 | 0.051194539249146756 | 0.010238907849829351 | 5 | -6.5361700549589266 | 9.4202794798815216e-06 | 1.0283171831955228e-14 | 6.4688436157782434e-16 | PASS |
| 48 | r99 | 20.6875 | 18.9375 | 0.084592145015105744 | 0.0030211480362537764 | 28 | -15.652475842498529 | 1.0626678392473952e-10 | 7.4783111220355879e-19 | 4.6220321327614681e-19 | PASS |
| 48 | radial_profile_L2 | null | null | 0.060063750148703944 | 0.023214227865607628 | 2.5873679924409489 | 0.51102559010583459 | 0.61677195402910201 | 1.85907590651432e-20 | 3.2946471860339307e-18 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1640625, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1640625, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1484375, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

#### no_activation: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 255.99987439227883 | 4.9065516088270256e-07 | 0 | null | -1.0352175095070719 | 0.31696944719393905 | 6.3310518845588904e-81 | 6.3307856304481416e-81 | PASS |
| 4 | r_K_ratio | 1 | 0.9999990186955926 | 9.8130440741306391e-07 | 0 | null | -1.0352112711393178 | 0.31697226570925685 | 2.0746026777977294e-76 | 2.0744281865575098e-76 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13.0625 | 13.0625 | 0 | 0.0047846889952153108 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 4 | r90 | 17.9375 | 17.6875 | 0.013937282229965157 | 0.0069686411149825784 | 2 | -2.2360679774997898 | 0.040968955955836141 | 1.1148111900685798e-16 | 3.8956480726576652e-17 | PASS |
| 4 | r99 | 19.25 | 19.375 | 0.0064935064935064939 | 0.00974025974025974 | 0.66666666666666663 | 1 | 0.33317013591547739 | 1.2617854706008803e-18 | 7.1507823999893397e-17 | PASS |
| 4 | radial_profile_L2 | null | null | 0.01390051111786134 | 0.026445474852675486 | 0.52562909894034393 | -7.7751682619068863 | 1.2210780923911863e-06 | 1.7063425300474121e-19 | 5.5915970512964533e-19 | PASS |
| 8 | total_mass | 256 | 255.99987439079101 | 4.9066097267125297e-07 | 0 | null | -1.0352297764902778 | 0.31696390498294125 | 6.331051438669818e-81 | 6.3307851814239781e-81 | PASS |
| 8 | r_K_ratio | 1 | 0.99999901869484664 | 9.8130515333721968e-07 | 0 | null | -1.0352120739000901 | 0.3169719030182524 | 2.0746022011006527e-76 | 2.0744277097679295e-76 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 8 | r90 | 18.0625 | 18 | 0.0034602076124567475 | 0.0034602076124567475 | 1 | -1 | 0.33317013591547739 | 1.4301148022486204e-22 | 1.9701065876166452e-19 | PASS |
| 8 | r99 | 20.0625 | 20 | 0.0031152647975077881 | 0.0031152647975077881 | 1 | -1 | 0.33317013591547739 | 6.5149803197211003e-25 | 5.0766352092358305e-21 | PASS |
| 8 | radial_profile_L2 | null | null | 0.067252396223438593 | 0.031674679515173401 | 2.1232226261743889 | -1.1211565613575198 | 0.2798500299122274 | 1.1233733444798528e-17 | 3.6201109868020346e-15 | PASS |
| 12 | total_mass | 256 | 255.99987438921451 | 4.9066713076612034e-07 | 0 | null | -1.035242772400685 | 0.31695803351978358 | 6.3310511518897591e-81 | 6.3307848913142876e-81 | PASS |
| 12 | r_K_ratio | 1 | 0.99999901869489838 | 9.8130510160082673e-07 | 0 | null | -1.0352120282386914 | 0.31697192364827376 | 2.0746019330439565e-76 | 2.0744274417429762e-76 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 12 | r90 | 18.6875 | 18 | 0.03678929765886288 | 0 | null | -5.7445626465380286 | 3.8761814875216441e-05 | 1.2411344634416193e-17 | 3.0051609001276733e-16 | PASS |
| 12 | r99 | 20.5625 | 20.125 | 0.021276595744680851 | 0.015197568389057751 | 1.3999999999999999 | -3.415650255319866 | 0.0038327022878305922 | 3.2453744019125243e-19 | 3.5829058514479523e-17 | PASS |
| 12 | radial_profile_L2 | null | null | 0.072893691384749507 | 0.033001760468462815 | 2.2087819058746359 | -0.25736842332843862 | 0.80039152143350745 | 4.3685958447085025e-20 | 2.3969087135186269e-17 | PASS |
| 16 | total_mass | 256 | 256.00565975007765 | 2.2108398740866564e-05 | 0 | null | 29.0549503061726 | 1.3336686824631748e-14 | 7.6725742637134177e-78 | 7.6871276199674024e-78 | PASS |
| 16 | r_K_ratio | 1 | 1.0000442168399475 | 4.4216839947416875e-05 | 0 | null | 29.054978697547448 | 1.3336494557056482e-14 | 2.5117676850929536e-73 | 2.5213053845686074e-73 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.0047619047619047623 | 1 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 16 | r90 | 19.0625 | 18 | 0.055737704918032788 | 0.0032786885245901639 | 17 | -17 | 3.2769117377966068e-11 | 2.2623006381054966e-21 | 5.3331248280097439e-21 | PASS |
| 16 | r99 | 21.1875 | 20.8125 | 0.017699115044247787 | 0.014749262536873156 | 1.2 | -2.4227185592617446 | 0.028528068449741654 | 4.6950640203719188e-18 | 3.0093626661365215e-16 | PASS |
| 16 | radial_profile_L2 | null | null | 0.08775376074242433 | 0.029422120612765671 | 2.9825776971476938 | 0.20757506469588111 | 0.83835262262231791 | 5.9173726605820711e-20 | 1.2294617205940524e-16 | PASS |
| 20 | total_mass | 256 | 259.38930481728698 | 0.013239471942527239 | 0 | null | 48.255029981106148 | 7.1535576417554098e-18 | 9.9612836898674566e-40 | 3.1001936901933769e-39 | PASS |
| 20 | r_K_ratio | 1 | 1.0264789234962126 | 0.026478923496212517 | 0 | null | 48.255010551841188 | 7.1536005859322785e-18 | 1.9075543162085774e-35 | 1.8535898820419411e-34 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0 | null | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 20 | r90 | 19.375 | 18.125 | 0.064516129032258063 | 0 | null | -11.180339887498949 | 1.129731174656377e-08 | 6.9886341424657286e-17 | 1.1645945213682297e-17 | PASS |
| 20 | r99 | 21.9375 | 21 | 0.042735042735042736 | 0 | null | -6.5361700549589266 | 9.4202794798815216e-06 | 3.6683425778815099e-19 | 5.1861867803112176e-17 | PASS |
| 20 | radial_profile_L2 | null | null | 0.092421060738694855 | 0.040888641583758337 | 2.2603113519771729 | 1.0369297345499009 | 0.31619654317416779 | 1.1226962030547471e-19 | 3.5490613354533042e-16 | PASS |
| 24 | total_mass | 284.625 | 298.98754644036325 | 0.050461296233160362 | 0.0083443126921387799 | 6.047387974889797 | 19.730067038543012 | 3.8414531453545918e-12 | 1.6655638915222936e-27 | 4.0384969275713948e-22 | PASS |
| 24 | r_K_ratio | 1.2236328125 | 1.3358250453136638 | 0.091687826337742723 | 0.015163607342378291 | 6.0465708632206123 | 19.727446316178806 | 3.8488472047269843e-12 | 2.9675877127982947e-24 | 2.7669535033326392e-17 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.5625 | 13.0625 | 0.036866359447004608 | 0.018433179723502304 | 2 | -3.8729833462074166 | 0.0015017735386984323 | 9.6072681659644597e-15 | 7.6467248765302002e-14 | PASS |
| 24 | r90 | 19.9375 | 18.25 | 0.084639498432601878 | 0.003134796238244514 | 27 | -14.100290132411525 | 4.6333644836412128e-10 | 7.7245411315955695e-15 | 9.481589102890017e-19 | PASS |
| 24 | r99 | 22.875 | 21.25 | 0.071038251366120214 | 0.0054644808743169399 | 13 | -6.7890285822722154 | 6.1053729318667233e-06 | 1.2337518345797697e-14 | 7.416708381213179e-15 | PASS |
| 24 | radial_profile_L2 | null | null | 0.1311571822925125 | 0.031998880089473554 | 4.0988053933693243 | 4.1375477785895995 | 0.00087694673719733397 | 3.3174571066196259e-19 | 3.9701746118309766e-14 | PASS |
| 28 | total_mass | 360.0625 | 355.79428911997218 | 0.011854083332831953 | 0.0034716195105016492 | 3.4145687040222441 | -4.9981373855468707 | 0.00015894216488354686 | 4.9346299105874915e-27 | 1.1189772562667494e-23 | PASS |
| 28 | r_K_ratio | 1.81298828125 | 1.7783820031380952 | 0.019087976723183724 | 0.0053864799353622404 | 3.543682878659058 | -5.1899130998290026 | 0.00010985368152231689 | 4.8401929701778768e-24 | 5.9660093150193663e-21 | PASS |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.875 | 13.0625 | 0.058558558558558557 | 0.022522522522522521 | 2.6000000000000001 | -8.0622577482985509 | 7.8258249630702749e-07 | 2.2778649613890186e-15 | 2.5893641111623486e-16 | PASS |
| 28 | r90 | 20.5625 | 18.75 | 0.088145896656534953 | 0.0060790273556231003 | 14.5 | -11.066875266096559 | 1.2961577489987425e-08 | 1.5680800107416445e-13 | 1.9477056267575937e-16 | PASS |
| 28 | r99 | 23.625 | 21.9375 | 0.071428571428571425 | 0 | null | -8.5098307000225546 | 3.9910679004152207e-07 | 1.2402986697836266e-16 | 5.9545457575235637e-16 | PASS |
| 28 | radial_profile_L2 | null | null | 0.16812366489946665 | 0.025528598153441336 | 6.5856990614583761 | 8.2179631045461328 | 6.1744432441080883e-07 | 2.0668183461795087e-20 | 1.1633368762049021e-13 | PASS |
| 32 | total_mass | 385.5625 | 379.14312167383497 | 0.016649384538602752 | 0.0024315124007132437 | 6.8473368812426925 | -5.396509588004732 | 7.4173402428140978e-05 | 1.4730519212317234e-25 | 2.2336216332001333e-22 | PASS |
| 32 | r_K_ratio | 1.8212005047589639 | 1.9391831914335063 | 0.064782920038866071 | 0.0014984898683687194 | 43.232137504799965 | 14.726328885627707 | 2.5168461149834619e-10 | 2.6740318924724625e-25 | 4.4405753482401013e-18 | PASS |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.022421524663677129 | 2.7999999999999998 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 32 | r90 | 20.8125 | 18.9375 | 0.090090090090090086 | 0.006006006006006006 | 15 | -21.957751641341996 | 8.1207854767429677e-13 | 2.5303126792687863e-18 | 2.9422487852212645e-20 | PASS |
| 32 | r99 | 24.375 | 22 | 0.097435897435897437 | 0.0051282051282051282 | 19 | -15.343884208657718 | 1.4090734244550085e-10 | 8.3953741872725587e-18 | 3.7020237911275972e-18 | PASS |
| 32 | radial_profile_L2 | null | null | 0.17938424840949052 | 0.026253620259367155 | 6.8327433183424358 | 10.818675251643851 | 1.7576228289281234e-08 | 2.8179852830823343e-21 | 5.7756288571017507e-14 | PASS |
| 36 | total_mass | 404.5 | 389.43826184084367 | 0.037235446623377903 | 0.001854140914709518 | 20.082317545541816 | -11.824525086757328 | 5.2854372155451087e-09 | 1.4574793360891194e-25 | 2.4656946513950525e-22 | PASS |
| 36 | r_K_ratio | 1.6452597661040846 | 1.8931472248902017 | 0.15066767199510719 | 0.0025243036803777237 | 59.686824991104899 | 18.653271937262037 | 8.6419882044172036e-12 | 9.373935059396823e-23 | 6.9902926261365534e-12 | PASS |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0089686098654708519 | 7 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 36 | r90 | 21.3125 | 19 | 0.10850439882697947 | 0 | null | -19.322619811082461 | 5.195040620846464e-12 | 7.7864380968352111e-16 | 1.5335258890741899e-18 | PASS |
| 36 | r99 | 24.625 | 22.25 | 0.096446700507614211 | 0.0025380710659898475 | 38 | -11.783299785974805 | 5.5430657455111653e-09 | 3.3901457391067851e-15 | 5.0500172350014425e-17 | PASS |
| 36 | radial_profile_L2 | null | null | 0.19113782511778182 | 0.021200415960106734 | 9.0157582510385588 | 12.937967455874144 | 1.5352977834313781e-09 | 2.3205167683510447e-21 | 1.873559415950491e-13 | PASS |
| 40 | total_mass | 421.125 | 406.62557094179493 | 0.034430226318088571 | 0.0014841199168892847 | 23.199086493128078 | -12.123072892031722 | 3.7595679941645351e-09 | 1.454436824969024e-26 | 8.6049179868181217e-23 | PASS |
| 40 | r_K_ratio | 1.4862818219349081 | 1.6995080743478872 | 0.1434628677187155 | 0.0052088165665525443 | 27.542315204558335 | 16.515053506398107 | 4.9544479197942529e-11 | 1.9779308949787462e-22 | 1.6280039013798181e-11 | PASS |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0044843049327354259 | 14 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 40 | r90 | 21.5 | 19.0625 | 0.11337209302325581 | 0 | null | -19.030051422496395 | 6.4760803247711866e-12 | 6.3605411260372649e-15 | 2.3244786420838592e-18 | PASS |
| 40 | r99 | 25 | 22.6875 | 0.092499999999999999 | 0.0074999999999999997 | 12.333333333333334 | -15.363413772941893 | 1.3839300528290396e-10 | 5.4380328870719136e-17 | 4.8751207713627678e-19 | PASS |
| 40 | radial_profile_L2 | null | null | 0.20294210478022673 | 0.019135549966259355 | 10.605501547541785 | 12.331063058139474 | 2.9774061335968384e-09 | 2.6639476681283848e-22 | 9.5858061539418275e-14 | PASS |
| 44 | total_mass | 433 | 439.12414900465899 | 0.014143531188588849 | 0.0072170900692840644 | 1.9597276814908708 | 4.3124362386774742 | 0.00061619879125071012 | 5.6783101501054603e-26 | 3.3785584158268306e-21 | PASS |
| 44 | r_K_ratio | 1.3851122690599202 | 1.4826505247964381 | 0.070419025168780944 | 0.023883569861476831 | 2.948429634983663 | 10.150700193299443 | 4.1027820504537489e-08 | 1.7459206768544517e-22 | 6.5239376359162711e-15 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 14.0625 | 13.0625 | 0.071111111111111111 | 0.0088888888888888889 | 8 | -7.7459666924148332 | 1.2783477850604916e-06 | 1.25956553381152e-13 | 5.7635238603637946e-15 | PASS |
| 44 | r90 | 21.875 | 19.1875 | 0.12285714285714286 | 0.0057142857142857143 | 21.5 | -22.456017618285021 | 5.8544703070058476e-13 | 4.6040517988789706e-14 | 9.7976516786557512e-20 | PASS |
| 44 | r99 | 25.4375 | 22.9375 | 0.098280098280098274 | 0.0024570024570024569 | 40 | -15.811388300841896 | 9.2073297610083179e-11 | 1.32240574860815e-17 | 2.0753132400270518e-18 | PASS |
| 44 | radial_profile_L2 | null | null | 0.20578903863170117 | 0.022666905945975712 | 9.0788323347782285 | 13.563001661788283 | 7.9734125142345473e-10 | 3.9888272112335183e-23 | 2.1411534627969286e-14 | PASS |
| 48 | total_mass | 468.375 | 495.63025101823672 | 0.058191088376272695 | 0.0066720042700827327 | 8.7216803258357523 | 13.932439671958022 | 5.478965023639955e-10 | 6.2262438154028136e-25 | 9.8983062404198524e-19 | PASS |
| 48 | r_K_ratio | 1.4233832881828432 | 1.4168010867127527 | 0.00462433521928831 | 0.023969344303723984 | 0.19292706386506586 | -0.56934483893448984 | 0.5775492092836918 | 5.7033749883260893e-20 | 1.7002683649562875e-15 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 14.4375 | 13.1875 | 0.086580086580086577 | 0.030303030303030304 | 2.8571428571428572 | -11.180339887498949 | 1.129731174656377e-08 | 5.7520127976418164e-14 | 2.7873182537487998e-16 | PASS |
| 48 | r90 | 22.375 | 19.5 | 0.12849162011173185 | 0 | null | -15.998991903725805 | 7.7865353236989104e-11 | 1.0061554522861548e-11 | 4.8264550809103502e-17 | PASS |
| 48 | r99 | 26.1875 | 23.1875 | 0.11455847255369929 | 0.0071599045346062056 | 16 | -23.2379000772445 | 3.5513304925948273e-13 | 1.7995348202025614e-17 | 8.6169881800138773e-21 | PASS |
| 48 | radial_profile_L2 | null | null | 0.23341720836842744 | 0.017773720507827693 | 13.1327151378142 | 12.748408916387797 | 1.8827621113569947e-09 | 3.1119166670442979e-21 | 6.3693627289909649e-11 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1875, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

#### no_division: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 255.99987439227883 | 4.9065516088270256e-07 | 0 | null | -1.0352175095070719 | 0.31696944719393905 | 6.3310518845588904e-81 | 6.3307856304481416e-81 | PASS |
| 4 | r_K_ratio | 1 | 0.9999990186955926 | 9.8130440741306391e-07 | 0 | null | -1.0352112711393178 | 0.31697226570925685 | 2.0746026777977294e-76 | 2.0744281865575098e-76 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13.0625 | 13.0625 | 0 | 0.0047846889952153108 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 4 | r90 | 17.9375 | 17.6875 | 0.013937282229965157 | 0.0069686411149825784 | 2 | -2.2360679774997898 | 0.040968955955836141 | 1.1148111900685798e-16 | 3.8956480726576652e-17 | PASS |
| 4 | r99 | 19.25 | 19.375 | 0.0064935064935064939 | 0.00974025974025974 | 0.66666666666666663 | 1 | 0.33317013591547739 | 1.2617854706008803e-18 | 7.1507823999893397e-17 | PASS |
| 4 | radial_profile_L2 | null | null | 0.01390051111786134 | 0.026445474852675486 | 0.52562909894034393 | -7.7751682619068863 | 1.2210780923911863e-06 | 1.7063425300474121e-19 | 5.5915970512964533e-19 | PASS |
| 8 | total_mass | 256 | 255.99987439079101 | 4.9066097267125297e-07 | 0 | null | -1.0352297764902778 | 0.31696390498294125 | 6.331051438669818e-81 | 6.3307851814239781e-81 | PASS |
| 8 | r_K_ratio | 1 | 0.99999901869484664 | 9.8130515333721968e-07 | 0 | null | -1.0352120739000901 | 0.3169719030182524 | 2.0746022011006527e-76 | 2.0744277097679295e-76 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 8 | r90 | 18.0625 | 18 | 0.0034602076124567475 | 0.0034602076124567475 | 1 | -1 | 0.33317013591547739 | 1.4301148022486204e-22 | 1.9701065876166452e-19 | PASS |
| 8 | r99 | 20.0625 | 20 | 0.0031152647975077881 | 0.0031152647975077881 | 1 | -1 | 0.33317013591547739 | 6.5149803197211003e-25 | 5.0766352092358305e-21 | PASS |
| 8 | radial_profile_L2 | null | null | 0.067252396223438593 | 0.031674679515173401 | 2.1232226261743889 | -1.1211565613575198 | 0.2798500299122274 | 1.1233733444798528e-17 | 3.6201109868020346e-15 | PASS |
| 12 | total_mass | 256 | 255.99987438904748 | 4.9066778327194749e-07 | 0 | null | -1.0352441528078298 | 0.31695740986592646 | 6.3310508120942005e-81 | 6.3307845511789133e-81 | PASS |
| 12 | r_K_ratio | 1 | 0.99999901869359353 | 9.8130640650839762e-07 | 0 | null | -1.0352134085283748 | 0.31697130002772445 | 2.0746018219831146e-76 | 2.0744273304594725e-76 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 12 | r90 | 18.6875 | 18 | 0.03678929765886288 | 0 | null | -5.7445626465380286 | 3.8761814875216441e-05 | 1.2411344634416193e-17 | 3.0051609001276733e-16 | PASS |
| 12 | r99 | 20.5625 | 20.125 | 0.021276595744680851 | 0.015197568389057751 | 1.3999999999999999 | -3.415650255319866 | 0.0038327022878305922 | 3.2453744019125243e-19 | 3.5829058514479523e-17 | PASS |
| 12 | radial_profile_L2 | null | null | 0.072893691384777443 | 0.033001760468462815 | 2.2087819058754823 | -0.2573684233257047 | 0.80039152143557923 | 4.3685958447461105e-20 | 2.3969087135451209e-17 | PASS |
| 16 | total_mass | 256 | 255.99987438712913 | 4.9067527684715229e-07 | 0 | null | -1.0352599704717151 | 0.31695026367018198 | 6.3310501526611944e-81 | 6.3307838877075447e-81 | PASS |
| 16 | r_K_ratio | 1 | 0.99999901869249208 | 9.813075079051492e-07 | 0 | null | -1.0352145828292734 | 0.3169707694744538 | 2.0746014493280099e-76 | 2.0744269576397917e-76 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.0047619047619047623 | 1 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 16 | r90 | 19.0625 | 18 | 0.055737704918032788 | 0.0032786885245901639 | 17 | -17 | 3.2769117377966068e-11 | 2.2623006381054966e-21 | 5.3331248280097439e-21 | PASS |
| 16 | r99 | 21.1875 | 20.8125 | 0.017699115044247787 | 0.014749262536873156 | 1.2 | -2.4227185592617446 | 0.028528068449741654 | 4.6950640203719188e-18 | 3.0093626661365215e-16 | PASS |
| 16 | radial_profile_L2 | null | null | 0.087754755358269959 | 0.029422120612765671 | 2.9826115021836639 | 0.20764347161909227 | 0.83830016142821528 | 5.9183589942778144e-20 | 1.2297761744011539e-16 | PASS |
| 20 | total_mass | 256 | 255.99987438505411 | 4.9068338240504383e-07 | 0 | null | -1.0352770727410276 | 0.31694253723903631 | 6.3310500985700103e-81 | 6.3307838292200258e-81 | PASS |
| 20 | r_K_ratio | 1 | 0.9999990186911859 | 9.8130881408947657e-07 | 0 | null | -1.0352159724376802 | 0.31697014164532411 | 2.0746010986317065e-76 | 2.0744266067408639e-76 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0 | null | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 20 | r90 | 19.375 | 18.0625 | 0.067741935483870974 | 0 | null | -10.966892325208963 | 1.464381425614242e-08 | 1.4279603269384565e-16 | 3.1631854909452851e-17 | PASS |
| 20 | r99 | 21.9375 | 21 | 0.042735042735042736 | 0 | null | -6.5361700549589266 | 9.4202794798815216e-06 | 3.6683425778815099e-19 | 5.1861867803112176e-17 | PASS |
| 20 | radial_profile_L2 | null | null | 0.093009594440941076 | 0.040888641583758337 | 2.2747049263159203 | 1.0690581454497692 | 0.30194617124546341 | 1.1805805401977152e-19 | 3.9362528250669236e-16 | PASS |
| 24 | total_mass | 256 | 255.99987438281661 | 4.9069212255253847e-07 | 0 | null | -1.0352955151313261 | 0.31693420552315876 | 6.3310499326699688e-81 | 6.3307836585843207e-81 | PASS |
| 24 | r_K_ratio | 1 | 0.99999901868982044 | 9.8131017959440792e-07 | 0 | null | -1.035217407255228 | 0.31696949339151714 | 2.0746012700957044e-76 | 2.0744267779475894e-76 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.3125 | 13.0625 | 0.018779342723004695 | 0.0046948356807511738 | 4 | -2.2360679774997898 | 0.040968955955836141 | 3.9747717054621687e-16 | 3.6766171447656853e-14 | PASS |
| 24 | r90 | 19.875 | 18.125 | 0.088050314465408799 | 0.0031446540880503146 | 28 | -15.652475842498529 | 1.0626678392473952e-10 | 1.2944120883335444e-15 | 8.8704215647820619e-19 | PASS |
| 24 | r99 | 22.3125 | 21.0625 | 0.056022408963585436 | 0.011204481792717087 | 5 | -8.6602540378443873 | 3.1998825433805415e-07 | 3.6206560128381799e-18 | 1.0967931003213516e-17 | PASS |
| 24 | radial_profile_L2 | null | null | 0.11673751969066783 | 0.022394649227381803 | 5.2127416020400785 | 3.787816485567411 | 0.0017872015439631119 | 3.2425496813676196e-19 | 9.700628894214686e-15 | PASS |
| 28 | total_mass | 256 | 255.9998743805387 | 4.9070102071946398e-07 | 0 | null | -1.0353142948070588 | 0.31692572159463139 | 6.3310494104320913e-81 | 6.3307831315397614e-81 | PASS |
| 28 | r_K_ratio | 1 | 0.99999901868788166 | 9.8131211836299803e-07 | 0 | null | -1.0352194748632058 | 0.31696855924324929 | 2.0746005988649749e-76 | 2.0744261064286434e-76 | PASS |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.625 | 13.0625 | 0.041284403669724773 | 0.01834862385321101 | 2.25 | -4.3915503282683988 | 0.00052573078993947925 | 1.1808157431972854e-14 | 4.8233148568867511e-14 | PASS |
| 28 | r90 | 20 | 18.375 | 0.081250000000000003 | 0 | null | -13 | 1.4369304143959042e-09 | 1.9652481154345611e-14 | 8.8277683744579181e-19 | PASS |
| 28 | r99 | 22.5625 | 21.625 | 0.041551246537396121 | 0.01662049861495845 | 2.5 | -5.5141096657035584 | 5.9461335275977702e-05 | 3.3590340125425255e-17 | 1.6099581824988478e-16 | PASS |
| 28 | radial_profile_L2 | null | null | 0.13978619672757553 | 0.018351607168479032 | 7.6171092506641154 | 5.4207173322569888 | 7.0862977256131812e-05 | 8.5769360995731561e-18 | 2.2252589260061686e-12 | PASS |
| 32 | total_mass | 256 | 255.99987437811777 | 4.9071047747428764e-07 | 0 | null | -1.0353342445328497 | 0.31691670926278337 | 6.3310496677110311e-81 | 6.3307833836764088e-81 | PASS |
| 32 | r_K_ratio | 1 | 0.99999901868601937 | 9.813139805817106e-07 | 0 | null | -1.0352214614519035 | 0.31696766170148438 | 2.074599935572742e-76 | 2.0744254428610281e-76 | PASS |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.6875 | 13.0625 | 0.045662100456621002 | 0 | null | -5 | 0.0001583695146220272 | 1.1779091896870362e-14 | 2.3594520889093847e-14 | PASS |
| 32 | r90 | 20.125 | 18.75 | 0.068322981366459631 | 0 | null | -11 | 1.4062516106729139e-08 | 2.0258582724161173e-15 | 5.5280594142612277e-18 | PASS |
| 32 | r99 | 23.1875 | 22 | 0.051212938005390833 | 0.0080862533692722376 | 6.333333333333333 | -8.7331326357615122 | 2.8777724060903971e-07 | 1.1912985689352152e-19 | 7.1361375974270531e-18 | PASS |
| 32 | radial_profile_L2 | null | null | 0.15090808493996466 | 0.027516323573992162 | 5.4843113228469136 | 7.6685834092872005 | 1.4441975658565469e-06 | 1.069863440748573e-20 | 9.8190951434487498e-15 | PASS |
| 36 | total_mass | 256 | 255.99987437556251 | 4.9072045889969607e-07 | 0 | null | -1.0353553206898731 | 0.31690718826448627 | 6.331048141512166e-81 | 6.3307818521254335e-81 | PASS |
| 36 | r_K_ratio | 1 | 0.9999990186842691 | 9.8131573089688118e-07 | 0 | null | -1.035223309200515 | 0.31696682688938971 | 2.0745998972012456e-76 | 2.0744254041815161e-76 | PASS |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0 | null | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 36 | r90 | 20.3125 | 19 | 0.064615384615384616 | 0.0030769230769230769 | 21 | -10.966892325208963 | 1.464381425614242e-08 | 2.8767227697700644e-17 | 2.1976378153234503e-17 | PASS |
| 36 | r99 | 23.8125 | 22 | 0.076115485564304461 | 0.0052493438320209973 | 14.5 | -8.6913130282879933 | 3.0581888927239612e-07 | 2.2591980500354052e-16 | 9.9911192352252567e-16 | PASS |
| 36 | radial_profile_L2 | null | null | 0.16069714407800084 | 0.024317400523262839 | 6.6083191714621217 | 9.8012618933072062 | 6.5030784298316884e-08 | 8.573087776919444e-22 | 2.2800715018192439e-15 | PASS |
| 40 | total_mass | 256 | 255.9998743730192 | 4.9073039378594308e-07 | 0 | null | -1.0353762665236637 | 0.31689772634366264 | 6.3310495625183636e-81 | 6.3307832676809295e-81 | PASS |
| 40 | r_K_ratio | 1 | 0.99999901868181862 | 9.8131818140889671e-07 | 0 | null | -1.0352259095778034 | 0.31696565204264782 | 2.0745994390357899e-76 | 2.0744249456188915e-76 | PASS |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 14 | 13.0625 | 0.066964285714285712 | 0.004464285714285714 | 15 | -15 | 1.9412759304364392e-10 | 4.0487734815626444e-17 | 1.1177178271930511e-20 | PASS |
| 40 | r90 | 20.6875 | 19 | 0.081570996978851965 | 0.0030211480362537764 | 27 | -14.100290132411525 | 4.6333644836412128e-10 | 9.12499431131016e-17 | 7.658417637110565e-18 | PASS |
| 40 | r99 | 23.8125 | 22.25 | 0.065616797900262466 | 0.007874015748031496 | 8.3333333333333339 | -7.0059809900903849 | 4.2376047084641228e-06 | 1.2550586234197359e-15 | 2.329922511027581e-15 | PASS |
| 40 | radial_profile_L2 | null | null | 0.17199804296440277 | 0.030353446894415725 | 5.6665077795842063 | 9.3102227955353438 | 1.2688299753401375e-07 | 4.9235453567793169e-23 | 4.6514123933925176e-16 | PASS |
| 44 | total_mass | 256 | 255.99987437034028 | 4.9074085824857283e-07 | 0 | null | -1.0353983711860357 | 0.31688774116346385 | 6.3310471779658376e-81 | 6.3307808775502395e-81 | PASS |
| 44 | r_K_ratio | 1 | 0.99999901868063235 | 9.8131936767525962e-07 | 0 | null | -1.0352271735700516 | 0.31696508097389442 | 2.0745990616009301e-76 | 2.0744245680048234e-76 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 14 | 13.0625 | 0.066964285714285712 | 0.004464285714285714 | 15 | -15 | 1.9412759304364392e-10 | 4.0487734815626444e-17 | 1.1177178271930511e-20 | PASS |
| 44 | r90 | 20.75 | 19 | 0.084337349397590355 | 0.0030120481927710845 | 28 | -15.652475842498529 | 1.0626678392473952e-10 | 4.0567440289856516e-17 | 2.3449936241056999e-18 | PASS |
| 44 | r99 | 24.25 | 22.6875 | 0.064432989690721643 | 0 | null | -8.5917924577933835 | 3.5372277238368397e-07 | 8.1352456598561994e-17 | 7.0640630267479829e-17 | PASS |
| 44 | radial_profile_L2 | null | null | 0.18402034179826446 | 0.02352354020837873 | 7.8228166410393953 | 10.43009788287987 | 2.8633323824094639e-08 | 4.6321364753720373e-22 | 1.671793096074856e-14 | PASS |
| 48 | total_mass | 256 | 255.99987436776033 | 4.907509362564455e-07 | 0 | null | -1.0354196343172224 | 0.31687813633677431 | 6.3310471929774437e-81 | 6.3307808870923723e-81 | PASS |
| 48 | r_K_ratio | 1 | 0.99999901867900443 | 9.8132099562303621e-07 | 0 | null | -1.0352289312075473 | 0.31696428687860334 | 2.0745978515295287e-76 | 2.0744233577457725e-76 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 14.125 | 13.0625 | 0.075221238938053103 | 0 | null | -17 | 3.2769117377966068e-11 | 2.8517534137653066e-18 | 2.0670478007634922e-19 | PASS |
| 48 | r90 | 21.0625 | 19.0625 | 0.094955489614243327 | 0 | null | -21.908902300206645 | 8.3886833819515926e-13 | 4.8178999609525139e-17 | 1.2775397916135412e-20 | PASS |
| 48 | r99 | 24.625 | 23 | 0.065989847715736044 | 0.0025380710659898475 | 26 | -9.0429084673232811 | 1.845613926270835e-07 | 7.7618792499716411e-18 | 1.000515420358179e-16 | PASS |
| 48 | radial_profile_L2 | null | null | 0.19551846661883626 | 0.01345602297659409 | 14.53018228037574 | 13.029636004631103 | 1.3923171074832122e-09 | 1.0973254564531445e-22 | 1.5977754381756672e-14 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.203125, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.203125, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1796875, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

#### no_exchange: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 255.99987439227883 | 4.9065516088270256e-07 | 0 | null | -1.0352175095070719 | 0.31696944719393905 | 6.3310518845588904e-81 | 6.3307856304481416e-81 | PASS |
| 4 | r_K_ratio | 1 | 0.9999990186955926 | 9.8130440741306391e-07 | 0 | null | -1.0352112711393178 | 0.31697226570925685 | 2.0746026777977294e-76 | 2.0744281865575098e-76 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13.0625 | 13.0625 | 0 | 0.0047846889952153108 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 4 | r90 | 17.9375 | 17.6875 | 0.013937282229965157 | 0.0069686411149825784 | 2 | -2.2360679774997898 | 0.040968955955836141 | 1.1148111900685798e-16 | 3.8956480726576652e-17 | PASS |
| 4 | r99 | 19.25 | 19.375 | 0.0064935064935064939 | 0.00974025974025974 | 0.66666666666666663 | 1 | 0.33317013591547739 | 1.2617854706008803e-18 | 7.1507823999893397e-17 | PASS |
| 4 | radial_profile_L2 | null | null | 0.01390051111786134 | 0.026445474852675486 | 0.52562909894034393 | -7.7751682619068863 | 1.2210780923911863e-06 | 1.706342530047388e-19 | 5.5915970512963349e-19 | PASS |
| 8 | total_mass | 256 | 255.99987439079101 | 4.9066097267125297e-07 | 0 | null | -1.0352297764902778 | 0.31696390498294125 | 6.331051438669818e-81 | 6.3307851814239781e-81 | PASS |
| 8 | r_K_ratio | 1 | 0.99999901869484664 | 9.8130515333721968e-07 | 0 | null | -1.0352120739000901 | 0.3169719030182524 | 2.0746022011006527e-76 | 2.0744277097679295e-76 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 8 | r90 | 18.0625 | 18 | 0.0034602076124567475 | 0.0034602076124567475 | 1 | -1 | 0.33317013591547739 | 1.4301148022486204e-22 | 1.9701065876166452e-19 | PASS |
| 8 | r99 | 20.0625 | 20 | 0.0031152647975077881 | 0.0031152647975077881 | 1 | -1 | 0.33317013591547739 | 6.5149803197211003e-25 | 5.0766352092358305e-21 | PASS |
| 8 | radial_profile_L2 | null | null | 0.066248819732666941 | 0.030175108540665053 | 2.1954790864592142 | -1.2820130111540684 | 0.21929912435332344 | 8.195558640542962e-18 | 2.4247462922964637e-15 | PASS |
| 12 | total_mass | 256 | 255.99987438921451 | 4.9066713076612034e-07 | 0 | null | -1.035242772400685 | 0.31695803351978358 | 6.3310511518897591e-81 | 6.3307848913142876e-81 | PASS |
| 12 | r_K_ratio | 1 | 0.99999901869489838 | 9.8130510160082673e-07 | 0 | null | -1.0352120282386914 | 0.31697192364827376 | 2.0746019330439565e-76 | 2.0744274417429762e-76 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 12 | r90 | 18.6875 | 18 | 0.03678929765886288 | 0 | null | -5.7445626465380286 | 3.8761814875216441e-05 | 1.2411344634416193e-17 | 3.0051609001276733e-16 | PASS |
| 12 | r99 | 20.625 | 20.125 | 0.024242424242424242 | 0.012121212121212121 | 2 | -3.8729833462074166 | 0.0015017735386984323 | 4.6661893853944489e-19 | 3.0166728576555097e-17 | PASS |
| 12 | radial_profile_L2 | null | null | 0.071941926817638585 | 0.03591453299254424 | 2.0031424836450897 | -0.7838739412885124 | 0.4453164250255045 | 2.436187300926997e-19 | 1.2180140854272305e-16 | PASS |
| 16 | total_mass | 256 | 256.00565975007765 | 2.2108398740866564e-05 | 0 | null | 29.0549503061726 | 1.3336686824631748e-14 | 7.6725742637134177e-78 | 7.6871276199674024e-78 | PASS |
| 16 | r_K_ratio | 1 | 1.0000442168399475 | 4.4216839947416875e-05 | 0 | null | 29.054978697547448 | 1.3336494557056482e-14 | 2.5117676850929536e-73 | 2.5213053845686074e-73 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13.0625 | 13.0625 | 0 | 0.0047846889952153108 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 16 | r90 | 19.125 | 18 | 0.058823529411764705 | 0.0098039215686274508 | 6 | -13.174650984805199 | 1.1942408086569779e-09 | 2.909519278220455e-19 | 4.6398892759804035e-19 | PASS |
| 16 | r99 | 21.125 | 20.8125 | 0.014792899408284023 | 0.017751479289940829 | 0.83333333333333337 | -2.0761369963434992 | 0.055487124843606725 | 3.6381493848699621e-18 | 1.9476262532447893e-16 | PASS |
| 16 | radial_profile_L2 | null | null | 0.086926891747487897 | 0.032450865897919859 | 2.6787233357942544 | 0.014402595107158356 | 0.98869860927762288 | 1.3703854872143146e-19 | 2.627664689266668e-16 | PASS |
| 20 | total_mass | 256 | 259.38930481728698 | 0.013239471942527239 | 0 | null | 48.255029981106148 | 7.1535576417554098e-18 | 9.9612836898674566e-40 | 3.1001936901933769e-39 | PASS |
| 20 | r_K_ratio | 1 | 1.0264789234962126 | 0.026478923496212517 | 0 | null | 48.255010551841188 | 7.1536005859322785e-18 | 1.9075543162085774e-35 | 1.8535898820419411e-34 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13.25 | 13.0625 | 0.014150943396226415 | 0.009433962264150943 | 1.5 | -1.8605210188381269 | 0.082530706346961857 | 6.9722898961555865e-17 | 1.1503841893528265e-14 | PASS |
| 20 | r90 | 19.4375 | 18.125 | 0.067524115755627015 | 0.0032154340836012861 | 21 | -10.966892325208963 | 1.464381425614242e-08 | 2.3619300207992757e-16 | 2.4732173233321653e-17 | PASS |
| 20 | r99 | 21.9375 | 21 | 0.042735042735042736 | 0.0028490028490028491 | 15 | -6.5361700549589266 | 9.4202794798815216e-06 | 3.6683425778815099e-19 | 5.1861867803112176e-17 | PASS |
| 20 | radial_profile_L2 | null | null | 0.093679748498929849 | 0.041094034480661526 | 2.2796434977202025 | 1.1187438414798105 | 0.28084589490878142 | 3.0281526611303621e-19 | 1.0651740086217456e-15 | PASS |
| 24 | total_mass | 284.8125 | 298.98754644036325 | 0.049769748309372898 | 0.008558262014483212 | 5.8154036678413403 | 20.0361625079031 | 3.0736708183358354e-12 | 7.7789323162985945e-28 | 3.426427010845314e-22 | PASS |
| 24 | r_K_ratio | 1.22509765625 | 1.3358250453136638 | 0.090382500120519446 | 0.015544041450777202 | 5.8146075077534176 | 20.033500066322517 | 3.0795956745168783e-12 | 1.4010009851764736e-24 | 2.2381843671310074e-17 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.4375 | 13.0625 | 0.027906976744186046 | 0.013953488372093023 | 2 | -3 | 0.0089727374772233335 | 3.3710806931812163e-15 | 9.4193184812143518e-14 | PASS |
| 24 | r90 | 19.875 | 18.25 | 0.081761006289308172 | 0.0062893081761006293 | 13 | -13 | 1.4369304143959042e-09 | 7.5808473075664949e-15 | 3.5959184396777153e-18 | PASS |
| 24 | r99 | 22.875 | 21.25 | 0.071038251366120214 | 0 | null | -7.3441248772584169 | 2.4289023849674602e-06 | 4.9867352619398684e-15 | 1.9987441079264612e-15 | PASS |
| 24 | radial_profile_L2 | null | null | 0.12342551612005506 | 0.038196622833248067 | 3.2313201263599649 | 2.4870774880217148 | 0.025138508828478491 | 8.3372720494273021e-19 | 4.6527836070042218e-14 | PASS |
| 28 | total_mass | 360.625 | 355.79428911997218 | 0.013395385455882996 | 0.0057192374350086657 | 2.342162850922572 | -6.0011359669318205 | 2.4276682640427044e-05 | 1.0865789711141528e-27 | 7.4058362084402783e-24 | PASS |
| 28 | r_K_ratio | 1.8173828125 | 1.7783820031380952 | 0.021459875758511864 | 0.0088662009672219235 | 2.4204138658539747 | -6.2046064801154497 | 1.6858623535707178e-05 | 1.0925975605988282e-24 | 3.780937455022612e-21 | PASS |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.75 | 13.0625 | 0.050000000000000003 | 0.022727272727272728 | 2.2000000000000002 | -5.7445626465380286 | 3.8761814875216441e-05 | 9.3244950117320688e-15 | 8.5085128194990869e-15 | PASS |
| 28 | r90 | 20.4375 | 18.75 | 0.082568807339449546 | 0.0030581039755351682 | 27 | -11.211139780254895 | 1.088573529846995e-08 | 3.0364248363029417e-14 | 8.4111151283276818e-17 | PASS |
| 28 | r99 | 23.625 | 21.9375 | 0.071428571428571425 | 0.010582010582010581 | 6.75 | -9.5859666337058052 | 8.6904680420919338e-08 | 9.7834698306601055e-18 | 1.5675213484275571e-16 | PASS |
| 28 | radial_profile_L2 | null | null | 0.16542603291152005 | 0.022094871907015326 | 7.4870781603850691 | 6.6942300207160192 | 7.1759868216089277e-06 | 2.6089760311678888e-21 | 1.1373539624592964e-14 | PASS |
| 32 | total_mass | 386.4375 | 379.14312167383497 | 0.018875958793246057 | 0.0035581432961345623 | 5.3050024190291083 | -9.8397143431062375 | 6.1779713220304871e-08 | 2.816241675263828e-28 | 1.2403157484159302e-25 | PASS |
| 32 | r_K_ratio | 1.8277699278440873 | 1.9391831914335063 | 0.060955846735499364 | 0.0010973701559265133 | 55.547206570452182 | 12.944742767156743 | 1.5242144419604462e-09 | 3.5859363142369425e-25 | 1.4315633171656693e-17 | PASS |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.6875 | 13.0625 | 0.045662100456621002 | 0 | null | -5 | 0.0001583695146220272 | 1.1779091896870362e-14 | 2.3594520889093847e-14 | PASS |
| 32 | r90 | 20.875 | 18.9375 | 0.092814371257485026 | 0.0029940119760479044 | 31 | -31 | 5.1219766458245595e-15 | 4.0733458333999718e-20 | 2.7256983675835709e-22 | PASS |
| 32 | r99 | 24.5 | 22 | 0.10204081632653061 | 0.01020408163265306 | 10 | -19.364916731037081 | 5.0333996972544363e-12 | 7.3215372064694177e-19 | 1.9251631615291149e-19 | PASS |
| 32 | radial_profile_L2 | null | null | 0.17790357500801085 | 0.039054048176036631 | 4.5553171391121392 | 10.194805911092661 | 3.8743810984394743e-08 | 6.7219953237852673e-22 | 1.1942699611155942e-14 | PASS |
| 36 | total_mass | 404.75 | 389.43826184084367 | 0.037830112808292432 | 0.0021618282890673254 | 17.499129324750129 | -13.602287848182099 | 7.6583666500409407e-10 | 1.8586717845383549e-26 | 4.2233180049179681e-23 | PASS |
| 36 | r_K_ratio | 1.6521944851880093 | 1.8931472248902017 | 0.14583800022475779 | 0.0065422591158011299 | 22.291688183447814 | 20.193452609479298 | 2.7442973758045753e-12 | 1.9030176432677248e-23 | 1.014820077938134e-12 | PASS |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.8125 | 13.0625 | 0.054298642533936653 | 0.0045248868778280547 | 12 | -6.7082039324993694 | 7.0065583016648931e-06 | 5.5716931266528315e-15 | 2.0366123748919692e-15 | PASS |
| 36 | r90 | 21.1875 | 19 | 0.10324483775811209 | 0.0029498525073746312 | 35 | -21.706078553111482 | 9.6059095504530756e-13 | 3.8511104687127252e-17 | 1.6037390378163572e-19 | PASS |
| 36 | r99 | 24.8125 | 22.25 | 0.10327455919395466 | 0 | null | -12.593049895096 | 2.2296860033886611e-09 | 4.9424412867465559e-15 | 4.56913906542305e-17 | PASS |
| 36 | radial_profile_L2 | null | null | 0.18882425299079525 | 0.022438997325317463 | 8.4150040330789952 | 11.420321507114894 | 8.4786085811858253e-09 | 2.0842269236263657e-21 | 1.2824564343926165e-13 | PASS |
| 40 | total_mass | 421.375 | 406.62557094179493 | 0.035003094768804623 | 0.0022248590922574903 | 15.732724328752052 | -11.89707781301192 | 4.8623705407176048e-09 | 1.8686176541286275e-26 | 1.6141577258952962e-22 | PASS |
| 40 | r_K_ratio | 1.4911902774804242 | 1.6995080743478872 | 0.13969900422060516 | 0.0094160123977586416 | 14.836323309627192 | 18.963542625932902 | 6.8118131341903514e-12 | 2.2258678112555633e-23 | 1.138963136627173e-12 | PASS |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 14 | 13.0625 | 0.066964285714285712 | 0.004464285714285714 | 15 | -15 | 1.9412759304364392e-10 | 4.0487734815626444e-17 | 1.1177178271930511e-20 | PASS |
| 40 | r90 | 21.75 | 19.0625 | 0.1235632183908046 | 0.014367816091954023 | 8.5999999999999996 | -22.456017618285021 | 5.8544703070058476e-13 | 9.3244950117320688e-15 | 3.5906941929051801e-19 | PASS |
| 40 | r99 | 25.25 | 22.6875 | 0.10148514851485149 | 0.0074257425742574254 | 13.666666666666666 | -16.291747991900039 | 6.0154840401030353e-11 | 1.1696589951698737e-16 | 7.8093950089525235e-19 | PASS |
| 40 | radial_profile_L2 | null | null | 0.19946607524360083 | 0.019535136155509093 | 10.210631431271059 | 11.446172163323196 | 8.2228699855027335e-09 | 7.1734444791093529e-22 | 1.6351820551398424e-13 | PASS |
| 44 | total_mass | 432.9375 | 439.12414900465899 | 0.014289935625024331 | 0.0051970550021654396 | 2.7496217798484315 | 4.9361790769697187 | 0.00017926066775840707 | 1.0981696089920881e-26 | 4.9990598504938405e-22 | PASS |
| 44 | r_K_ratio | 1.3857999867500879 | 1.4826505247964381 | 0.069887818568594079 | 0.026930385991516315 | 2.5951287363876006 | 9.1974342045595563 | 1.4847491192254041e-07 | 9.7527730764142976e-22 | 1.9621691026177413e-14 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 14 | 13.0625 | 0.066964285714285712 | 0.004464285714285714 | 15 | -15 | 1.9412759304364392e-10 | 4.0487734815626444e-17 | 1.1177178271930511e-20 | PASS |
| 44 | r90 | 21.875 | 19.1875 | 0.12285714285714286 | 0 | null | -22.456017618285021 | 5.8544703070058476e-13 | 4.6040517988789706e-14 | 9.7976516786557512e-20 | PASS |
| 44 | r99 | 25.625 | 22.9375 | 0.1048780487804878 | 0.012195121951219513 | 8.5999999999999996 | -22.456017618285021 | 5.8544703070058476e-13 | 3.6838863499834168e-19 | 2.3186371219820026e-20 | PASS |
| 44 | radial_profile_L2 | null | null | 0.20535251306838187 | 0.028784123839218512 | 7.1342283758725369 | 14.358272126628426 | 3.5929928341317586e-10 | 3.6781928483533143e-24 | 1.9359477112293512e-15 | PASS |
| 48 | total_mass | 469.1875 | 495.63025101823672 | 0.056358600811480961 | 0.0035966431330757961 | 15.669778381177318 | 12.268563975436868 | 3.1924773173098637e-09 | 3.0886071593648527e-24 | 3.1704042354516734e-18 | PASS |
| 48 | r_K_ratio | 1.4261711036759821 | 1.4168010867127527 | 0.0065700510542374248 | 0.028301913590415343 | 0.23214158411049712 | -0.74957653649342315 | 0.46510423827223757 | 2.1351143976585334e-19 | 4.3418128465380576e-15 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 14.375 | 13.1875 | 0.082608695652173908 | 0.0043478260869565218 | 19 | -8.7331326357615122 | 2.8777724060903971e-07 | 9.3936200515851003e-13 | 4.2247526237449935e-15 | PASS |
| 48 | r90 | 22.5 | 19.5 | 0.13333333333333333 | 0.0027777777777777779 | 48 | -18.973665961010276 | 6.7595384058058026e-12 | 4.3019444621415003e-12 | 3.8744608228872494e-18 | PASS |
| 48 | r99 | 26.1875 | 23.1875 | 0.11455847255369929 | 0 | null | -16.431676725154986 | 5.3252598643225297e-11 | 7.3932535647552742e-16 | 3.6474139884078859e-18 | PASS |
| 48 | radial_profile_L2 | null | null | 0.23448128323377226 | 0.023791703208090063 | 9.855590462898844 | 13.077557604554258 | 1.3232590908885747e-09 | 1.8490557852920067e-24 | 5.7468066170384352e-14 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1875, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

#### enlarged_boundary: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 255.99987439227883 | 4.9065516088270256e-07 | 0 | null | -1.0352175095070719 | 0.31696944719393905 | 6.3310518845588904e-81 | 6.3307856304481416e-81 | PASS |
| 4 | r_K_ratio | 1 | 0.9999990186955926 | 9.8130440741306391e-07 | 0 | null | -1.0352112711393178 | 0.31697226570925685 | 2.0746026777977294e-76 | 2.0744281865575098e-76 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13.0625 | 13.0625 | 0 | 0.0047846889952153108 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 4 | r90 | 17.9375 | 17.6875 | 0.013937282229965157 | 0.0069686411149825784 | 2 | -2.2360679774997898 | 0.040968955955836141 | 1.1148111900685798e-16 | 3.8956480726576652e-17 | PASS |
| 4 | r99 | 19.25 | 19.375 | 0.0064935064935064939 | 0.00974025974025974 | 0.66666666666666663 | 1 | 0.33317013591547739 | 1.2617854706008803e-18 | 7.1507823999893397e-17 | PASS |
| 4 | radial_profile_L2 | null | null | 0.01390051111786134 | 0.026445474852675486 | 0.52562909894034393 | -7.7751682619068863 | 1.2210780923911863e-06 | 1.7063425300474121e-19 | 5.5915970512964533e-19 | PASS |
| 8 | total_mass | 256 | 255.99987439079101 | 4.9066097267125297e-07 | 0 | null | -1.0352297764902778 | 0.31696390498294125 | 6.331051438669818e-81 | 6.3307851814239781e-81 | PASS |
| 8 | r_K_ratio | 1 | 0.99999901869484664 | 9.8130515333721968e-07 | 0 | null | -1.0352120739000901 | 0.3169719030182524 | 2.0746022011006527e-76 | 2.0744277097679295e-76 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 8 | r90 | 18.0625 | 18 | 0.0034602076124567475 | 0.0034602076124567475 | 1 | -1 | 0.33317013591547739 | 1.4301148022486204e-22 | 1.9701065876166452e-19 | PASS |
| 8 | r99 | 20.0625 | 20 | 0.0031152647975077881 | 0.0031152647975077881 | 1 | -1 | 0.33317013591547739 | 6.5149803197211003e-25 | 5.0766352092358305e-21 | PASS |
| 8 | radial_profile_L2 | null | null | 0.067252396223438593 | 0.031674679515173401 | 2.1232226261743889 | -1.1211565613575198 | 0.2798500299122274 | 1.1233733444798528e-17 | 3.6201109868020346e-15 | PASS |
| 12 | total_mass | 256 | 255.99987438921451 | 4.9066713076612034e-07 | 0 | null | -1.035242772400685 | 0.31695803351978358 | 6.3310511518897591e-81 | 6.3307848913142876e-81 | PASS |
| 12 | r_K_ratio | 1 | 0.99999901869489838 | 9.8130510160082673e-07 | 0 | null | -1.0352120282386914 | 0.31697192364827376 | 2.0746019330439565e-76 | 2.0744274417429762e-76 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 12 | r90 | 18.6875 | 18 | 0.03678929765886288 | 0 | null | -5.7445626465380286 | 3.8761814875216441e-05 | 1.2411344634416193e-17 | 3.0051609001276733e-16 | PASS |
| 12 | r99 | 20.5625 | 20.125 | 0.021276595744680851 | 0.015197568389057751 | 1.3999999999999999 | -3.415650255319866 | 0.0038327022878305922 | 3.2453744019125243e-19 | 3.5829058514479523e-17 | PASS |
| 12 | radial_profile_L2 | null | null | 0.072893691384749507 | 0.033001760468462815 | 2.2087819058746359 | -0.25736842332843862 | 0.80039152143350745 | 4.3685958447085025e-20 | 2.3969087135186269e-17 | PASS |
| 16 | total_mass | 256 | 256.00565975007765 | 2.2108398740866564e-05 | 0 | null | 29.0549503061726 | 1.3336686824631748e-14 | 7.6725742637134177e-78 | 7.6871276199674024e-78 | PASS |
| 16 | r_K_ratio | 1 | 1.0000442168399475 | 4.4216839947416875e-05 | 0 | null | 29.054978697547448 | 1.3336494557056482e-14 | 2.5117676850929536e-73 | 2.5213053845686074e-73 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.0047619047619047623 | 1 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 16 | r90 | 19.0625 | 18 | 0.055737704918032788 | 0.0032786885245901639 | 17 | -17 | 3.2769117377966068e-11 | 2.2623006381054966e-21 | 5.3331248280097439e-21 | PASS |
| 16 | r99 | 21.1875 | 20.8125 | 0.017699115044247787 | 0.014749262536873156 | 1.2 | -2.4227185592617446 | 0.028528068449741654 | 4.6950640203719188e-18 | 3.0093626661365215e-16 | PASS |
| 16 | radial_profile_L2 | null | null | 0.08775376074242433 | 0.029422120612765671 | 2.9825776971476938 | 0.20757506469588111 | 0.83835262262231791 | 5.9173726605820711e-20 | 1.2294617205940524e-16 | PASS |
| 20 | total_mass | 256 | 259.38930481728698 | 0.013239471942527239 | 0 | null | 48.255029981106148 | 7.1535576417554098e-18 | 9.9612836898674566e-40 | 3.1001936901933769e-39 | PASS |
| 20 | r_K_ratio | 1 | 1.0264789234962126 | 0.026478923496212517 | 0 | null | 48.255010551841188 | 7.1536005859322785e-18 | 1.9075543162085774e-35 | 1.8535898820419411e-34 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0 | null | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 20 | r90 | 19.375 | 18.125 | 0.064516129032258063 | 0 | null | -11.180339887498949 | 1.129731174656377e-08 | 6.9886341424657286e-17 | 1.1645945213682297e-17 | PASS |
| 20 | r99 | 21.9375 | 21 | 0.042735042735042736 | 0 | null | -6.5361700549589266 | 9.4202794798815216e-06 | 3.6683425778815099e-19 | 5.1861867803112176e-17 | PASS |
| 20 | radial_profile_L2 | null | null | 0.092421060738694855 | 0.040888641583758337 | 2.2603113519771729 | 1.0369297345499009 | 0.31619654317416779 | 1.1226962030547471e-19 | 3.5490613354533042e-16 | PASS |
| 24 | total_mass | 284.625 | 298.98754644036325 | 0.050461296233160362 | 0.0083443126921387799 | 6.047387974889797 | 19.730067038543012 | 3.8414531453545918e-12 | 1.6655638915222936e-27 | 4.0384969275713948e-22 | PASS |
| 24 | r_K_ratio | 1.2236328125 | 1.3358250453136638 | 0.091687826337742723 | 0.015163607342378291 | 6.0465708632206123 | 19.727446316178806 | 3.8488472047269843e-12 | 2.9675877127982947e-24 | 2.7669535033326392e-17 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.5625 | 13.0625 | 0.036866359447004608 | 0.018433179723502304 | 2 | -3.8729833462074166 | 0.0015017735386984323 | 9.6072681659644597e-15 | 7.6467248765302002e-14 | PASS |
| 24 | r90 | 19.9375 | 18.25 | 0.084639498432601878 | 0.003134796238244514 | 27 | -14.100290132411525 | 4.6333644836412128e-10 | 7.7245411315955695e-15 | 9.481589102890017e-19 | PASS |
| 24 | r99 | 22.875 | 21.25 | 0.071038251366120214 | 0.0054644808743169399 | 13 | -6.7890285822722154 | 6.1053729318667233e-06 | 1.2337518345797697e-14 | 7.416708381213179e-15 | PASS |
| 24 | radial_profile_L2 | null | null | 0.1311571822925125 | 0.031998880089473554 | 4.0988053933693243 | 4.1375477785895995 | 0.00087694673719733397 | 3.3174571066196259e-19 | 3.9701746118309766e-14 | PASS |
| 28 | total_mass | 360.0625 | 355.79428911997218 | 0.011854083332831953 | 0.0034716195105016492 | 3.4145687040222441 | -4.9981373855468707 | 0.00015894216488354686 | 4.9346299105874915e-27 | 1.1189772562667494e-23 | PASS |
| 28 | r_K_ratio | 1.81298828125 | 1.7783820031380952 | 0.019087976723183724 | 0.0053864799353622404 | 3.543682878659058 | -5.1899130998290026 | 0.00010985368152231689 | 4.8401929701778768e-24 | 5.9660093150193663e-21 | PASS |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.875 | 13.0625 | 0.058558558558558557 | 0.022522522522522521 | 2.6000000000000001 | -8.0622577482985509 | 7.8258249630702749e-07 | 2.2778649613890186e-15 | 2.5893641111623486e-16 | PASS |
| 28 | r90 | 20.5625 | 18.75 | 0.088145896656534953 | 0.0060790273556231003 | 14.5 | -11.066875266096559 | 1.2961577489987425e-08 | 1.5680800107416445e-13 | 1.9477056267575937e-16 | PASS |
| 28 | r99 | 23.625 | 21.9375 | 0.071428571428571425 | 0 | null | -8.5098307000225546 | 3.9910679004152207e-07 | 1.2402986697836266e-16 | 5.9545457575235637e-16 | PASS |
| 28 | radial_profile_L2 | null | null | 0.16812366489946665 | 0.025528598153441336 | 6.5856990614583761 | 8.2179631045461328 | 6.1744432441080883e-07 | 2.0668183461795087e-20 | 1.1633368762049021e-13 | PASS |
| 32 | total_mass | 385.5625 | 379.14312167383497 | 0.016649384538602752 | 0.0024315124007132437 | 6.8473368812426925 | -5.396509588004732 | 7.4173402428140978e-05 | 1.4730519212317234e-25 | 2.2336216332001333e-22 | PASS |
| 32 | r_K_ratio | 1.8212005047589639 | 1.9391831914335063 | 0.064782920038866071 | 0.0014984898683687194 | 43.232137504799965 | 14.726328885627707 | 2.5168461149834619e-10 | 2.6740318924724625e-25 | 4.4405753482401013e-18 | PASS |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.022421524663677129 | 2.7999999999999998 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 32 | r90 | 20.8125 | 18.9375 | 0.090090090090090086 | 0.006006006006006006 | 15 | -21.957751641341996 | 8.1207854767429677e-13 | 2.5303126792687863e-18 | 2.9422487852212645e-20 | PASS |
| 32 | r99 | 24.375 | 22 | 0.097435897435897437 | 0.0051282051282051282 | 19 | -15.343884208657718 | 1.4090734244550085e-10 | 8.3953741872725587e-18 | 3.7020237911275972e-18 | PASS |
| 32 | radial_profile_L2 | null | null | 0.17938424840949052 | 0.026253620259367155 | 6.8327433183424358 | 10.818675251643851 | 1.7576228289281234e-08 | 2.8179852830823343e-21 | 5.7756288571017507e-14 | PASS |
| 36 | total_mass | 404.5 | 389.43826184084367 | 0.037235446623377903 | 0.001854140914709518 | 20.082317545541816 | -11.824525086757328 | 5.2854372155451087e-09 | 1.4574793360891194e-25 | 2.4656946513950525e-22 | PASS |
| 36 | r_K_ratio | 1.6452597661040846 | 1.8931472248902017 | 0.15066767199510719 | 0.0025243036803777237 | 59.686824991104899 | 18.653271937262037 | 8.6419882044172036e-12 | 9.373935059396823e-23 | 6.9902926261365534e-12 | PASS |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0089686098654708519 | 7 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 36 | r90 | 21.3125 | 19 | 0.10850439882697947 | 0 | null | -19.322619811082461 | 5.195040620846464e-12 | 7.7864380968352111e-16 | 1.5335258890741899e-18 | PASS |
| 36 | r99 | 24.625 | 22.25 | 0.096446700507614211 | 0.0025380710659898475 | 38 | -11.783299785974805 | 5.5430657455111653e-09 | 3.3901457391067851e-15 | 5.0500172350014425e-17 | PASS |
| 36 | radial_profile_L2 | null | null | 0.19113782511778182 | 0.021200415960106734 | 9.0157582510385588 | 12.937967455874144 | 1.5352977834313781e-09 | 2.3205167683510447e-21 | 1.873559415950491e-13 | PASS |
| 40 | total_mass | 421.125 | 406.62557094179493 | 0.034430226318088571 | 0.0014841199168892847 | 23.199086493128078 | -12.123072892031722 | 3.7595679941645351e-09 | 1.454436824969024e-26 | 8.6049179868181217e-23 | PASS |
| 40 | r_K_ratio | 1.4862818219349081 | 1.6995080743478872 | 0.1434628677187155 | 0.0052088165665525443 | 27.542315204558335 | 16.515053506398107 | 4.9544479197942529e-11 | 1.9779308949787462e-22 | 1.6280039013798181e-11 | PASS |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0044843049327354259 | 14 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 40 | r90 | 21.5 | 19.0625 | 0.11337209302325581 | 0 | null | -19.030051422496395 | 6.4760803247711866e-12 | 6.3605411260372649e-15 | 2.3244786420838592e-18 | PASS |
| 40 | r99 | 25 | 22.6875 | 0.092499999999999999 | 0.0074999999999999997 | 12.333333333333334 | -15.363413772941893 | 1.3839300528290396e-10 | 5.4380328870719136e-17 | 4.8751207713627678e-19 | PASS |
| 40 | radial_profile_L2 | null | null | 0.20294210478022673 | 0.019135549966259355 | 10.605501547541785 | 12.331063058139474 | 2.9774061335968384e-09 | 2.6639476681283848e-22 | 9.5858061539418275e-14 | PASS |
| 44 | total_mass | 433 | 439.12414900465899 | 0.014143531188588849 | 0.0072170900692840644 | 1.9597276814908708 | 4.3124362386774742 | 0.00061619879125071012 | 5.6783101501054603e-26 | 3.3785584158268306e-21 | PASS |
| 44 | r_K_ratio | 1.3851122690599202 | 1.4826505247964381 | 0.070419025168780944 | 0.023883569861476831 | 2.948429634983663 | 10.150700193299443 | 4.1027820504537489e-08 | 1.7459206768544517e-22 | 6.5239376359162711e-15 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 14.0625 | 13.0625 | 0.071111111111111111 | 0.0088888888888888889 | 8 | -7.7459666924148332 | 1.2783477850604916e-06 | 1.25956553381152e-13 | 5.7635238603637946e-15 | PASS |
| 44 | r90 | 21.875 | 19.1875 | 0.12285714285714286 | 0.0057142857142857143 | 21.5 | -22.456017618285021 | 5.8544703070058476e-13 | 4.6040517988789706e-14 | 9.7976516786557512e-20 | PASS |
| 44 | r99 | 25.4375 | 22.9375 | 0.098280098280098274 | 0.0024570024570024569 | 40 | -15.811388300841896 | 9.2073297610083179e-11 | 1.32240574860815e-17 | 2.0753132400270518e-18 | PASS |
| 44 | radial_profile_L2 | null | null | 0.20578903863170117 | 0.022666905945975712 | 9.0788323347782285 | 13.563001661788283 | 7.9734125142345473e-10 | 3.9888272112335183e-23 | 2.1411534627969286e-14 | PASS |
| 48 | total_mass | 468.375 | 495.63025101823672 | 0.058191088376272695 | 0.0066720042700827327 | 8.7216803258357523 | 13.932439671958022 | 5.478965023639955e-10 | 6.2262438154028136e-25 | 9.8983062404198524e-19 | PASS |
| 48 | r_K_ratio | 1.4233832881828432 | 1.4168010867127527 | 0.00462433521928831 | 0.023969344303723984 | 0.19292706386506586 | -0.56934483893448984 | 0.5775492092836918 | 5.7033749883260893e-20 | 1.7002683649562875e-15 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 14.4375 | 13.1875 | 0.086580086580086577 | 0.030303030303030304 | 2.8571428571428572 | -11.180339887498949 | 1.129731174656377e-08 | 5.7520127976418164e-14 | 2.7873182537487998e-16 | PASS |
| 48 | r90 | 22.375 | 19.5 | 0.12849162011173185 | 0 | null | -15.998991903725805 | 7.7865353236989104e-11 | 1.0061554522861548e-11 | 4.8264550809103502e-17 | PASS |
| 48 | r99 | 26.1875 | 23.1875 | 0.11455847255369929 | 0.0071599045346062056 | 16 | -23.2379000772445 | 3.5513304925948273e-13 | 1.7995348202025614e-17 | 8.6169881800138773e-21 | PASS |
| 48 | radial_profile_L2 | null | null | 0.23341720836842744 | 0.017773720507827693 | 13.1327151378142 | 12.748408916387797 | 1.8827621113569947e-09 | 3.1119166670442979e-21 | 6.3693627289909649e-11 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.10546875, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.10546875, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.09375, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

### ablation-regular/correction_report.json

Equivalence by intervention: `{"axial_normal_transport": true, "full": true, "legacy_growth_mean": true, "unconditional_lattice": true}`.

| Intervention | Hours | Metric | ABM effect | PDE effect | Interaction/error change | SE | Paired t | Paired p |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| legacy_growth_mean | 0 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 0 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 0 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 0 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 0 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 0 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 0 | radial_profile_L2 | null | null | 0 | null | null | null |
| legacy_growth_mean | 4 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 4 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 4 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 4 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 4 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 4 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 4 | radial_profile_L2 | null | null | 0 | null | null | null |
| legacy_growth_mean | 8 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 8 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 8 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 8 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 8 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 8 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 8 | radial_profile_L2 | null | null | 0 | null | null | null |
| legacy_growth_mean | 12 | total_mass | 0 | -5.2708060138684232e-11 | -5.2708060138684232e-11 | 1.3475507456281886e-12 | -39.113970519984598 | 1.6310076508290242e-16 |
| legacy_growth_mean | 12 | r_K_ratio | 0 | -4.1175396425785493e-13 | -4.1175396425785493e-13 | 1.0530922850605146e-14 | -39.099513888680136 | 1.6399936781436858e-16 |
| legacy_growth_mean | 12 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 12 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 12 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 12 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 12 | radial_profile_L2 | null | null | 8.4238171993433753e-15 | null | null | null |
| legacy_growth_mean | 16 | total_mass | 0 | -0.0013548832269307809 | -0.0013548832269307809 | 3.2856351343835961e-05 | -41.236569841616472 | 7.4333327120627247e-17 |
| legacy_growth_mean | 16 | r_K_ratio | 0 | -1.0585025768200529e-05 | -1.0585025768200529e-05 | 2.5669025227705708e-07 | -41.23657082535275 | 7.4333300740200638e-17 |
| legacy_growth_mean | 16 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 16 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 16 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 16 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 16 | radial_profile_L2 | null | null | 2.2873373650700302e-07 | null | null | null |
| legacy_growth_mean | 20 | total_mass | 0 | -0.51222337681729257 | -0.51222337681729257 | 0.010139020187334471 | -50.520007589801942 | 3.6093610640102084e-18 |
| legacy_growth_mean | 20 | r_K_ratio | 0 | -0.0040017611655629587 | -0.0040017611655629587 | 7.9211244743189593e-05 | -50.520114644543838 | 3.6092469702563372e-18 |
| legacy_growth_mean | 20 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 20 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 20 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 20 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 20 | radial_profile_L2 | null | null | 8.2508845083789639e-05 | null | null | null |
| legacy_growth_mean | 24 | total_mass | 0 | -3.3144358904289248 | -3.3144358904289248 | 0.036527247009110725 | -90.738726890702424 | 5.6869027420350231e-22 |
| legacy_growth_mean | 24 | r_K_ratio | 0 | -0.025903330024911034 | -0.025903330024911034 | 0.00028544052385554074 | -90.748607363194552 | 5.6776380588764706e-22 |
| legacy_growth_mean | 24 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 24 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 24 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 24 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 24 | radial_profile_L2 | null | null | 0.00031920650579114751 | null | null | null |
| legacy_growth_mean | 28 | total_mass | 0 | -2.6964674948193839 | -2.6964674948193839 | 0.020716804694583892 | -130.1584648101807 | 2.5559159288847309e-24 |
| legacy_growth_mean | 28 | r_K_ratio | 0 | -0.02167481050424698 | -0.02167481050424698 | 0.00015644300046953769 | -138.54765275016226 | 1.0022215487951639e-24 |
| legacy_growth_mean | 28 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 28 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 28 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 28 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 28 | radial_profile_L2 | null | null | 0.00015679904073243045 | null | null | null |
| legacy_growth_mean | 32 | total_mass | 0 | -0.42159100973291785 | -0.42159100973291785 | 0.033716067087915302 | -12.504157398714716 | 2.4581668771253039e-09 |
| legacy_growth_mean | 32 | r_K_ratio | 0 | -0.011759123973643607 | -0.011759123973643607 | 0.00024131215214583967 | -48.729928721272394 | 6.1816570341381747e-18 |
| legacy_growth_mean | 32 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 32 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 32 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 32 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 32 | radial_profile_L2 | null | null | 0.00010152356398987483 | null | null | null |
| legacy_growth_mean | 36 | total_mass | 0 | 1.6441626322730265 | 1.6441626322730265 | 0.031186346465206911 | 52.720591496901974 | 1.9104769476611809e-18 |
| legacy_growth_mean | 36 | r_K_ratio | 0 | -0.026433134569651212 | -0.026433134569651212 | 0.00033856360497263709 | -78.074353478684017 | 5.3974100218055031e-21 |
| legacy_growth_mean | 36 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 36 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 36 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 36 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 36 | radial_profile_L2 | null | null | 0.00022344298610388336 | null | null | null |
| legacy_growth_mean | 40 | total_mass | 0 | 4.1976020777141088 | 4.1976020777141088 | 0.048622320810658731 | 86.330763479186544 | 1.1986669298001418e-21 |
| legacy_growth_mean | 40 | r_K_ratio | 0 | -0.048358515863126075 | -0.048358515863126075 | 0.00032257347122876206 | -149.91473315805089 | 3.0735185594357952e-25 |
| legacy_growth_mean | 40 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 40 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 40 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 40 | r99 | 0 | -0.0625 | -0.0625 | 0.0625 | -1 | 0.33317013591547739 |
| legacy_growth_mean | 40 | radial_profile_L2 | null | null | 0.00049833821184336324 | null | null | null |
| legacy_growth_mean | 44 | total_mass | 0 | 4.5447654442418859 | 4.5447654442418859 | 0.047902951921343292 | 94.874433870054531 | 2.917450808057221e-22 |
| legacy_growth_mean | 44 | r_K_ratio | 0 | -0.058479572827443968 | -0.058479572827443968 | 0.0002427955711206686 | -240.85930627778961 | 2.5124761744591586e-28 |
| legacy_growth_mean | 44 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 44 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 44 | r90 | 0 | -0.0625 | -0.0625 | 0.0625 | -1 | 0.33317013591547739 |
| legacy_growth_mean | 44 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 44 | radial_profile_L2 | null | null | 0.0010484219444226173 | null | null | null |
| legacy_growth_mean | 48 | total_mass | 0 | -0.18772253662969263 | -0.18772253662969263 | 0.073642587113581515 | -2.5491029577785156 | 0.022237479317375783 |
| legacy_growth_mean | 48 | r_K_ratio | 0 | -0.069097011925742149 | -0.069097011925742149 | 0.00051639270403372808 | -133.80710336532775 | 1.6888503728577956e-24 |
| legacy_growth_mean | 48 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 48 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 48 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 48 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| legacy_growth_mean | 48 | radial_profile_L2 | null | null | 0.0014993497897471675 | null | null | null |
| axial_normal_transport | 0 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 0 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 0 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 0 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 0 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 0 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 0 | radial_profile_L2 | null | null | 0 | null | null | null |
| axial_normal_transport | 4 | total_mass | 0 | -0.00035076012503942877 | -0.00035076012503942877 | 0.00033153098305723145 | -1.0580010405207825 | 0.30679646898755764 |
| axial_normal_transport | 4 | r_K_ratio | 0 | -2.7403151370633538e-06 | -2.7403151370633538e-06 | 2.5900858049319083e-06 | -1.0580016815834388 | 0.30679618613906756 |
| axial_normal_transport | 4 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 4 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 4 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 4 | r99 | 0 | -0.125 | -0.125 | 0.085391256382996647 | -1.4638501094227998 | 0.16387561365565406 |
| axial_normal_transport | 4 | radial_profile_L2 | null | null | 0.0041384544810808173 | null | null | null |
| axial_normal_transport | 8 | total_mass | 0 | -0.00035075973581299991 | -0.00035075973581299991 | 0.00033153098256511452 | -1.0579998680639411 | 0.30679698629730845 |
| axial_normal_transport | 8 | r_K_ratio | 0 | -2.7403151192442743e-06 | -2.7403151192442743e-06 | 2.5900858176102967e-06 | -1.0580016695248284 | 0.3067961914595434 |
| axial_normal_transport | 8 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 8 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 8 | r90 | 0 | -0.0625 | -0.0625 | 0.0625 | -1 | 0.33317013591547739 |
| axial_normal_transport | 8 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 8 | radial_profile_L2 | null | null | 0.0090948955471488041 | null | null | null |
| axial_normal_transport | 12 | total_mass | 0 | -0.00035075934268036235 | -0.00035075934268036235 | 0.00033153098277032792 | -1.057998681599406 | 0.30679750978816245 |
| axial_normal_transport | 12 | r_K_ratio | 0 | -2.7403145138188423e-06 | -2.7403145138188423e-06 | 2.590085817996387e-06 | -1.0580014356198699 | 0.30679629466263969 |
| axial_normal_transport | 12 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 12 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 12 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 12 | r99 | 0 | -0.125 | -0.125 | 0.085391256382996647 | -1.4638501094227998 | 0.16387561365565406 |
| axial_normal_transport | 12 | radial_profile_L2 | null | null | 0.013119959471610682 | null | null | null |
| axial_normal_transport | 16 | total_mass | 0 | -0.00039219681264412998 | -0.00039219681264412998 | 0.00033092801907639132 | -1.1851423573583697 | 0.25439828708500056 |
| axial_normal_transport | 16 | r_K_ratio | 0 | -3.0640482741178809e-06 | -3.0640482741178809e-06 | 2.5853751598559205e-06 | -1.1851464815221 | 0.2543967056370986 |
| axial_normal_transport | 16 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 16 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 16 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 16 | r99 | 0 | -0.5 | -0.5 | 0.12909944487358058 | -3.8729833462074166 | 0.0015017735386984323 |
| axial_normal_transport | 16 | radial_profile_L2 | null | null | 0.017071354155039867 | null | null | null |
| axial_normal_transport | 20 | total_mass | 0 | -0.018494736910870557 | -0.018494736910870557 | 0.00073574725814329616 | -25.137350776601156 | 1.1236361349340641e-13 |
| axial_normal_transport | 20 | r_K_ratio | 0 | -0.0001444899968784108 | -0.0001444899968784108 | 5.7480247732708369e-06 | -25.137330226952152 | 1.1236496127717266e-13 |
| axial_normal_transport | 20 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 20 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 20 | r90 | 0 | -0.125 | -0.125 | 0.085391256382996647 | -1.4638501094227998 | 0.16387561365565406 |
| axial_normal_transport | 20 | r99 | 0 | -0.1875 | -0.1875 | 0.10077822185373186 | -1.8605210188381269 | 0.082530706346961857 |
| axial_normal_transport | 20 | radial_profile_L2 | null | null | 0.019644369774697176 | null | null | null |
| axial_normal_transport | 24 | total_mass | 0 | -0.13563347423269434 | -0.13563347423269434 | 0.0044377170273770788 | -30.563795166737066 | 6.3153805129390129e-15 |
| axial_normal_transport | 24 | r_K_ratio | 0 | -0.001059517578861191 | -0.001059517578861191 | 3.4669307213585628e-05 | -30.560679287124749 | 6.3249019814719375e-15 |
| axial_normal_transport | 24 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 24 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 24 | r90 | 0 | -0.1875 | -0.1875 | 0.10077822185373186 | -1.8605210188381269 | 0.082530706346961857 |
| axial_normal_transport | 24 | r99 | 0 | -0.25 | -0.25 | 0.11180339887498948 | -2.2360679774997898 | 0.040968955955836141 |
| axial_normal_transport | 24 | radial_profile_L2 | null | null | 0.023596225171872126 | null | null | null |
| axial_normal_transport | 28 | total_mass | 0 | -0.26582728880722328 | -0.26582728880722328 | 0.04571197779970948 | -5.81526553000979 | 3.404204010011615e-05 |
| axial_normal_transport | 28 | r_K_ratio | 0 | -0.002066484123058307 | -0.002066484123058307 | 0.00035686091973195045 | -5.7907268876920135 | 3.560842276931662e-05 |
| axial_normal_transport | 28 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 28 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 28 | r90 | 0 | -0.5625 | -0.5625 | 0.12808688457449499 | -4.3915503282683988 | 0.00052573078993947925 |
| axial_normal_transport | 28 | r99 | 0 | -0.625 | -0.625 | 0.125 | -5 | 0.0001583695146220272 |
| axial_normal_transport | 28 | radial_profile_L2 | null | null | 0.028232846580911436 | null | null | null |
| axial_normal_transport | 32 | total_mass | 0 | -0.51433679146442657 | -0.51433679146442657 | 0.090921937546508458 | -5.6569053117827872 | 4.5574618568156902e-05 |
| axial_normal_transport | 32 | r_K_ratio | 0 | -0.003797101436177619 | -0.003797101436177619 | 0.00069801599780435549 | -5.4398487257048451 | 6.8355316068437243e-05 |
| axial_normal_transport | 32 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 32 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 32 | r90 | 0 | -0.5625 | -0.5625 | 0.12808688457449499 | -4.3915503282683988 | 0.00052573078993947925 |
| axial_normal_transport | 32 | r99 | 0 | -0.1875 | -0.1875 | 0.10077822185373186 | -1.8605210188381269 | 0.082530706346961857 |
| axial_normal_transport | 32 | radial_profile_L2 | null | null | 0.030792673325043152 | null | null | null |
| axial_normal_transport | 36 | total_mass | 0 | -0.56186039433823609 | -0.56186039433823609 | 0.087736068745032939 | -6.403984158112225 | 1.1858882720198623e-05 |
| axial_normal_transport | 36 | r_K_ratio | 0 | -0.0029081854352815406 | -0.0029081854352815406 | 0.0005801819006442735 | -5.0125407773873905 | 0.00015456895210970596 |
| axial_normal_transport | 36 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 36 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 36 | r90 | 0 | -0.375 | -0.375 | 0.125 | -3 | 0.0089727374772233335 |
| axial_normal_transport | 36 | r99 | 0 | -0.25 | -0.25 | 0.11180339887498948 | -2.2360679774997898 | 0.040968955955836141 |
| axial_normal_transport | 36 | radial_profile_L2 | null | null | 0.034065816151554767 | null | null | null |
| axial_normal_transport | 40 | total_mass | 0 | -0.752958249057464 | -0.752958249057464 | 0.10188973941803413 | -7.3899320320000053 | 2.2551350870716943e-06 |
| axial_normal_transport | 40 | r_K_ratio | 0 | -0.0013850207764742883 | -0.0013850207764742883 | 0.00044946046487713537 | -3.0815185866300787 | 0.0075986728714294623 |
| axial_normal_transport | 40 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 40 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 40 | r90 | 0 | -0.25 | -0.25 | 0.11180339887498948 | -2.2360679774997898 | 0.040968955955836141 |
| axial_normal_transport | 40 | r99 | 0 | -0.6875 | -0.6875 | 0.11967838846954226 | -5.7445626465380286 | 3.8761814875216441e-05 |
| axial_normal_transport | 40 | radial_profile_L2 | null | null | 0.036057848160039041 | null | null | null |
| axial_normal_transport | 44 | total_mass | 0 | -2.0296856718791219 | -2.0296856718791219 | 0.17449065893544929 | -11.632059184497546 | 6.6082522441796848e-09 |
| axial_normal_transport | 44 | r_K_ratio | 0 | -0.0053517658392852052 | -0.0053517658392852052 | 0.0005731968981010074 | -9.336697140224457 | 1.223123230424243e-07 |
| axial_normal_transport | 44 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 44 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 44 | r90 | 0 | -0.25 | -0.25 | 0.11180339887498948 | -2.2360679774997898 | 0.040968955955836141 |
| axial_normal_transport | 44 | r99 | 0 | -0.5625 | -0.5625 | 0.12808688457449499 | -4.3915503282683988 | 0.00052573078993947925 |
| axial_normal_transport | 44 | radial_profile_L2 | null | null | 0.037838711507482053 | null | null | null |
| axial_normal_transport | 48 | total_mass | 0 | -5.9689576632359014 | -5.9689576632359014 | 0.21705181668952589 | -27.500150674960654 | 2.9998985472537727e-14 |
| axial_normal_transport | 48 | r_K_ratio | 0 | -0.021005036176118402 | -0.021005036176118402 | 0.0008680565316125968 | -24.197774466482212 | 1.964100434721749e-13 |
| axial_normal_transport | 48 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 48 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| axial_normal_transport | 48 | r90 | 0 | -0.5 | -0.5 | 0.12909944487358058 | -3.8729833462074166 | 0.0015017735386984323 |
| axial_normal_transport | 48 | r99 | 0 | -0.25 | -0.25 | 0.11180339887498948 | -2.2360679774997898 | 0.040968955955836141 |
| axial_normal_transport | 48 | radial_profile_L2 | null | null | 0.03983809105932945 | null | null | null |
| unconditional_lattice | 0 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 0 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 0 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 0 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 0 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 0 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 0 | radial_profile_L2 | null | null | 0 | null | null | null |
| unconditional_lattice | 4 | total_mass | 0 | -7.3649122089847197e-05 | -7.3649122089847197e-05 | 6.9757823737571226e-05 | -1.0557829666090937 | 0.30777626383338164 |
| unconditional_lattice | 4 | r_K_ratio | 0 | -5.7538366293741205e-07 | -5.7538366293741205e-07 | 5.4498298777884508e-07 | -1.0557827966015403 | 0.30777633901860568 |
| unconditional_lattice | 4 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 4 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 4 | r90 | 0 | -0.125 | -0.125 | 0.085391256382996647 | -1.4638501094227998 | 0.16387561365565406 |
| unconditional_lattice | 4 | r99 | 0 | -0.1875 | -0.1875 | 0.10077822185373186 | -1.8605210188381269 | 0.082530706346961857 |
| unconditional_lattice | 4 | radial_profile_L2 | null | null | 0.0031394465752887235 | null | null | null |
| unconditional_lattice | 8 | total_mass | 0 | -7.364913758323155e-05 | -7.364913758323155e-05 | 6.9757823707152145e-05 | -1.0557831891719471 | 0.30777616540578667 |
| unconditional_lattice | 8 | r_K_ratio | 0 | -5.7538339069684863e-07 | -5.7538339069684863e-07 | 5.4498299589257252e-07 | -1.0557822813434508 | 0.30777656688970845 |
| unconditional_lattice | 8 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 8 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 8 | r90 | 0 | -0.0625 | -0.0625 | 0.0625 | -1 | 0.33317013591547739 |
| unconditional_lattice | 8 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 8 | radial_profile_L2 | null | null | 0.0080946165967747535 | null | null | null |
| unconditional_lattice | 12 | total_mass | 0 | -7.3649132346531587e-05 | -7.3649132346531587e-05 | 6.9757824649919534e-05 | -1.0557830998334685 | 0.30777620491539481 |
| unconditional_lattice | 12 | r_K_ratio | 0 | -5.7538359272274464e-07 | -5.7538359272274464e-07 | 5.4498299671627405e-07 | -1.0557826504489967 | 0.30777640365404857 |
| unconditional_lattice | 12 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 12 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 12 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 12 | r99 | 0 | -0.125 | -0.125 | 0.085391256382996647 | -1.4638501094227998 | 0.16387561365565406 |
| unconditional_lattice | 12 | radial_profile_L2 | null | null | 0.0124759786334738 | null | null | null |
| unconditional_lattice | 16 | total_mass | 0 | -0.00012408566427879464 | -0.00012408566427879464 | 6.9014728206546279e-05 | -1.7979591820956371 | 0.092337656853166555 |
| unconditional_lattice | 16 | r_K_ratio | 0 | -9.6941891775115252e-07 | -9.6941891775115252e-07 | 5.3917755189007974e-07 | -1.7979586026029226 | 0.092337752283329827 |
| unconditional_lattice | 16 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 16 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 16 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 16 | r99 | 0 | -0.5 | -0.5 | 0.12909944487358058 | -3.8729833462074166 | 0.0015017735386984323 |
| unconditional_lattice | 16 | radial_profile_L2 | null | null | 0.016548448976341126 | null | null | null |
| unconditional_lattice | 20 | total_mass | 0 | -0.020773233234738342 | -0.020773233234738342 | 0.00080241197201444577 | -25.888488655755452 | 7.2926645996163482e-14 |
| unconditional_lattice | 20 | r_K_ratio | 0 | -0.00016229065317086011 | -0.00016229065317086011 | 6.2688424310752575e-06 | -25.88845627485701 | 7.2927986047190727e-14 |
| unconditional_lattice | 20 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 20 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 20 | r90 | 0 | -0.125 | -0.125 | 0.085391256382996647 | -1.4638501094227998 | 0.16387561365565406 |
| unconditional_lattice | 20 | r99 | 0 | -0.125 | -0.125 | 0.085391256382996647 | -1.4638501094227998 | 0.16387561365565406 |
| unconditional_lattice | 20 | radial_profile_L2 | null | null | 0.019306855116084293 | null | null | null |
| unconditional_lattice | 24 | total_mass | 0 | -0.148108547517662 | -0.148108547517662 | 0.0046407629662026939 | -31.914697776269293 | 3.3316959022739196e-15 |
| unconditional_lattice | 24 | r_K_ratio | 0 | -0.0011569331288104334 | -0.0011569331288104334 | 3.6255086700029898e-05 | -31.910918828653067 | 3.3375379972149031e-15 |
| unconditional_lattice | 24 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 24 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 24 | r90 | 0 | -0.1875 | -0.1875 | 0.10077822185373186 | -1.8605210188381269 | 0.082530706346961857 |
| unconditional_lattice | 24 | r99 | 0 | -0.25 | -0.25 | 0.11180339887498948 | -2.2360679774997898 | 0.040968955955836141 |
| unconditional_lattice | 24 | radial_profile_L2 | null | null | 0.023403317530079737 | null | null | null |
| unconditional_lattice | 28 | total_mass | 0 | -0.27627176812544008 | -0.27627176812544008 | 0.046895535117581585 | -5.8912168809406156 | 2.9632660885928902e-05 |
| unconditional_lattice | 28 | r_K_ratio | 0 | -0.0021452148152148048 | -0.0021452148152148048 | 0.00036609149132009873 | -5.8597778590245833 | 3.1380996508865633e-05 |
| unconditional_lattice | 28 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 28 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 28 | r90 | 0 | -0.5625 | -0.5625 | 0.12808688457449499 | -4.3915503282683988 | 0.00052573078993947925 |
| unconditional_lattice | 28 | r99 | 0 | -0.5625 | -0.5625 | 0.12808688457449499 | -4.3915503282683988 | 0.00052573078993947925 |
| unconditional_lattice | 28 | radial_profile_L2 | null | null | 0.028148143228199723 | null | null | null |
| unconditional_lattice | 32 | total_mass | 0 | -0.51910037703153478 | -0.51910037703153478 | 0.091396398424900549 | -5.6796590016407924 | 4.3694488240740845e-05 |
| unconditional_lattice | 32 | r_K_ratio | 0 | -0.0037960087885948873 | -0.0037960087885948873 | 0.00070138912133355146 | -5.41212954854152 | 7.2019469621084509e-05 |
| unconditional_lattice | 32 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 32 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 32 | r90 | 0 | -0.625 | -0.625 | 0.125 | -5 | 0.0001583695146220272 |
| unconditional_lattice | 32 | r99 | 0 | -0.1875 | -0.1875 | 0.10077822185373186 | -1.8605210188381269 | 0.082530706346961857 |
| unconditional_lattice | 32 | radial_profile_L2 | null | null | 0.030796706655422612 | null | null | null |
| unconditional_lattice | 36 | total_mass | 0 | -0.57000568013064878 | -0.57000568013064878 | 0.087409396836365338 | -6.5211030022061278 | 9.6695706739500783e-06 |
| unconditional_lattice | 36 | r_K_ratio | 0 | -0.0027963734513067745 | -0.0027963734513067745 | 0.00057567085366532259 | -4.8575908151370477 | 0.00020894781076550773 |
| unconditional_lattice | 36 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 36 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 36 | r90 | 0 | -0.5 | -0.5 | 0.12909944487358058 | -3.8729833462074166 | 0.0015017735386984323 |
| unconditional_lattice | 36 | r99 | 0 | -0.25 | -0.25 | 0.11180339887498948 | -2.2360679774997898 | 0.040968955955836141 |
| unconditional_lattice | 36 | radial_profile_L2 | null | null | 0.034145538645972334 | null | null | null |
| unconditional_lattice | 40 | total_mass | 0 | -0.77649524179394191 | -0.77649524179394191 | 0.10165464054568971 | -7.638561679286429 | 1.5145008244465357e-06 |
| unconditional_lattice | 40 | r_K_ratio | 0 | -0.0011959233571199601 | -0.0011959233571199601 | 0.00044279615085408234 | -2.7008440674410039 | 0.016430951372654389 |
| unconditional_lattice | 40 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 40 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 40 | r90 | 0 | -0.25 | -0.25 | 0.11180339887498948 | -2.2360679774997898 | 0.040968955955836141 |
| unconditional_lattice | 40 | r99 | 0 | -0.6875 | -0.6875 | 0.11967838846954226 | -5.7445626465380286 | 3.8761814875216441e-05 |
| unconditional_lattice | 40 | radial_profile_L2 | null | null | 0.036200238705557636 | null | null | null |
| unconditional_lattice | 44 | total_mass | 0 | -2.0708397336941395 | -2.0708397336941395 | 0.17574872303084202 | -11.782957497396605 | 5.5452596998036987e-09 |
| unconditional_lattice | 44 | r_K_ratio | 0 | -0.005181250044684102 | -0.005181250044684102 | 0.00057578205031847535 | -8.9986307176791289 | 1.9652926840326025e-07 |
| unconditional_lattice | 44 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 44 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 44 | r90 | 0 | -0.25 | -0.25 | 0.11180339887498948 | -2.2360679774997898 | 0.040968955955836141 |
| unconditional_lattice | 44 | r99 | 0 | -0.5625 | -0.5625 | 0.12808688457449499 | -4.3915503282683988 | 0.00052573078993947925 |
| unconditional_lattice | 44 | radial_profile_L2 | null | null | 0.038035325088707106 | null | null | null |
| unconditional_lattice | 48 | total_mass | 0 | -6.043381428239968 | -6.043381428239968 | 0.23084121131155894 | -26.179820292501461 | 6.1868934000244605e-14 |
| unconditional_lattice | 48 | r_K_ratio | 0 | -0.021006566001453605 | -0.021006566001453605 | 0.00093002230110981424 | -22.587163744768322 | 5.3774814417931102e-13 |
| unconditional_lattice | 48 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 48 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| unconditional_lattice | 48 | r90 | 0 | -0.5 | -0.5 | 0.12909944487358058 | -3.8729833462074166 | 0.0015017735386984323 |
| unconditional_lattice | 48 | r99 | 0 | -0.25 | -0.25 | 0.11180339887498948 | -2.2360679774997898 | 0.040968955955836141 |
| unconditional_lattice | 48 | radial_profile_L2 | null | null | 0.040061878843109411 | null | null | null |

#### full: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 255.99987439227883 | 4.9065516088270256e-07 | 0 | null | -1.0352175095070719 | 0.31696944719393905 | 6.3310518845588904e-81 | 6.3307856304481416e-81 | PASS |
| 4 | r_K_ratio | 1 | 0.9999990186955926 | 9.8130440741306391e-07 | 0 | null | -1.0352112711393178 | 0.31697226570925685 | 2.0746026777977294e-76 | 2.0744281865575098e-76 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13.0625 | 13.0625 | 0 | 0.0047846889952153108 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 4 | r90 | 17.9375 | 17.6875 | 0.013937282229965157 | 0.0069686411149825784 | 2 | -2.2360679774997898 | 0.040968955955836141 | 1.1148111900685798e-16 | 3.8956480726576652e-17 | PASS |
| 4 | r99 | 19.25 | 19.375 | 0.0064935064935064939 | 0.00974025974025974 | 0.66666666666666663 | 1 | 0.33317013591547739 | 1.2617854706008803e-18 | 7.1507823999893397e-17 | PASS |
| 4 | radial_profile_L2 | null | null | 0.01390051111786134 | 0.026445474852675486 | 0.52562909894034393 | -7.7751682619068863 | 1.2210780923911863e-06 | 1.7063425300474121e-19 | 5.5915970512964533e-19 | PASS |
| 8 | total_mass | 256 | 255.99987439079101 | 4.9066097267125297e-07 | 0 | null | -1.0352297764902778 | 0.31696390498294125 | 6.331051438669818e-81 | 6.3307851814239781e-81 | PASS |
| 8 | r_K_ratio | 1 | 0.99999901869484664 | 9.8130515333721968e-07 | 0 | null | -1.0352120739000901 | 0.3169719030182524 | 2.0746022011006527e-76 | 2.0744277097679295e-76 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 8 | r90 | 18.0625 | 18 | 0.0034602076124567475 | 0.0034602076124567475 | 1 | -1 | 0.33317013591547739 | 1.4301148022486204e-22 | 1.9701065876166452e-19 | PASS |
| 8 | r99 | 20.0625 | 20 | 0.0031152647975077881 | 0.0031152647975077881 | 1 | -1 | 0.33317013591547739 | 6.5149803197211003e-25 | 5.0766352092358305e-21 | PASS |
| 8 | radial_profile_L2 | null | null | 0.067252396223438593 | 0.031674679515173401 | 2.1232226261743889 | -1.1211565613575198 | 0.2798500299122274 | 1.1233733444798528e-17 | 3.6201109868020346e-15 | PASS |
| 12 | total_mass | 256 | 255.99987438921451 | 4.9066713076612034e-07 | 0 | null | -1.035242772400685 | 0.31695803351978358 | 6.3310511518897591e-81 | 6.3307848913142876e-81 | PASS |
| 12 | r_K_ratio | 1 | 0.99999901869489838 | 9.8130510160082673e-07 | 0 | null | -1.0352120282386914 | 0.31697192364827376 | 2.0746019330439565e-76 | 2.0744274417429762e-76 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 12 | r90 | 18.6875 | 18 | 0.03678929765886288 | 0 | null | -5.7445626465380286 | 3.8761814875216441e-05 | 1.2411344634416193e-17 | 3.0051609001276733e-16 | PASS |
| 12 | r99 | 20.5625 | 20.125 | 0.021276595744680851 | 0.015197568389057751 | 1.3999999999999999 | -3.415650255319866 | 0.0038327022878305922 | 3.2453744019125243e-19 | 3.5829058514479523e-17 | PASS |
| 12 | radial_profile_L2 | null | null | 0.072893691384749507 | 0.033001760468462815 | 2.2087819058746359 | -0.25736842332843862 | 0.80039152143350745 | 4.3685958447085025e-20 | 2.3969087135186269e-17 | PASS |
| 16 | total_mass | 256 | 256.00565975007765 | 2.2108398740866564e-05 | 0 | null | 29.0549503061726 | 1.3336686824631748e-14 | 7.6725742637134177e-78 | 7.6871276199674024e-78 | PASS |
| 16 | r_K_ratio | 1 | 1.0000442168399475 | 4.4216839947416875e-05 | 0 | null | 29.054978697547448 | 1.3336494557056482e-14 | 2.5117676850929536e-73 | 2.5213053845686074e-73 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.0047619047619047623 | 1 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 16 | r90 | 19.0625 | 18 | 0.055737704918032788 | 0.0032786885245901639 | 17 | -17 | 3.2769117377966068e-11 | 2.2623006381054966e-21 | 5.3331248280097439e-21 | PASS |
| 16 | r99 | 21.1875 | 20.8125 | 0.017699115044247787 | 0.014749262536873156 | 1.2 | -2.4227185592617446 | 0.028528068449741654 | 4.6950640203719188e-18 | 3.0093626661365215e-16 | PASS |
| 16 | radial_profile_L2 | null | null | 0.08775376074242433 | 0.029422120612765671 | 2.9825776971476938 | 0.20757506469588111 | 0.83835262262231791 | 5.9173726605820711e-20 | 1.2294617205940524e-16 | PASS |
| 20 | total_mass | 256 | 259.38930481728698 | 0.013239471942527239 | 0 | null | 48.255029981106148 | 7.1535576417554098e-18 | 9.9612836898674566e-40 | 3.1001936901933769e-39 | PASS |
| 20 | r_K_ratio | 1 | 1.0264789234962126 | 0.026478923496212517 | 0 | null | 48.255010551841188 | 7.1536005859322785e-18 | 1.9075543162085774e-35 | 1.8535898820419411e-34 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0 | null | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 20 | r90 | 19.375 | 18.125 | 0.064516129032258063 | 0 | null | -11.180339887498949 | 1.129731174656377e-08 | 6.9886341424657286e-17 | 1.1645945213682297e-17 | PASS |
| 20 | r99 | 21.9375 | 21 | 0.042735042735042736 | 0 | null | -6.5361700549589266 | 9.4202794798815216e-06 | 3.6683425778815099e-19 | 5.1861867803112176e-17 | PASS |
| 20 | radial_profile_L2 | null | null | 0.092421060738694855 | 0.040888641583758337 | 2.2603113519771729 | 1.0369297345499009 | 0.31619654317416779 | 1.1226962030547471e-19 | 3.5490613354533042e-16 | PASS |
| 24 | total_mass | 284.625 | 298.98754644036325 | 0.050461296233160362 | 0.0083443126921387799 | 6.047387974889797 | 19.730067038543012 | 3.8414531453545918e-12 | 1.6655638915222936e-27 | 4.0384969275713948e-22 | PASS |
| 24 | r_K_ratio | 1.2236328125 | 1.3358250453136638 | 0.091687826337742723 | 0.015163607342378291 | 6.0465708632206123 | 19.727446316178806 | 3.8488472047269843e-12 | 2.9675877127982947e-24 | 2.7669535033326392e-17 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.5625 | 13.0625 | 0.036866359447004608 | 0.018433179723502304 | 2 | -3.8729833462074166 | 0.0015017735386984323 | 9.6072681659644597e-15 | 7.6467248765302002e-14 | PASS |
| 24 | r90 | 19.9375 | 18.25 | 0.084639498432601878 | 0.003134796238244514 | 27 | -14.100290132411525 | 4.6333644836412128e-10 | 7.7245411315955695e-15 | 9.481589102890017e-19 | PASS |
| 24 | r99 | 22.875 | 21.25 | 0.071038251366120214 | 0.0054644808743169399 | 13 | -6.7890285822722154 | 6.1053729318667233e-06 | 1.2337518345797697e-14 | 7.416708381213179e-15 | PASS |
| 24 | radial_profile_L2 | null | null | 0.1311571822925125 | 0.031998880089473554 | 4.0988053933693243 | 4.1375477785895995 | 0.00087694673719733397 | 3.3174571066196259e-19 | 3.9701746118309766e-14 | PASS |
| 28 | total_mass | 360.0625 | 355.79428911997218 | 0.011854083332831953 | 0.0034716195105016492 | 3.4145687040222441 | -4.9981373855468707 | 0.00015894216488354686 | 4.9346299105874915e-27 | 1.1189772562667494e-23 | PASS |
| 28 | r_K_ratio | 1.81298828125 | 1.7783820031380952 | 0.019087976723183724 | 0.0053864799353622404 | 3.543682878659058 | -5.1899130998290026 | 0.00010985368152231689 | 4.8401929701778768e-24 | 5.9660093150193663e-21 | PASS |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.875 | 13.0625 | 0.058558558558558557 | 0.022522522522522521 | 2.6000000000000001 | -8.0622577482985509 | 7.8258249630702749e-07 | 2.2778649613890186e-15 | 2.5893641111623486e-16 | PASS |
| 28 | r90 | 20.5625 | 18.75 | 0.088145896656534953 | 0.0060790273556231003 | 14.5 | -11.066875266096559 | 1.2961577489987425e-08 | 1.5680800107416445e-13 | 1.9477056267575937e-16 | PASS |
| 28 | r99 | 23.625 | 21.9375 | 0.071428571428571425 | 0 | null | -8.5098307000225546 | 3.9910679004152207e-07 | 1.2402986697836266e-16 | 5.9545457575235637e-16 | PASS |
| 28 | radial_profile_L2 | null | null | 0.16812366489946665 | 0.025528598153441336 | 6.5856990614583761 | 8.2179631045461328 | 6.1744432441080883e-07 | 2.0668183461795087e-20 | 1.1633368762049021e-13 | PASS |
| 32 | total_mass | 385.5625 | 379.14312167383497 | 0.016649384538602752 | 0.0024315124007132437 | 6.8473368812426925 | -5.396509588004732 | 7.4173402428140978e-05 | 1.4730519212317234e-25 | 2.2336216332001333e-22 | PASS |
| 32 | r_K_ratio | 1.8212005047589639 | 1.9391831914335063 | 0.064782920038866071 | 0.0014984898683687194 | 43.232137504799965 | 14.726328885627707 | 2.5168461149834619e-10 | 2.6740318924724625e-25 | 4.4405753482401013e-18 | PASS |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.022421524663677129 | 2.7999999999999998 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 32 | r90 | 20.8125 | 18.9375 | 0.090090090090090086 | 0.006006006006006006 | 15 | -21.957751641341996 | 8.1207854767429677e-13 | 2.5303126792687863e-18 | 2.9422487852212645e-20 | PASS |
| 32 | r99 | 24.375 | 22 | 0.097435897435897437 | 0.0051282051282051282 | 19 | -15.343884208657718 | 1.4090734244550085e-10 | 8.3953741872725587e-18 | 3.7020237911275972e-18 | PASS |
| 32 | radial_profile_L2 | null | null | 0.17938424840949052 | 0.026253620259367155 | 6.8327433183424358 | 10.818675251643851 | 1.7576228289281234e-08 | 2.8179852830823343e-21 | 5.7756288571017507e-14 | PASS |
| 36 | total_mass | 404.5 | 389.43826184084367 | 0.037235446623377903 | 0.001854140914709518 | 20.082317545541816 | -11.824525086757328 | 5.2854372155451087e-09 | 1.4574793360891194e-25 | 2.4656946513950525e-22 | PASS |
| 36 | r_K_ratio | 1.6452597661040846 | 1.8931472248902017 | 0.15066767199510719 | 0.0025243036803777237 | 59.686824991104899 | 18.653271937262037 | 8.6419882044172036e-12 | 9.373935059396823e-23 | 6.9902926261365534e-12 | PASS |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0089686098654708519 | 7 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 36 | r90 | 21.3125 | 19 | 0.10850439882697947 | 0 | null | -19.322619811082461 | 5.195040620846464e-12 | 7.7864380968352111e-16 | 1.5335258890741899e-18 | PASS |
| 36 | r99 | 24.625 | 22.25 | 0.096446700507614211 | 0.0025380710659898475 | 38 | -11.783299785974805 | 5.5430657455111653e-09 | 3.3901457391067851e-15 | 5.0500172350014425e-17 | PASS |
| 36 | radial_profile_L2 | null | null | 0.19113782511778182 | 0.021200415960106734 | 9.0157582510385588 | 12.937967455874144 | 1.5352977834313781e-09 | 2.3205167683510447e-21 | 1.873559415950491e-13 | PASS |
| 40 | total_mass | 421.125 | 406.62557094179493 | 0.034430226318088571 | 0.0014841199168892847 | 23.199086493128078 | -12.123072892031722 | 3.7595679941645351e-09 | 1.454436824969024e-26 | 8.6049179868181217e-23 | PASS |
| 40 | r_K_ratio | 1.4862818219349081 | 1.6995080743478872 | 0.1434628677187155 | 0.0052088165665525443 | 27.542315204558335 | 16.515053506398107 | 4.9544479197942529e-11 | 1.9779308949787462e-22 | 1.6280039013798181e-11 | PASS |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0044843049327354259 | 14 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 40 | r90 | 21.5 | 19.0625 | 0.11337209302325581 | 0 | null | -19.030051422496395 | 6.4760803247711866e-12 | 6.3605411260372649e-15 | 2.3244786420838592e-18 | PASS |
| 40 | r99 | 25 | 22.6875 | 0.092499999999999999 | 0.0074999999999999997 | 12.333333333333334 | -15.363413772941893 | 1.3839300528290396e-10 | 5.4380328870719136e-17 | 4.8751207713627678e-19 | PASS |
| 40 | radial_profile_L2 | null | null | 0.20294210478022673 | 0.019135549966259355 | 10.605501547541785 | 12.331063058139474 | 2.9774061335968384e-09 | 2.6639476681283848e-22 | 9.5858061539418275e-14 | PASS |
| 44 | total_mass | 433 | 439.12414900465899 | 0.014143531188588849 | 0.0072170900692840644 | 1.9597276814908708 | 4.3124362386774742 | 0.00061619879125071012 | 5.6783101501054603e-26 | 3.3785584158268306e-21 | PASS |
| 44 | r_K_ratio | 1.3851122690599202 | 1.4826505247964381 | 0.070419025168780944 | 0.023883569861476831 | 2.948429634983663 | 10.150700193299443 | 4.1027820504537489e-08 | 1.7459206768544517e-22 | 6.5239376359162711e-15 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 14.0625 | 13.0625 | 0.071111111111111111 | 0.0088888888888888889 | 8 | -7.7459666924148332 | 1.2783477850604916e-06 | 1.25956553381152e-13 | 5.7635238603637946e-15 | PASS |
| 44 | r90 | 21.875 | 19.1875 | 0.12285714285714286 | 0.0057142857142857143 | 21.5 | -22.456017618285021 | 5.8544703070058476e-13 | 4.6040517988789706e-14 | 9.7976516786557512e-20 | PASS |
| 44 | r99 | 25.4375 | 22.9375 | 0.098280098280098274 | 0.0024570024570024569 | 40 | -15.811388300841896 | 9.2073297610083179e-11 | 1.32240574860815e-17 | 2.0753132400270518e-18 | PASS |
| 44 | radial_profile_L2 | null | null | 0.20578903863170117 | 0.022666905945975712 | 9.0788323347782285 | 13.563001661788283 | 7.9734125142345473e-10 | 3.9888272112335183e-23 | 2.1411534627969286e-14 | PASS |
| 48 | total_mass | 468.375 | 495.63025101823672 | 0.058191088376272695 | 0.0066720042700827327 | 8.7216803258357523 | 13.932439671958022 | 5.478965023639955e-10 | 6.2262438154028136e-25 | 9.8983062404198524e-19 | PASS |
| 48 | r_K_ratio | 1.4233832881828432 | 1.4168010867127527 | 0.00462433521928831 | 0.023969344303723984 | 0.19292706386506586 | -0.56934483893448984 | 0.5775492092836918 | 5.7033749883260893e-20 | 1.7002683649562875e-15 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 14.4375 | 13.1875 | 0.086580086580086577 | 0.030303030303030304 | 2.8571428571428572 | -11.180339887498949 | 1.129731174656377e-08 | 5.7520127976418164e-14 | 2.7873182537487998e-16 | PASS |
| 48 | r90 | 22.375 | 19.5 | 0.12849162011173185 | 0 | null | -15.998991903725805 | 7.7865353236989104e-11 | 1.0061554522861548e-11 | 4.8264550809103502e-17 | PASS |
| 48 | r99 | 26.1875 | 23.1875 | 0.11455847255369929 | 0.0071599045346062056 | 16 | -23.2379000772445 | 3.5513304925948273e-13 | 1.7995348202025614e-17 | 8.6169881800138773e-21 | PASS |
| 48 | radial_profile_L2 | null | null | 0.23341720836842744 | 0.017773720507827693 | 13.1327151378142 | 12.748408916387797 | 1.8827621113569947e-09 | 3.1119166670442979e-21 | 6.3693627289909649e-11 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1875, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

#### legacy_growth_mean: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 255.99987439227883 | 4.9065516088270256e-07 | 0 | null | -1.0352175095070719 | 0.31696944719393905 | 6.3310518845588904e-81 | 6.3307856304481416e-81 | PASS |
| 4 | r_K_ratio | 1 | 0.9999990186955926 | 9.8130440741306391e-07 | 0 | null | -1.0352112711393178 | 0.31697226570925685 | 2.0746026777977294e-76 | 2.0744281865575098e-76 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13.0625 | 13.0625 | 0 | 0.0047846889952153108 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 4 | r90 | 17.9375 | 17.6875 | 0.013937282229965157 | 0.0069686411149825784 | 2 | -2.2360679774997898 | 0.040968955955836141 | 1.1148111900685798e-16 | 3.8956480726576652e-17 | PASS |
| 4 | r99 | 19.25 | 19.375 | 0.0064935064935064939 | 0.00974025974025974 | 0.66666666666666663 | 1 | 0.33317013591547739 | 1.2617854706008803e-18 | 7.1507823999893397e-17 | PASS |
| 4 | radial_profile_L2 | null | null | 0.01390051111786134 | 0.026445474852675486 | 0.52562909894034393 | -7.7751682619068863 | 1.2210780923911863e-06 | 1.7063425300474121e-19 | 5.5915970512964533e-19 | PASS |
| 8 | total_mass | 256 | 255.99987439079101 | 4.9066097267125297e-07 | 0 | null | -1.0352297764902778 | 0.31696390498294125 | 6.331051438669818e-81 | 6.3307851814239781e-81 | PASS |
| 8 | r_K_ratio | 1 | 0.99999901869484664 | 9.8130515333721968e-07 | 0 | null | -1.0352120739000901 | 0.3169719030182524 | 2.0746022011006527e-76 | 2.0744277097679295e-76 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 8 | r90 | 18.0625 | 18 | 0.0034602076124567475 | 0.0034602076124567475 | 1 | -1 | 0.33317013591547739 | 1.4301148022486204e-22 | 1.9701065876166452e-19 | PASS |
| 8 | r99 | 20.0625 | 20 | 0.0031152647975077881 | 0.0031152647975077881 | 1 | -1 | 0.33317013591547739 | 6.5149803197211003e-25 | 5.0766352092358305e-21 | PASS |
| 8 | radial_profile_L2 | null | null | 0.067252396223438593 | 0.031674679515173401 | 2.1232226261743889 | -1.1211565613575198 | 0.2798500299122274 | 1.1233733444798528e-17 | 3.6201109868020346e-15 | PASS |
| 12 | total_mass | 256 | 255.99987438916182 | 4.9066733665698026e-07 | 0 | null | -1.0352432081050333 | 0.31695783667291111 | 6.3310510325227097e-81 | 6.3307847718405203e-81 | PASS |
| 12 | r_K_ratio | 1 | 0.9999990186944866 | 9.8130551335479099e-07 | 0 | null | -1.0352124639278295 | 0.316971726802036 | 2.0746018935233379e-76 | 2.0744274021525138e-76 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 12 | r90 | 18.6875 | 18 | 0.03678929765886288 | 0 | null | -5.7445626465380286 | 3.8761814875216441e-05 | 1.2411344634416193e-17 | 3.0051609001276733e-16 | PASS |
| 12 | r99 | 20.5625 | 20.125 | 0.021276595744680851 | 0.015197568389057751 | 1.3999999999999999 | -3.415650255319866 | 0.0038327022878305922 | 3.2453744019125243e-19 | 3.5829058514479523e-17 | PASS |
| 12 | radial_profile_L2 | null | null | 0.07289369138475793 | 0.033001760468462815 | 2.2087819058748908 | -0.25736842332757609 | 0.80039152143416259 | 4.3685958447204601e-20 | 2.3969087135269472e-17 | PASS |
| 16 | total_mass | 256 | 256.0043048668507 | 1.6815886135668201e-05 | 0 | null | 25.281090489138435 | 1.0334452942718925e-13 | 1.0205672990545397e-78 | 1.0220393662766436e-78 | PASS |
| 16 | r_K_ratio | 1 | 1.0000336318141791 | 3.3631814179216346e-05 | 0 | null | 25.281122633742193 | 1.0334260097149979e-13 | 3.3417845866988837e-74 | 3.3514319378064038e-74 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.0047619047619047623 | 1 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 16 | r90 | 19.0625 | 18 | 0.055737704918032788 | 0.0032786885245901639 | 17 | -17 | 3.2769117377966068e-11 | 2.2623006381054966e-21 | 5.3331248280097439e-21 | PASS |
| 16 | r99 | 21.1875 | 20.8125 | 0.017699115044247787 | 0.014749262536873156 | 1.2 | -2.4227185592617446 | 0.028528068449741654 | 4.6950640203719188e-18 | 3.0093626661365215e-16 | PASS |
| 16 | radial_profile_L2 | null | null | 0.087753989476160837 | 0.029422120612765671 | 2.9825854713574294 | 0.20759079081507811 | 0.83834056220867414 | 5.9176043809817526e-20 | 1.2295350420357534e-16 | PASS |
| 20 | total_mass | 256 | 258.87708144046968 | 0.01123859937683469 | 0 | null | 47.865717680034862 | 8.0718266803171438e-18 | 1.046401918589985e-40 | 2.7428326630653689e-40 | PASS |
| 20 | r_K_ratio | 1 | 1.0224771623306497 | 0.022477162330649558 | 0 | null | 47.865676358154623 | 8.0719305657427456e-18 | 2.1652465987556308e-36 | 1.4905853712380573e-35 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0 | null | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 20 | r90 | 19.375 | 18.125 | 0.064516129032258063 | 0 | null | -11.180339887498949 | 1.129731174656377e-08 | 6.9886341424657286e-17 | 1.1645945213682297e-17 | PASS |
| 20 | r99 | 21.9375 | 21 | 0.042735042735042736 | 0 | null | -6.5361700549589266 | 9.4202794798815216e-06 | 3.6683425778815099e-19 | 5.1861867803112176e-17 | PASS |
| 20 | radial_profile_L2 | null | null | 0.092503569583778644 | 0.040888641583758337 | 2.262329243545294 | 1.0412800493061503 | 0.31423893386661061 | 1.1299200850275354e-19 | 3.5986830288696835e-16 | PASS |
| 24 | total_mass | 284.625 | 295.67311054993434 | 0.038816374351987148 | 0.0083443126921387799 | 4.6518360210249865 | 15.265995292494356 | 1.5142857579525539e-10 | 1.6107229299610676e-27 | 2.416144962238214e-22 | PASS |
| 24 | r_K_ratio | 1.2236328125 | 1.3099217152887528 | 0.070518624465828197 | 0.015163607342378291 | 4.6505177081938287 | 15.261699878099721 | 1.5203257829393643e-10 | 3.8496103631446465e-24 | 9.0668949975115367e-18 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.5625 | 13.0625 | 0.036866359447004608 | 0.018433179723502304 | 2 | -3.8729833462074166 | 0.0015017735386984323 | 9.6072681659644597e-15 | 7.6467248765302002e-14 | PASS |
| 24 | r90 | 19.9375 | 18.25 | 0.084639498432601878 | 0.003134796238244514 | 27 | -14.100290132411525 | 4.6333644836412128e-10 | 7.7245411315955695e-15 | 9.481589102890017e-19 | PASS |
| 24 | r99 | 22.875 | 21.25 | 0.071038251366120214 | 0.0054644808743169399 | 13 | -6.7890285822722154 | 6.1053729318667233e-06 | 1.2337518345797697e-14 | 7.416708381213179e-15 | PASS |
| 24 | radial_profile_L2 | null | null | 0.13147638879830365 | 0.031998880089473554 | 4.1087809457917404 | 4.1522696561222547 | 0.00085120132006816891 | 3.4084957453872834e-19 | 4.2064121648031338e-14 | PASS |
| 28 | total_mass | 360.0625 | 353.09782162515279 | 0.019342970664390734 | 0.0034716195105016492 | 5.5717426998777508 | -8.161096187531987 | 6.7304060788605287e-07 | 8.1307857338724048e-27 | 7.7339332178511622e-24 | PASS |
| 28 | r_K_ratio | 1.81298828125 | 1.7567071926338482 | 0.031043272148095605 | 0.0053864799353622404 | 5.7631834742939487 | -8.4487206340487813 | 4.3691968704139455e-07 | 9.7875648471756778e-24 | 3.4731767895840761e-21 | PASS |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.875 | 13.0625 | 0.058558558558558557 | 0.022522522522522521 | 2.6000000000000001 | -8.0622577482985509 | 7.8258249630702749e-07 | 2.2778649613890186e-15 | 2.5893641111623486e-16 | PASS |
| 28 | r90 | 20.5625 | 18.75 | 0.088145896656534953 | 0.0060790273556231003 | 14.5 | -11.066875266096559 | 1.2961577489987425e-08 | 1.5680800107416445e-13 | 1.9477056267575937e-16 | PASS |
| 28 | r99 | 23.625 | 21.9375 | 0.071428571428571425 | 0 | null | -8.5098307000225546 | 3.9910679004152207e-07 | 1.2402986697836266e-16 | 5.9545457575235637e-16 | PASS |
| 28 | radial_profile_L2 | null | null | 0.16828046394019908 | 0.025528598153441336 | 6.5918411551131078 | 8.2277534081340384 | 6.083713433054982e-07 | 2.1126918139331655e-20 | 1.2090508852925571e-13 | PASS |
| 32 | total_mass | 385.5625 | 378.72153066410203 | 0.017742828558010546 | 0.0024315124007132437 | 7.2970339582911379 | -5.6687269916952081 | 4.4587525660575328e-05 | 2.1920019954162635e-25 | 2.4846359304662547e-22 | PASS |
| 32 | r_K_ratio | 1.8212005047589639 | 1.9274240674598628 | 0.058326121930741277 | 0.0014984898683687194 | 38.923267458749017 | 13.326733304203259 | 1.018274515780266e-09 | 3.3753979921218085e-25 | 2.9307059479931345e-18 | PASS |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.022421524663677129 | 2.7999999999999998 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 32 | r90 | 20.8125 | 18.9375 | 0.090090090090090086 | 0.006006006006006006 | 15 | -21.957751641341996 | 8.1207854767429677e-13 | 2.5303126792687863e-18 | 2.9422487852212645e-20 | PASS |
| 32 | r99 | 24.375 | 22 | 0.097435897435897437 | 0.0051282051282051282 | 19 | -15.343884208657718 | 1.4090734244550085e-10 | 8.3953741872725587e-18 | 3.7020237911275972e-18 | PASS |
| 32 | radial_profile_L2 | null | null | 0.1794857719734804 | 0.026253620259367155 | 6.8366103493647055 | 10.823536860150758 | 1.7470767308970091e-08 | 2.865706124766142e-21 | 5.9396786326554099e-14 | PASS |
| 36 | total_mass | 404.5 | 391.08242447311665 | 0.033170767680799351 | 0.001854140914709518 | 17.890100702511116 | -10.530655744457547 | 2.5202341654704299e-08 | 1.2559770490901061e-25 | 2.8705589055420151e-22 | PASS |
| 36 | r_K_ratio | 1.6452597661040846 | 1.8667140903205504 | 0.13460143424090509 | 0.0025243036803777237 | 53.322203381157387 | 16.810212895812448 | 3.8473914145124943e-11 | 1.2955422517217467e-22 | 2.0669110244284775e-12 | PASS |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0089686098654708519 | 7 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 36 | r90 | 21.3125 | 19 | 0.10850439882697947 | 0 | null | -19.322619811082461 | 5.195040620846464e-12 | 7.7864380968352111e-16 | 1.5335258890741899e-18 | PASS |
| 36 | r99 | 24.625 | 22.25 | 0.096446700507614211 | 0.0025380710659898475 | 38 | -11.783299785974805 | 5.5430657455111653e-09 | 3.3901457391067851e-15 | 5.0500172350014425e-17 | PASS |
| 36 | radial_profile_L2 | null | null | 0.19136126810388571 | 0.021200415960106734 | 9.0262978077399136 | 12.934746003801171 | 1.5405975393385223e-09 | 2.3732601250895472e-21 | 1.9667646951546544e-13 | PASS |
| 40 | total_mass | 421.125 | 410.82317301950906 | 0.024462634563350408 | 0.0014841199168892847 | 16.482923168785504 | -8.7119968556387093 | 2.9675094272265218e-07 | 8.1079733217685101e-27 | 1.0927942833861877e-22 | PASS |
| 40 | r_K_ratio | 1.4862818219349081 | 1.6511495584847613 | 0.11092629548225309 | 0.0052088165665525443 | 21.295872884936255 | 12.947495404433035 | 1.5197358430371177e-09 | 4.227627470605292e-22 | 1.6955397343341467e-12 | PASS |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0044843049327354259 | 14 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 40 | r90 | 21.5 | 19.0625 | 0.11337209302325581 | 0 | null | -19.030051422496395 | 6.4760803247711866e-12 | 6.3605411260372649e-15 | 2.3244786420838592e-18 | PASS |
| 40 | r99 | 25 | 22.625 | 0.095000000000000001 | 0.0074999999999999997 | 12.666666666666666 | -15.343884208657718 | 1.4090734244550085e-10 | 1.1392578071103199e-16 | 5.9909545657122436e-19 | PASS |
| 40 | radial_profile_L2 | null | null | 0.20344044299207009 | 0.019135549966259355 | 10.631544081606499 | 12.376038368950912 | 2.8321744740201371e-09 | 2.6925111769781923e-22 | 1.0315292932798695e-13 | PASS |
| 44 | total_mass | 433 | 443.66891444890086 | 0.024639525286145166 | 0.0072170900692840644 | 3.4140526236482742 | 7.563079679294102 | 1.7076285752582351e-06 | 3.4983526725905944e-26 | 4.9312756871498478e-21 | PASS |
| 44 | r_K_ratio | 1.3851122690599202 | 1.424170951968994 | 0.028198929272053221 | 0.023883569861476831 | 1.1806831824390238 | 4.045431960035101 | 0.0010570582978787631 | 9.7167156177911812e-22 | 8.5326598152564426e-16 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 14.0625 | 13.0625 | 0.071111111111111111 | 0.0088888888888888889 | 8 | -7.7459666924148332 | 1.2783477850604916e-06 | 1.25956553381152e-13 | 5.7635238603637946e-15 | PASS |
| 44 | r90 | 21.875 | 19.125 | 0.12571428571428572 | 0.0057142857142857143 | 22 | -24.596747752497688 | 1.5455927408310337e-13 | 1.5653444810397694e-14 | 4.3711393676092139e-20 | PASS |
| 44 | r99 | 25.4375 | 22.9375 | 0.098280098280098274 | 0.0024570024570024569 | 40 | -15.811388300841896 | 9.2073297610083179e-11 | 1.32240574860815e-17 | 2.0753132400270518e-18 | PASS |
| 44 | radial_profile_L2 | null | null | 0.20683746057612379 | 0.022666905945975712 | 9.1250857558195211 | 13.648045699485785 | 7.308001845110556e-10 | 4.269046298900404e-23 | 2.6196580378288016e-14 | PASS |
| 48 | total_mass | 468.375 | 495.44252848160704 | 0.05779029299515779 | 0.0066720042700827327 | 8.6616091141142491 | 13.775742598573697 | 6.4172861336992929e-10 | 6.1730863039587142e-25 | 1.0516027242724572e-18 | PASS |
| 48 | r_K_ratio | 1.4233832881828432 | 1.3477040747870106 | 0.053168541477291088 | 0.023969344303723984 | 2.2181892338636309 | -6.4979728201011548 | 1.0065778784167947e-05 | 5.1558544665053705e-19 | 2.8948839113440579e-16 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 14.4375 | 13.1875 | 0.086580086580086577 | 0.030303030303030304 | 2.8571428571428572 | -11.180339887498949 | 1.129731174656377e-08 | 5.7520127976418164e-14 | 2.7873182537487998e-16 | PASS |
| 48 | r90 | 22.375 | 19.5 | 0.12849162011173185 | 0 | null | -15.998991903725805 | 7.7865353236989104e-11 | 1.0061554522861548e-11 | 4.8264550809103502e-17 | PASS |
| 48 | r99 | 26.1875 | 23.1875 | 0.11455847255369929 | 0.0071599045346062056 | 16 | -23.2379000772445 | 3.5513304925948273e-13 | 1.7995348202025614e-17 | 8.6169881800138773e-21 | PASS |
| 48 | radial_profile_L2 | null | null | 0.23491655815817461 | 0.017773720507827693 | 13.217072815717756 | 12.821834077613682 | 1.7391789565352464e-09 | 3.2349784172690698e-21 | 8.2290308324627074e-11 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1875, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

#### axial_normal_transport: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 255.99952363215377 | 1.8608118993179712e-06 | 0 | null | -1.051936204319238 | 0.30948091139523803 | 2.4045068684167232e-72 | 2.4041233845721485e-72 | PASS |
| 4 | r_K_ratio | 1 | 0.99999627838045557 | 3.7216195444764177e-06 | 0 | null | -1.0519350021930869 | 0.30948144517619536 | 7.8797163836281917e-68 | 7.8772031868028145e-68 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13.0625 | 13.0625 | 0 | 0.0047846889952153108 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 4 | r90 | 17.9375 | 17.6875 | 0.013937282229965157 | 0.0069686411149825784 | 2 | -2.2360679774997898 | 0.040968955955836141 | 1.1148111900685798e-16 | 3.8956480726576652e-17 | PASS |
| 4 | r99 | 19.25 | 19.25 | 0 | 0.00974025974025974 | 0 | 0 | 1 | 1.6302193524528563e-20 | 9.3136886168229159e-19 | PASS |
| 4 | radial_profile_L2 | null | null | 0.018038965598942157 | 0.026445474852675486 | 0.68211917915768328 | -6.9561301233640789 | 4.6059224602554363e-06 | 1.1097248239271743e-18 | 5.1697635186980713e-18 | PASS |
| 8 | total_mass | 256 | 255.99952363105518 | 1.8608161906907839e-06 | 0 | null | -1.0519386328313456 | 0.30947983306326909 | 2.4045067812667786e-72 | 2.4041232965518539e-72 | PASS |
| 8 | r_K_ratio | 1 | 0.99999627837972738 | 3.721620272581494e-06 | 0 | null | -1.0519352085026517 | 0.30948135356836082 | 7.8797163269185298e-68 | 7.8772031296198129e-68 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 8 | r90 | 18.0625 | 17.9375 | 0.006920415224913495 | 0.0034602076124567475 | 2 | -1.4638501094227998 | 0.16387561365565406 | 2.2451366252918472e-19 | 3.9072803388172492e-18 | PASS |
| 8 | r99 | 20.0625 | 20 | 0.0031152647975077881 | 0.0031152647975077881 | 1 | -1 | 0.33317013591547739 | 6.5149803197211003e-25 | 5.0766352092358305e-21 | PASS |
| 8 | radial_profile_L2 | null | null | 0.076347291770587397 | 0.031674679515173401 | 2.4103571982161358 | -0.25602421236708434 | 0.80140998523340368 | 4.1542646456033906e-17 | 2.940075560798739e-14 | PASS |
| 12 | total_mass | 256 | 255.99952362987185 | 1.8608208131112858e-06 | 0 | null | -1.0519412461377775 | 0.30947867268015405 | 2.4045067747373516e-72 | 2.4041232890709577e-72 | PASS |
| 12 | r_K_ratio | 1 | 0.99999627838038463 | 3.721619615419669e-06 | 0 | null | -1.051935024985529 | 0.30948143505564507 | 7.879716075788825e-68 | 7.8772028790139432e-68 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 12 | r90 | 18.6875 | 18 | 0.03678929765886288 | 0 | null | -5.7445626465380286 | 3.8761814875216441e-05 | 1.2411344634416193e-17 | 3.0051609001276733e-16 | PASS |
| 12 | r99 | 20.5625 | 20 | 0.027355623100303952 | 0.015197568389057751 | 1.8 | -4.3915503282683988 | 0.00052573078993947925 | 7.4736080236482832e-20 | 4.9956120358567663e-17 | PASS |
| 12 | radial_profile_L2 | null | null | 0.086013650856360188 | 0.033001760468462815 | 2.6063352268299949 | 0.83491828510420141 | 0.41686386867330738 | 5.7822343406062139e-19 | 1.0101692911077732e-15 | PASS |
| 16 | total_mass | 256 | 256.00526755326501 | 2.0576379941475431e-05 | 0 | null | 10.813974594958577 | 1.7678838910466133e-08 | 7.1724366191937527e-72 | 7.1850977337153706e-72 | PASS |
| 16 | r_K_ratio | 1 | 1.0000411527916733 | 4.1152791673298994e-05 | 0 | null | 10.813983010450089 | 1.7678654643032528e-08 | 2.3481923329700558e-67 | 2.3564899304942474e-67 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.0047619047619047623 | 1 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 16 | r90 | 19.0625 | 18 | 0.055737704918032788 | 0.0032786885245901639 | 17 | -17 | 3.2769117377966068e-11 | 2.2623006381054966e-21 | 5.3331248280097439e-21 | PASS |
| 16 | r99 | 21.1875 | 20.3125 | 0.041297935103244837 | 0.014749262536873156 | 2.7999999999999998 | -4.8692584054817667 | 0.00020423834070405623 | 2.6523458183220211e-16 | 6.1694084828562171e-16 | PASS |
| 16 | radial_profile_L2 | null | null | 0.1048251148974642 | 0.029422120612765671 | 3.5627994418588127 | 1.4430290693048549 | 0.16956995678679132 | 2.3449520353673083e-19 | 2.3042158871135368e-15 | PASS |
| 20 | total_mass | 256 | 259.37081008037609 | 0.013167226876469151 | 0 | null | 48.215782923809847 | 7.2408678647785973e-18 | 9.3167974770322864e-40 | 2.8816887094645916e-39 | PASS |
| 20 | r_K_ratio | 1 | 1.0263344334993341 | 0.026334433499334106 | 0 | null | 48.21576356862019 | 7.2409112021921303e-18 | 1.7890963770194272e-35 | 1.7169652989695532e-34 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0 | null | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 20 | r90 | 19.375 | 18 | 0.070967741935483872 | 0 | null | -11 | 1.4062516106729139e-08 | 1.8615816799604571e-16 | 6.2976132309167012e-17 | PASS |
| 20 | r99 | 21.9375 | 20.8125 | 0.05128205128205128 | 0 | null | -7.2681556777852343 | 2.7487748593004314e-06 | 1.4996759808931379e-17 | 4.7769350818866047e-17 | PASS |
| 20 | radial_profile_L2 | null | null | 0.11206543051339203 | 0.040888641583758337 | 2.740747214206948 | 2.5783059073659937 | 0.020985207826472028 | 4.8986791257649807e-19 | 9.3826858503428699e-15 | PASS |
| 24 | total_mass | 284.625 | 298.85191296613056 | 0.049984762287678779 | 0.0083443126921387799 | 5.9902791436339253 | 19.544713476455406 | 4.4036920802131299e-12 | 1.6721366288231709e-27 | 3.9588891326149724e-22 | PASS |
| 24 | r_K_ratio | 1.2236328125 | 1.3347655277348025 | 0.090821947645999823 | 0.015163607342378291 | 5.9894684421283042 | 19.542117936131536 | 4.412160873778029e-12 | 3.0139351588142226e-24 | 2.6427272380398998e-17 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.5625 | 13.0625 | 0.036866359447004608 | 0.018433179723502304 | 2 | -3.8729833462074166 | 0.0015017735386984323 | 9.6072681659644597e-15 | 7.6467248765302002e-14 | PASS |
| 24 | r90 | 19.9375 | 18.0625 | 0.094043887147335428 | 0.003134796238244514 | 30 | -21.957751641341996 | 8.1207854767429677e-13 | 3.8511104687127252e-17 | 1.1403537175415299e-20 | PASS |
| 24 | r99 | 22.875 | 21 | 0.081967213114754092 | 0.0054644808743169399 | 15 | -9.3026050941906355 | 1.2823132219716826e-07 | 3.6636365223864341e-16 | 8.7002119456432826e-16 | PASS |
| 24 | radial_profile_L2 | null | null | 0.15475340746438462 | 0.031998880089473554 | 4.8362132372030349 | 6.2080230349262155 | 1.6756528385413101e-05 | 8.254505345088776e-19 | 1.0321512154211331e-12 | PASS |
| 28 | total_mass | 360.0625 | 355.52846183116498 | 0.012592364294629484 | 0.0034716195105016492 | 3.6272305350680227 | -5.4277169865740404 | 6.9934605036360041e-05 | 2.8832917209595357e-27 | 9.0356243054192236e-24 | PASS |
| 28 | r_K_ratio | 1.81298828125 | 1.7763155190150368 | 0.020227798830380977 | 0.0053864799353622404 | 3.7552908528602287 | -5.6219983673717495 | 4.8624010335170369e-05 | 2.8851635547785351e-24 | 4.747337350319515e-21 | PASS |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.875 | 13.0625 | 0.058558558558558557 | 0.022522522522522521 | 2.6000000000000001 | -8.0622577482985509 | 7.8258249630702749e-07 | 2.2778649613890186e-15 | 2.5893641111623486e-16 | PASS |
| 28 | r90 | 20.5625 | 18.1875 | 0.11550151975683891 | 0.0060790273556231003 | 19 | -19 | 6.625544366525624e-12 | 3.1383809442133001e-14 | 1.7974071832715349e-18 | PASS |
| 28 | r99 | 23.625 | 21.3125 | 0.097883597883597878 | 0 | null | -10.593069184806563 | 2.3293939237721809e-08 | 1.1628481028638976e-14 | 4.8813247618385495e-16 | PASS |
| 28 | radial_profile_L2 | null | null | 0.19635651148037808 | 0.025528598153441336 | 7.6916292191276705 | 10.101880366381172 | 4.3722967238119133e-08 | 1.1850468520539632e-19 | 1.5717355376421298e-11 | PASS |
| 32 | total_mass | 385.5625 | 378.62878488237055 | 0.017983375244297488 | 0.0024315124007132437 | 7.3959627921380804 | -5.9833779902730022 | 2.5068494207002363e-05 | 8.9953736732024445e-26 | 1.5674202209301484e-22 | PASS |
| 32 | r_K_ratio | 1.8212005047589639 | 1.9353860899973288 | 0.062697975835163386 | 0.0014984898683687194 | 41.840773940912541 | 14.436562802220847 | 3.3287627952940745e-10 | 1.9670682299690936e-25 | 3.5355828040475571e-18 | PASS |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.022421524663677129 | 2.7999999999999998 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 32 | r90 | 20.8125 | 18.375 | 0.11711711711711711 | 0.006006006006006006 | 19.5 | -19.030051422496395 | 6.4760803247711866e-12 | 1.9103560717606373e-13 | 6.0923273587867951e-19 | PASS |
| 32 | r99 | 24.375 | 21.8125 | 0.10512820512820513 | 0.0051282051282051282 | 20.5 | -16.291747991900039 | 6.0154840401030353e-11 | 1.0891629608238304e-16 | 1.894984565766644e-18 | PASS |
| 32 | radial_profile_L2 | null | null | 0.21017692173453367 | 0.026253620259367155 | 8.0056357812040666 | 13.504841041681912 | 8.46533013782107e-10 | 1.1785520475251349e-20 | 9.4989822847516493e-12 | PASS |
| 36 | total_mass | 404.5 | 388.87640144650538 | 0.038624471084041039 | 0.001854140914709518 | 20.8314647379928 | -12.217703521842136 | 3.3796496177522227e-09 | 1.723490525880749e-25 | 2.4339445109858802e-22 | PASS |
| 36 | r_K_ratio | 1.6452597661040846 | 1.89023903945492 | 0.14890005724198646 | 0.0025243036803777237 | 58.986586439435783 | 18.704585438236595 | 8.3063867817411075e-12 | 6.9304230898565987e-23 | 5.2963912414831648e-12 | PASS |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0089686098654708519 | 7 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 36 | r90 | 21.3125 | 18.625 | 0.12609970674486803 | 0 | null | -17.854778168554091 | 1.6220160419503157e-11 | 2.1663433333220131e-12 | 5.0379176673465418e-18 | PASS |
| 36 | r99 | 24.625 | 22 | 0.1065989847715736 | 0.0025380710659898475 | 42 | -16.959029914832215 | 3.3919667718492142e-11 | 1.4370664950489062e-17 | 2.2632691711460884e-18 | PASS |
| 36 | radial_profile_L2 | null | null | 0.22520364126933659 | 0.021200415960106734 | 10.622604843843961 | 16.557179848562978 | 4.7776253506237389e-11 | 1.0041190496382848e-20 | 6.0215508883229342e-11 | PASS |
| 40 | total_mass | 421.125 | 405.8726126927375 | 0.036218194852508191 | 0.0014841199168892847 | 24.403819691620022 | -12.735344945237809 | 1.9096007744660321e-09 | 1.7174725205750857e-26 | 8.0611989384802897e-23 | PASS |
| 40 | r_K_ratio | 1.4862818219349081 | 1.6981230535714129 | 0.14253099816609513 | 0.0052088165665525443 | 27.363412849154958 | 16.55512964898778 | 4.786073731169937e-11 | 1.6295438955927735e-22 | 1.3944986306706897e-11 | PASS |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0044843049327354259 | 14 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 40 | r90 | 21.5 | 18.8125 | 0.125 | 0 | null | -22.456017618285021 | 5.8544703070058476e-13 | 2.3582850886443684e-14 | 3.6683425778814834e-19 | PASS |
| 40 | r99 | 25 | 22 | 0.12 | 0.0074999999999999997 | 16 | -23.2379000772445 | 3.5513304925948273e-13 | 2.2385921293450308e-18 | 7.4091209429095752e-20 | PASS |
| 40 | radial_profile_L2 | null | null | 0.23899995294026577 | 0.019135549966259355 | 12.489839767431874 | 15.13068812307853 | 1.7174123746578507e-10 | 8.0050480277539244e-22 | 4.0379933091621905e-11 | PASS |
| 44 | total_mass | 433 | 437.09446333277987 | 0.0094560354105770184 | 0.0072170900692840644 | 1.3102282664895517 | 2.9271288709690571 | 0.010405804071278588 | 6.1224620907034733e-26 | 2.1971081179948403e-21 | PASS |
| 44 | r_K_ratio | 1.3851122690599202 | 1.4772987589571529 | 0.066555247510586155 | 0.023883569861476831 | 2.7866540846532706 | 9.6486126743403862 | 7.9833914895397646e-08 | 1.9549751240246803e-22 | 4.9066701258011471e-15 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 14.0625 | 13.0625 | 0.071111111111111111 | 0.0088888888888888889 | 8 | -7.7459666924148332 | 1.2783477850604916e-06 | 1.25956553381152e-13 | 5.7635238603637946e-15 | PASS |
| 44 | r90 | 21.875 | 18.9375 | 0.13428571428571429 | 0.0057142857142857143 | 23.5 | -47 | 1.0595740345741668e-17 | 1.9641369866418272e-18 | 2.4448502950007348e-23 | PASS |
| 44 | r99 | 25.4375 | 22.375 | 0.12039312039312039 | 0.0024570024570024569 | 49 | -21.35148884619916 | 1.2208711374192415e-12 | 9.3308929176533803e-17 | 1.3651351700783063e-19 | PASS |
| 44 | radial_profile_L2 | null | null | 0.24362775013918322 | 0.022666905945975712 | 10.748169631966773 | 15.747172026832581 | 9.7549904941400463e-11 | 1.3965905219659489e-22 | 1.5709829498195762e-11 | PASS |
| 48 | total_mass | 468.375 | 489.66129335500079 | 0.04544711685081574 | 0.0066720042700827327 | 6.8116138736002627 | 11.246211669170751 | 1.0436335293909497e-08 | 6.2154756329781409e-25 | 3.4422212994246096e-19 | PASS |
| 48 | r_K_ratio | 1.4233832881828432 | 1.3957960505366345 | 0.019381453945147738 | 0.023969344303723984 | 0.80859341413592811 | -2.4360973075961376 | 0.027789558624606973 | 6.7448200236472116e-20 | 7.5204583909274865e-16 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 14.4375 | 13.1875 | 0.086580086580086577 | 0.030303030303030304 | 2.8571428571428572 | -11.180339887498949 | 1.129731174656377e-08 | 5.7520127976418164e-14 | 2.7873182537487998e-16 | PASS |
| 48 | r90 | 22.375 | 19 | 0.15083798882681565 | 0 | null | -21.804467033355703 | 8.9935030956328787e-13 | 3.0419385484799993e-12 | 6.5100670518587811e-18 | PASS |
| 48 | r99 | 26.1875 | 22.9375 | 0.12410501193317422 | 0.0071599045346062056 | 17.333333333333332 | -29.068883707497267 | 1.3242683256898944e-14 | 1.1102156212733309e-18 | 1.7102375145581811e-21 | PASS |
| 48 | radial_profile_L2 | null | null | 0.27325529942775689 | 0.017773720507827693 | 15.374119296374275 | 15.399492908167378 | 1.338725462620493e-10 | 2.8633870301260636e-21 | 4.4060863891113401e-08 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1796875, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

#### unconditional_lattice: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 255.99980074315673 | 7.7834704404616817e-07 | 0 | null | -1.0427566367979151 | 0.31357648155307283 | 5.7565801952177943e-78 | 5.7561961551276069e-78 | PASS |
| 4 | r_K_ratio | 1 | 0.99999844331192966 | 1.556688070350476e-06 | 0 | null | -1.0427526135594427 | 0.31357828515080266 | 1.8863789089188973e-73 | 1.886127225406674e-73 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13.0625 | 13.0625 | 0 | 0.0047846889952153108 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 4 | r90 | 17.9375 | 17.5625 | 0.020905923344947737 | 0.0069686411149825784 | 3 | -3 | 0.0089727374772233335 | 8.2904737133335303e-16 | 1.4323125621009057e-16 | PASS |
| 4 | r99 | 19.25 | 19.1875 | 0.003246753246753247 | 0.00974025974025974 | 0.33333333333333331 | -0.56493268286603204 | 0.58047096825964761 | 1.3011513422083216e-19 | 1.5001754470731887e-17 | PASS |
| 4 | radial_profile_L2 | null | null | 0.017039957693150063 | 0.026445474852675486 | 0.64434304122265107 | -7.0417721746323876 | 3.9922936176335822e-06 | 4.5476958508649654e-19 | 1.9471703971948076e-18 | PASS |
| 8 | total_mass | 256 | 255.99980074165342 | 7.7835291635575121e-07 | 0 | null | -1.0427645073829865 | 0.31357295323082313 | 5.7565799140089526e-78 | 5.7561958710401534e-78 | PASS |
| 8 | r_K_ratio | 1 | 0.99999844331145593 | 1.5566885440340683e-06 | 0 | null | -1.0427529353217537 | 0.31357814090608693 | 1.8863787878293894e-73 | 1.8861271042567716e-73 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 8 | r90 | 18.0625 | 17.9375 | 0.006920415224913495 | 0.0034602076124567475 | 2 | -1.4638501094227998 | 0.16387561365565406 | 2.2451366252918472e-19 | 3.9072803388172492e-18 | PASS |
| 8 | r99 | 20.0625 | 20 | 0.0031152647975077881 | 0.0031152647975077881 | 1 | -1 | 0.33317013591547739 | 6.5149803197211003e-25 | 5.0766352092358305e-21 | PASS |
| 8 | radial_profile_L2 | null | null | 0.075347012820213347 | 0.031674679515173401 | 2.3787774327477949 | -0.31234157362855258 | 0.75908091817351897 | 4.0884232909808705e-17 | 2.6497653398934852e-14 | PASS |
| 12 | total_mass | 256 | 255.99980074008218 | 7.7835905399475935e-07 | 0 | null | -1.0427727267120674 | 0.31356926860008083 | 5.7565801892955368e-78 | 5.7561961432799547e-78 | PASS |
| 12 | r_K_ratio | 1 | 0.99999844331130572 | 1.5566886943235714e-06 | 0 | null | -1.0427530410583707 | 0.31357809350480892 | 1.8863786504081206e-73 | 1.8861269668295003e-73 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 12 | r90 | 18.6875 | 18 | 0.03678929765886288 | 0 | null | -5.7445626465380286 | 3.8761814875216441e-05 | 1.2411344634416193e-17 | 3.0051609001276733e-16 | PASS |
| 12 | r99 | 20.5625 | 20 | 0.027355623100303952 | 0.015197568389057751 | 1.8 | -4.3915503282683988 | 0.00052573078993947925 | 7.4736080236482832e-20 | 4.9956120358567663e-17 | PASS |
| 12 | radial_profile_L2 | null | null | 0.085369670018223306 | 0.033001760468462815 | 2.5868216969759654 | 0.81243910103514283 | 0.42924588748345838 | 5.9083442979289233e-19 | 9.7389679190910955e-16 | PASS |
| 16 | total_mass | 256 | 256.0055356644134 | 2.1623689114777522e-05 | 0 | null | 22.396708536987919 | 6.0847965339237275e-13 | 2.729459805817243e-76 | 2.7345234390170231e-76 | PASS |
| 16 | r_K_ratio | 1 | 1.0000432474210297 | 4.3247421029665722e-05 | 0 | null | 22.396731097951328 | 6.0847071062715499e-13 | 8.9356075207146653e-72 | 8.9687925884318808e-72 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.0047619047619047623 | 1 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 16 | r90 | 19.0625 | 18 | 0.055737704918032788 | 0.0032786885245901639 | 17 | -17 | 3.2769117377966068e-11 | 2.2623006381054966e-21 | 5.3331248280097439e-21 | PASS |
| 16 | r99 | 21.1875 | 20.3125 | 0.041297935103244837 | 0.014749262536873156 | 2.7999999999999998 | -4.8692584054817667 | 0.00020423834070405623 | 2.6523458183220211e-16 | 6.1694084828562171e-16 | PASS |
| 16 | radial_profile_L2 | null | null | 0.10430220971876546 | 0.029422120612765671 | 3.5450269235015917 | 1.4356531909708659 | 0.17162602402341365 | 2.3554436438387767e-19 | 2.2048217119050509e-15 | PASS |
| 20 | total_mass | 256 | 259.36853158405222 | 0.013158326500204043 | 0 | null | 48.262380330371911 | 7.1373309945550849e-18 | 9.0934518498722605e-40 | 2.8104600100457734e-39 | PASS |
| 20 | r_K_ratio | 1 | 1.0263166328430418 | 0.026316632843041657 | 0 | null | 48.2623611500671 | 7.1373732859300102e-18 | 1.7468055223545254e-35 | 1.6738096572253818e-34 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0 | null | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 20 | r90 | 19.375 | 18 | 0.070967741935483872 | 0 | null | -11 | 1.4062516106729139e-08 | 1.8615816799604571e-16 | 6.2976132309167012e-17 | PASS |
| 20 | r99 | 21.9375 | 20.875 | 0.04843304843304843 | 0 | null | -7.4076593956201169 | 2.1914358166347718e-06 | 2.2187163001862335e-18 | 2.5025837161557336e-17 | PASS |
| 20 | radial_profile_L2 | null | null | 0.11172791585477915 | 0.040888641583758337 | 2.7324927296963413 | 2.5662481789804699 | 0.021493848326032357 | 5.8325466938142778e-19 | 1.0799140687632701e-14 | PASS |
| 24 | total_mass | 284.625 | 298.83943789284558 | 0.049940932429848414 | 0.0083443126921387799 | 5.9850264811981502 | 19.51489154846363 | 4.5020476151587431e-12 | 1.6848860681757457e-27 | 3.9864012147987753e-22 | PASS |
| 24 | r_K_ratio | 1.2236328125 | 1.3346681121848534 | 0.090742335895682161 | 0.015163607342378291 | 5.9842182566994611 | 19.512304953866831 | 4.5106879462100665e-12 | 3.0401380552475736e-24 | 2.6546119183177501e-17 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.5625 | 13.0625 | 0.036866359447004608 | 0.018433179723502304 | 2 | -3.8729833462074166 | 0.0015017735386984323 | 9.6072681659644597e-15 | 7.6467248765302002e-14 | PASS |
| 24 | r90 | 19.9375 | 18.0625 | 0.094043887147335428 | 0.003134796238244514 | 30 | -21.957751641341996 | 8.1207854767429677e-13 | 3.8511104687127252e-17 | 1.1403537175415299e-20 | PASS |
| 24 | r99 | 22.875 | 21 | 0.081967213114754092 | 0.0054644808743169399 | 15 | -9.3026050941906355 | 1.2823132219716826e-07 | 3.6636365223864341e-16 | 8.7002119456432826e-16 | PASS |
| 24 | radial_profile_L2 | null | null | 0.15456049982259223 | 0.031998880089473554 | 4.8301846624137612 | 6.1958453211923112 | 1.7123406501834216e-05 | 9.3457644651446481e-19 | 1.1416916692626425e-12 | PASS |
| 28 | total_mass | 360.0625 | 355.51801735184677 | 0.012621371701172005 | 0.0034716195105016492 | 3.6355861185225962 | -5.4394314584362968 | 6.8409016326756843e-05 | 2.8620586008648878e-27 | 9.0689878187386422e-24 | PASS |
| 28 | r_K_ratio | 1.81298828125 | 1.7762367883228802 | 0.02027122475484543 | 0.0053864799353622404 | 3.7633528757370542 | -5.6332155989104216 | 4.7621630760113733e-05 | 2.8660878181670006e-24 | 4.7625306473508407e-21 | PASS |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.875 | 13.0625 | 0.058558558558558557 | 0.022522522522522521 | 2.6000000000000001 | -8.0622577482985509 | 7.8258249630702749e-07 | 2.2778649613890186e-15 | 2.5893641111623486e-16 | PASS |
| 28 | r90 | 20.5625 | 18.1875 | 0.11550151975683891 | 0.0060790273556231003 | 19 | -19 | 6.625544366525624e-12 | 3.1383809442133001e-14 | 1.7974071832715349e-18 | PASS |
| 28 | r99 | 23.625 | 21.375 | 0.095238095238095233 | 0 | null | -10.50973574618056 | 2.5878480681431328e-08 | 7.6461076184081445e-15 | 4.0573331411132318e-16 | PASS |
| 28 | radial_profile_L2 | null | null | 0.19627180812766637 | 0.025528598153441336 | 7.6883112401222196 | 10.044477668252421 | 4.7133395665176025e-08 | 1.424529372320197e-19 | 1.8556330397462601e-11 | PASS |
| 32 | total_mass | 385.5625 | 378.62402129680345 | 0.017995730142834322 | 0.0024315124007132437 | 7.4010439500763292 | -5.9909530859460842 | 2.4727486265851006e-05 | 8.8740576032537707e-26 | 1.5571690359595192e-22 | PASS |
| 32 | r_K_ratio | 1.8212005047589639 | 1.9353871826449114 | 0.062698575795233608 | 0.0014984898683687194 | 41.841174317373465 | 14.438009332940839 | 3.3240794979143557e-10 | 1.9529226481982057e-25 | 3.5361982730036503e-18 | PASS |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.022421524663677129 | 2.7999999999999998 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 32 | r90 | 20.8125 | 18.3125 | 0.12012012012012012 | 0.006006006006006006 | 20 | -19.364916731037081 | 5.0333996972544363e-12 | 2.4953178164135283e-13 | 6.8701573430607085e-19 | PASS |
| 32 | r99 | 24.375 | 21.8125 | 0.10512820512820513 | 0.0051282051282051282 | 20.5 | -16.291747991900039 | 6.0154840401030353e-11 | 1.0891629608238304e-16 | 1.894984565766644e-18 | PASS |
| 32 | radial_profile_L2 | null | null | 0.21018095506491313 | 0.026253620259367155 | 8.0057894106974317 | 13.412464775332387 | 9.3140173054568233e-10 | 1.4607628044708794e-20 | 1.1676257930585241e-11 | PASS |
| 36 | total_mass | 404.5 | 388.86825616071297 | 0.038644607760907321 | 0.001854140914709518 | 20.842325119049349 | -12.225825152194394 | 3.3489968577194366e-09 | 1.7172537470951892e-25 | 2.4289733195670448e-22 | PASS |
| 36 | r_K_ratio | 1.6452597661040846 | 1.8903508514388949 | 0.14896801732116566 | 0.0025243036803777237 | 59.013508746647659 | 18.714486324847833 | 8.2432438344289491e-12 | 6.9031866547462477e-23 | 5.3183462254339121e-12 | PASS |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0089686098654708519 | 7 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 36 | r90 | 21.3125 | 18.5 | 0.13196480938416422 | 0 | null | -20.683735189961475 | 1.93756360331637e-12 | 1.2662618431714426e-12 | 8.362044513795371e-19 | PASS |
| 36 | r99 | 24.625 | 22 | 0.1065989847715736 | 0.0025380710659898475 | 42 | -16.959029914832215 | 3.3919667718492142e-11 | 1.4370664950489062e-17 | 2.2632691711460884e-18 | PASS |
| 36 | radial_profile_L2 | null | null | 0.22528336376375416 | 0.021200415960106734 | 10.626365265081334 | 16.462499317508072 | 5.1848520212984635e-11 | 1.1600148363176953e-20 | 6.9797343730044553e-11 | PASS |
| 40 | total_mass | 421.125 | 405.84907570000098 | 0.03627408560403441 | 0.0014841199168892847 | 24.441478879998385 | -12.745227200085626 | 1.8892615633971754e-09 | 1.7452483027066947e-26 | 8.1186936322662808e-23 | PASS |
| 40 | r_K_ratio | 1.4862818219349081 | 1.6983121509907673 | 0.142658226674554 | 0.0052088165665525443 | 27.387838456552206 | 16.570339045211455 | 4.7237755735314218e-11 | 1.6239908173181824e-22 | 1.4061513528213643e-11 | PASS |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0044843049327354259 | 14 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 40 | r90 | 21.5 | 18.8125 | 0.125 | 0 | null | -22.456017618285021 | 5.8544703070058476e-13 | 2.3582850886443684e-14 | 3.6683425778814834e-19 | PASS |
| 40 | r99 | 25 | 22 | 0.12 | 0.0074999999999999997 | 16 | -23.2379000772445 | 3.5513304925948273e-13 | 2.2385921293450308e-18 | 7.4091209429095752e-20 | PASS |
| 40 | radial_profile_L2 | null | null | 0.23914234348578436 | 0.019135549966259355 | 12.497280919934399 | 15.026895393764976 | 1.8927916382361989e-10 | 9.4044421932240575e-22 | 4.8080549724830386e-11 | PASS |
| 44 | total_mass | 433 | 437.05330927096486 | 0.0093609913879095421 | 0.0072170900692840644 | 1.297058966708746 | 2.9024963027503081 | 0.010939159915657146 | 5.8696913328327957e-26 | 2.1548953493993571e-21 | PASS |
| 44 | r_K_ratio | 1.3851122690599202 | 1.4774692747517539 | 0.066678353628703821 | 0.023883569861476831 | 2.7918085116853959 | 9.6651891056500698 | 7.8066528014100718e-08 | 1.9733876589946255e-22 | 4.9334668980457737e-15 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 14.0625 | 13.0625 | 0.071111111111111111 | 0.0088888888888888889 | 8 | -7.7459666924148332 | 1.2783477850604916e-06 | 1.25956553381152e-13 | 5.7635238603637946e-15 | PASS |
| 44 | r90 | 21.875 | 18.9375 | 0.13428571428571429 | 0.0057142857142857143 | 23.5 | -47 | 1.0595740345741668e-17 | 1.9641369866418272e-18 | 2.4448502950007348e-23 | PASS |
| 44 | r99 | 25.4375 | 22.375 | 0.12039312039312039 | 0.0024570024570024569 | 49 | -21.35148884619916 | 1.2208711374192415e-12 | 9.3308929176533803e-17 | 1.3651351700783063e-19 | PASS |
| 44 | radial_profile_L2 | null | null | 0.24382436372040828 | 0.022666905945975712 | 10.756843668983279 | 15.671948334281529 | 1.0440894992172389e-10 | 1.594444623288765e-22 | 1.8393228315232689e-11 | PASS |
| 48 | total_mass | 468.375 | 489.58686958999675 | 0.045288219033886856 | 0.0066720042700827327 | 6.7877982687989622 | 11.228643589852014 | 1.0658938277190579e-08 | 5.9781021731165872e-25 | 3.3492454759012593e-19 | PASS |
| 48 | r_K_ratio | 1.4233832881828432 | 1.3957945207112992 | 0.019382528726162734 | 0.023969344303723984 | 0.80863825395304523 | -2.4416599649049333 | 0.027487869696263403 | 6.4540051805283166e-20 | 7.3369478057528394e-16 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 14.4375 | 13.1875 | 0.086580086580086577 | 0.030303030303030304 | 2.8571428571428572 | -11.180339887498949 | 1.129731174656377e-08 | 5.7520127976418164e-14 | 2.7873182537487998e-16 | PASS |
| 48 | r90 | 22.375 | 19 | 0.15083798882681565 | 0 | null | -21.804467033355703 | 8.9935030956328787e-13 | 3.0419385484799993e-12 | 6.5100670518587811e-18 | PASS |
| 48 | r99 | 26.1875 | 22.9375 | 0.12410501193317422 | 0.0071599045346062056 | 17.333333333333332 | -29.068883707497267 | 1.3242683256898944e-14 | 1.1102156212733309e-18 | 1.7102375145581811e-21 | PASS |
| 48 | radial_profile_L2 | null | null | 0.27347908721153685 | 0.017773720507827693 | 15.386710232733456 | 15.382532721048291 | 1.3597767410060643e-10 | 2.9895194277448996e-21 | 4.7735071465391814e-08 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1796875, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

### ablation-spatial-birth/ablation_report.json

Equivalence by intervention: `{"full": true, "no_migration": true}`.

| Intervention | Hours | Metric | ABM effect | PDE effect | Interaction/error change | SE | Paired t | Paired p |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| no_migration | 0 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_migration | 4 | total_mass | 0 | 0.00012560772118597185 | 0.00012560772118597185 | 0.00012133461811883484 | 1.0352175095070719 | 0.31696944719393905 |
| no_migration | 4 | r_K_ratio | 0 | 9.8130440741306391e-07 | 9.8130440741306391e-07 | 9.4792670324490787e-07 | 1.0352112711393178 | 0.31697226570925685 |
| no_migration | 4 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 4 | r50 | -0.0625 | -0.0625 | 0 | 0 | null | 1 |
| no_migration | 4 | r90 | -0.5625 | -0.3125 | 0.25 | 0.11180339887498948 | 2.2360679774997898 | 0.040968955955836141 |
| no_migration | 4 | r99 | -0.3125 | -0.4375 | -0.125 | 0.125 | -1 | 0.33317013591547739 |
| no_migration | 4 | radial_profile_L2 | null | null | -0.01390051111786134 | null | null | null |
| no_migration | 8 | total_mass | 0 | 0.00012560920900384076 | 0.00012560920900384076 | 0.00012133461754712231 | 1.0352297764902778 | 0.31696390498294125 |
| no_migration | 8 | r_K_ratio | 0 | 9.8130515333721968e-07 | 9.8130515333721968e-07 | 9.4792668872207049e-07 | 1.0352120739000901 | 0.3169719030182524 |
| no_migration | 8 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 8 | r50 | -0.0625 | -0.0625 | 0 | 0 | null | 1 |
| no_migration | 8 | r90 | -0.6875 | -0.625 | 0.0625 | 0.0625 | 1 | 0.33317013591547739 |
| no_migration | 8 | r99 | -1.125 | -1.0625 | 0.0625 | 0.0625 | 1 | 0.33317013591547739 |
| no_migration | 8 | radial_profile_L2 | null | null | -0.067252396223438593 | null | null | null |
| no_migration | 12 | total_mass | 0 | 0.00012561093324059414 | 0.00012561093324059414 | 0.00012133461685530307 | 1.0352439929850421 | 0.31695748207218577 |
| no_migration | 12 | r_K_ratio | 0 | 9.8130625592746101e-07 | 9.8130625592746101e-07 | 9.4792667803676331e-07 | 1.0352132487291419 | 0.31697137222562988 |
| no_migration | 12 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 12 | r50 | -0.0625 | 0 | 0.0625 | 0.0625 | 1 | 0.33317013591547739 |
| no_migration | 12 | r90 | -1.3125 | -0.625 | 0.6875 | 0.11967838846954226 | 5.7445626465380286 | 3.8761814875216441e-05 |
| no_migration | 12 | r99 | -1.625 | -1.1875 | 0.4375 | 0.12808688457449499 | 3.415650255319866 | 0.0038327022878305922 |
| no_migration | 12 | radial_profile_L2 | null | null | -0.072893691384708734 | null | null | null |
| no_migration | 16 | total_mass | 0 | -0.00045255709874680861 | -0.00045255709874680861 | 0.00012653703696948577 | -3.5764793422177426 | 0.0027554243253045505 |
| no_migration | 16 | r_K_ratio | 0 | -3.5356453209484107e-06 | -3.5356453209484107e-06 | 9.8857058051810968e-07 | -3.576522901475967 | 0.0027551781529567632 |
| no_migration | 16 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 16 | r50 | -0.125 | 0 | 0.125 | 0.085391256382996647 | 1.4638501094227998 | 0.16387561365565406 |
| no_migration | 16 | r90 | -1.6875 | -0.625 | 1.0625 | 0.0625 | 17 | 3.2769117377966068e-11 |
| no_migration | 16 | r99 | -2.25 | -1.875 | 0.375 | 0.15478479684172258 | 2.4227185592617446 | 0.028528068449741654 |
| no_migration | 16 | radial_profile_L2 | null | null | -0.087752328239160862 | null | null | null |
| no_migration | 20 | total_mass | 0 | -0.2464161854248097 | -0.2464161854248097 | 0.0057522118254335289 | -42.838510281432129 | 4.2157438096982375e-17 |
| no_migration | 20 | r_K_ratio | 0 | -0.0019251243831664433 | -0.0019251243831664433 | 4.4939135201489528e-05 | -42.838483084619618 | 4.2157836501729852e-17 |
| no_migration | 20 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 20 | r50 | -0.125 | 0 | 0.125 | 0.085391256382996647 | 1.4638501094227998 | 0.16387561365565406 |
| no_migration | 20 | r90 | -2 | -0.75 | 1.25 | 0.11180339887498948 | 11.180339887498949 | 1.129731174656377e-08 |
| no_migration | 20 | r99 | -3 | -2.0625 | 0.9375 | 0.14343262065048754 | 6.5361700549589266 | 9.4202794798815216e-06 |
| no_migration | 20 | radial_profile_L2 | null | null | -0.091628362271557046 | null | null | null |
| no_migration | 24 | total_mass | -1.625 | -1.9044048627838492 | -0.27940486278384924 | 0.52477472118262747 | -0.53242820491464415 | 0.60222789667714349 |
| no_migration | 24 | r_K_ratio | -0.0126953125 | -0.014876560122035393 | -0.0021812476220353927 | 0.0040998096154116689 | -0.53203632037834758 | 0.60249266425193171 |
| no_migration | 24 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 24 | r50 | -0.5 | 0 | 0.5 | 0.12909944487358058 | 3.8729833462074166 | 0.0015017735386984323 |
| no_migration | 24 | r90 | -2.25 | -0.8125 | 1.4375 | 0.15728821740147395 | 9.1392732637488034 | 1.6109517538759161e-07 |
| no_migration | 24 | r99 | -3.8125 | -2.375 | 1.4375 | 0.22302372818454394 | 6.4455025108831538 | 1.1028806709730556e-05 |
| no_migration | 24 | radial_profile_L2 | null | null | -0.11830384975216335 | null | null | null |
| no_migration | 28 | total_mass | -5.8125 | -4.1506953596657148 | 1.6618046403342852 | 1.097777193034877 | 1.5137904584628119 | 0.15086290256705037 |
| no_migration | 28 | r_K_ratio | -0.04541015625 | -0.032289082658083598 | 0.013121073591916402 | 0.0085752357776827885 | 1.5301122828673983 | 0.14680254603397447 |
| no_migration | 28 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 28 | r50 | -0.75 | 0 | 0.75 | 0.11180339887498948 | 6.7082039324993694 | 7.0065583016648931e-06 |
| no_migration | 28 | r90 | -2.5625 | -1.25 | 1.3125 | 0.19830006723817989 | 6.6187572111286537 | 8.1682216137046505e-06 |
| no_migration | 28 | r99 | -3.6875 | -2.9375 | 0.75 | 0.19364916731037085 | 3.8729833462074166 | 0.0015017735386984323 |
| no_migration | 28 | radial_profile_L2 | null | null | -0.1281389136143764 | null | null | null |
| no_migration | 32 | total_mass | -8.4375 | -5.5707645019702206 | 2.8667354980297794 | 1.3277641434819116 | 2.1590698258442869 | 0.047451562897960299 |
| no_migration | 32 | r_K_ratio | -0.03243103546162622 | -0.040957820338540091 | -0.0085267848769138704 | 0.0082717381617152721 | -1.0308335092591605 | 0.31895462253661178 |
| no_migration | 32 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 32 | r50 | -0.8125 | -0.0625 | 0.75 | 0.14433756729740643 | 5.196152422706632 | 0.00010854971830166222 |
| no_migration | 32 | r90 | -2.8125 | -1.375 | 1.4375 | 0.15728821740147395 | 9.1392732637488034 | 1.6109517538759161e-07 |
| no_migration | 32 | r99 | -4.3125 | -3.0625 | 1.25 | 0.19364916731037085 | 6.4549722436790278 | 1.0848127703079008e-05 |
| no_migration | 32 | radial_profile_L2 | null | null | -0.13693182301987247 | null | null | null |
| no_migration | 36 | total_mass | -10.3125 | -6.3680669989901943 | 3.9444330010098057 | 1.2403963012577097 | 3.1799780417035399 | 0.0062133508906266162 |
| no_migration | 36 | r_K_ratio | -0.019031808192100444 | -0.034198472254671528 | -0.015166664062571084 | 0.0089352118516424217 | -1.697403969194444 | 0.11026477557069099 |
| no_migration | 36 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 36 | r50 | -0.875 | -0.0625 | 0.8125 | 0.13597640726733934 | 5.975301277099982 | 2.5437485863080221e-05 |
| no_migration | 36 | r90 | -3.3125 | -1.4375 | 1.875 | 0.1796988221070652 | 10.434125154603787 | 2.848681250543774e-08 |
| no_migration | 36 | r99 | -4.3125 | -3.5 | 0.8125 | 0.20854156260403664 | 3.8961058402670319 | 0.0014325630573340873 |
| no_migration | 36 | radial_profile_L2 | null | null | -0.14407542630053627 | null | null | null |
| no_migration | 40 | total_mass | -12.4375 | -8.4237867524657553 | 4.0137132475342447 | 1.1177352054277183 | 3.5909339063882637 | 0.0026749354760721091 |
| no_migration | 40 | r_K_ratio | -0.014536219382929247 | -0.022624410445823889 | -0.0080881910628946424 | 0.0086124813941057376 | -0.9391243583330221 | 0.36254442047530833 |
| no_migration | 40 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 40 | r50 | -0.875 | -0.0625 | 0.8125 | 0.13597640726733934 | 5.975301277099982 | 2.5437485863080221e-05 |
| no_migration | 40 | r90 | -3.4375 | -1.4375 | 2 | 0.20412414523193151 | 9.7979589711327115 | 6.5318351445822073e-08 |
| no_migration | 40 | r99 | -4.8125 | -3.9375 | 0.875 | 0.20155644370746373 | 4.3412157106222962 | 0.00058157801081101909 |
| no_migration | 40 | radial_profile_L2 | null | null | -0.15309848732055903 | null | null | null |
| no_migration | 44 | total_mass | -14.125 | -13.948353768478906 | 0.17664623152109371 | 1.4220855101137102 | 0.12421632191932611 | 0.90279330055230234 |
| no_migration | 44 | r_K_ratio | -0.020820605221394506 | -0.030432474283540067 | -0.0096118690621455616 | 0.0092023289176463886 | -1.0445039672200631 | 0.31279387333665098 |
| no_migration | 44 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 44 | r50 | -1 | -0.0625 | 0.9375 | 0.17001838135919303 | 5.5141096657035584 | 5.9461335275977702e-05 |
| no_migration | 44 | r90 | -3.8125 | -1.5 | 2.3125 | 0.1505199322349037 | 15.363413772941893 | 1.3839300528290396e-10 |
| no_migration | 44 | r99 | -5.1875 | -3.9375 | 1.25 | 0.19364916731037085 | 6.4549722436790278 | 1.0848127703079008e-05 |
| no_migration | 44 | radial_profile_L2 | null | null | -0.16002028372515748 | null | null | null |
| no_migration | 48 | total_mass | -25.4375 | -27.985431308373265 | -2.5479313083732649 | 1.443552635937482 | -1.7650421917026737 | 0.097901039872220361 |
| no_migration | 48 | r_K_ratio | -0.064987486008190057 | -0.082833326733961557 | -0.0178458407257715 | 0.01001767176744763 | -1.7814359603756895 | 0.095094387248461948 |
| no_migration | 48 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 48 | r50 | -1.3125 | -0.1875 | 1.125 | 0.085391256382996647 | 13.174650984805199 | 1.1942408086569779e-09 |
| no_migration | 48 | r90 | -4.0625 | -1.875 | 2.1875 | 0.20854156260403664 | 10.489515723795854 | 2.655026371957866e-08 |
| no_migration | 48 | r99 | -5.5 | -4.0625 | 1.4375 | 0.12808688457449499 | 11.22285083890813 | 1.0733435740340221e-08 |
| no_migration | 48 | radial_profile_L2 | null | null | -0.18520397585495668 | null | null | null |

#### full: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 255.99987439227883 | 4.9065516088270256e-07 | 0 | null | -1.0352175095070719 | 0.31696944719393905 | 6.3310518845588904e-81 | 6.3307856304481416e-81 | PASS |
| 4 | r_K_ratio | 1 | 0.9999990186955926 | 9.8130440741306391e-07 | 0 | null | -1.0352112711393178 | 0.31697226570925685 | 2.0746026777977294e-76 | 2.0744281865575098e-76 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13.0625 | 13.0625 | 0 | 0.0047846889952153108 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 4 | r90 | 17.9375 | 17.6875 | 0.013937282229965157 | 0.0069686411149825784 | 2 | -2.2360679774997898 | 0.040968955955836141 | 1.1148111900685798e-16 | 3.8956480726576652e-17 | PASS |
| 4 | r99 | 19.25 | 19.375 | 0.0064935064935064939 | 0.00974025974025974 | 0.66666666666666663 | 1 | 0.33317013591547739 | 1.2617854706008803e-18 | 7.1507823999893397e-17 | PASS |
| 4 | radial_profile_L2 | null | null | 0.01390051111786134 | 0.026445474852675486 | 0.52562909894034393 | -7.7751682619068863 | 1.2210780923911863e-06 | 1.7063425300474121e-19 | 5.5915970512964533e-19 | PASS |
| 8 | total_mass | 256 | 255.99987439079101 | 4.9066097267125297e-07 | 0 | null | -1.0352297764902778 | 0.31696390498294125 | 6.331051438669818e-81 | 6.3307851814239781e-81 | PASS |
| 8 | r_K_ratio | 1 | 0.99999901869484664 | 9.8130515333721968e-07 | 0 | null | -1.0352120739000901 | 0.3169719030182524 | 2.0746022011006527e-76 | 2.0744277097679295e-76 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 8 | r90 | 18.0625 | 18 | 0.0034602076124567475 | 0.0034602076124567475 | 1 | -1 | 0.33317013591547739 | 1.4301148022486204e-22 | 1.9701065876166452e-19 | PASS |
| 8 | r99 | 20.0625 | 20 | 0.0031152647975077881 | 0.0031152647975077881 | 1 | -1 | 0.33317013591547739 | 6.5149803197211003e-25 | 5.0766352092358305e-21 | PASS |
| 8 | radial_profile_L2 | null | null | 0.067252396223438593 | 0.031674679515173401 | 2.1232226261743889 | -1.1211565613575198 | 0.2798500299122274 | 1.1233733444798528e-17 | 3.6201109868020346e-15 | PASS |
| 12 | total_mass | 256 | 255.99987438921454 | 4.9066713065509804e-07 | 0 | null | -1.0352427723769593 | 0.31695803353050223 | 6.3310511325783142e-81 | 6.3307848720036545e-81 | PASS |
| 12 | r_K_ratio | 1 | 0.99999901869489871 | 9.8130510128163762e-07 | 0 | null | -1.0352120281216213 | 0.31697192370116678 | 2.0746019264410209e-76 | 2.0744274351406548e-76 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 12 | r90 | 18.6875 | 18 | 0.03678929765886288 | 0 | null | -5.7445626465380286 | 3.8761814875216441e-05 | 1.2411344634416193e-17 | 3.0051609001276733e-16 | PASS |
| 12 | r99 | 20.5625 | 20.125 | 0.021276595744680851 | 0.015197568389057751 | 1.3999999999999999 | -3.415650255319866 | 0.0038327022878305922 | 3.2453744019125243e-19 | 3.5829058514479523e-17 | PASS |
| 12 | radial_profile_L2 | null | null | 0.072893691384743886 | 0.033001760468462815 | 2.2087819058744653 | -0.25736842332884102 | 0.80039152143320558 | 4.3685958447070436e-20 | 2.3969087135166486e-17 | PASS |
| 16 | total_mass | 256 | 256.00565990929204 | 2.2109020672067548e-05 | 0 | null | 29.055621071987371 | 1.333214515984822e-14 | 7.673154674705337e-78 | 7.6877095417039995e-78 | PASS |
| 16 | r_K_ratio | 1 | 1.0000442180838114 | 4.4218083811317643e-05 | 0 | null | 29.055649465249196 | 1.3331952949265112e-14 | 2.5119576257366105e-73 | 2.5214963152902277e-73 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.0047619047619047623 | 1 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 16 | r90 | 19.0625 | 18 | 0.055737704918032788 | 0.0032786885245901639 | 17 | -17 | 3.2769117377966068e-11 | 2.2623006381054966e-21 | 5.3331248280097439e-21 | PASS |
| 16 | r99 | 21.1875 | 20.8125 | 0.017699115044247787 | 0.014749262536873156 | 1.2 | -2.4227185592617446 | 0.028528068449741654 | 4.6950640203719188e-18 | 3.0093626661365215e-16 | PASS |
| 16 | radial_profile_L2 | null | null | 0.087753521628794309 | 0.029422120612765671 | 2.9825695701457975 | 0.2075579969651119 | 0.8383657119697715 | 5.9175067186833957e-20 | 1.2294627722773073e-16 | PASS |
| 20 | total_mass | 256 | 259.38935697918242 | 0.013239675699931355 | 0 | null | 48.255478718939976 | 7.1525658805689808e-18 | 9.9621099804360735e-40 | 3.100505078257499e-39 | PASS |
| 20 | r_K_ratio | 1 | 1.0264793310107851 | 0.026479331010785034 | 0 | null | 48.255459289196338 | 7.1526088194555449e-18 | 1.9076976227082329e-35 | 1.8537942546079051e-34 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0 | null | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 20 | r90 | 19.375 | 18.125 | 0.064516129032258063 | 0 | null | -11.180339887498949 | 1.129731174656377e-08 | 6.9886341424657286e-17 | 1.1645945213682297e-17 | PASS |
| 20 | r99 | 21.9375 | 21 | 0.042735042735042736 | 0 | null | -6.5361700549589266 | 9.4202794798815216e-06 | 3.6683425778815099e-19 | 5.1861867803112176e-17 | PASS |
| 20 | radial_profile_L2 | null | null | 0.09226338129301527 | 0.040888641583758337 | 2.256455037862247 | 1.0231584195297536 | 0.32245164838009688 | 1.1251404688142578e-19 | 3.5060034693780733e-16 | PASS |
| 24 | total_mass | 284.625 | 298.98832661654501 | 0.050464037300114242 | 0.0083443126921387799 | 6.0477164701242172 | 19.731155400039221 | 3.8383869057078653e-12 | 1.6655364014943582e-27 | 4.0388732636059276e-22 | PASS |
| 24 | r_K_ratio | 1.2236328125 | 1.3358311402728638 | 0.091692807373832841 | 0.015163607342378291 | 6.0468993494427652 | 19.728534634782985 | 3.8457748046898245e-12 | 2.9673416993492329e-24 | 2.7676270594412945e-17 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.5625 | 13.0625 | 0.036866359447004608 | 0.018433179723502304 | 2 | -3.8729833462074166 | 0.0015017735386984323 | 9.6072681659644597e-15 | 7.6467248765302002e-14 | PASS |
| 24 | r90 | 19.9375 | 18.375 | 0.078369905956112859 | 0.003134796238244514 | 25 | -12.198750911856663 | 3.4523439761073071e-09 | 1.4966769468733605e-14 | 2.8116770770607522e-18 | PASS |
| 24 | r99 | 22.875 | 21.375 | 0.065573770491803282 | 0.0054644808743169399 | 12 | -6.7082039324993694 | 7.0065583016648931e-06 | 3.4512116296941117e-15 | 3.5841850263120592e-15 | PASS |
| 24 | radial_profile_L2 | null | null | 0.12927642609054224 | 0.031998880089473554 | 4.0400297050730032 | 3.9655118418184734 | 0.0012436052777004001 | 3.3218560868433396e-19 | 3.309702997651699e-14 | PASS |
| 28 | total_mass | 360.0625 | 355.79786995269785 | 0.011844138301828507 | 0.0034716195105016492 | 3.4117040378417016 | -4.9934353141941079 | 0.00016039731723934463 | 4.9448840042546594e-27 | 1.120423960280407e-23 | PASS |
| 28 | r_K_ratio | 1.81298828125 | 1.7784099485173064 | 0.019072562735404398 | 0.0053864799353622404 | 3.5408212718278262 | -5.1851935898516457 | 0.00011085077445508441 | 4.8490283073764659e-24 | 5.9749820292080538e-21 | PASS |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.875 | 13.0625 | 0.058558558558558557 | 0.022522522522522521 | 2.6000000000000001 | -8.0622577482985509 | 7.8258249630702749e-07 | 2.2778649613890186e-15 | 2.5893641111623486e-16 | PASS |
| 28 | r90 | 20.5625 | 18.8125 | 0.085106382978723402 | 0.0060790273556231003 | 14 | -12.124355652982143 | 3.7541255191393749e-09 | 1.4192741870443792e-14 | 4.4765474968068944e-17 | PASS |
| 28 | r99 | 23.625 | 21.9375 | 0.071428571428571425 | 0 | null | -8.5098307000225546 | 3.9910679004152207e-07 | 1.2402986697836266e-16 | 5.9545457575235637e-16 | PASS |
| 28 | radial_profile_L2 | null | null | 0.16405669091957067 | 0.025528598153441336 | 6.4263885519093931 | 7.831324694430629 | 1.1183966370800388e-06 | 1.966819189326678e-20 | 7.1409093913155926e-14 | PASS |
| 32 | total_mass | 385.5625 | 379.16568384801553 | 0.016590866985208501 | 0.0024315124007132437 | 6.8232705621167495 | -5.3686155293505395 | 7.8187578527019736e-05 | 1.5288902634240917e-25 | 2.2786753490476405e-22 | PASS |
| 32 | r_K_ratio | 1.8212005047589639 | 1.9393568638136029 | 0.064878281521384168 | 0.0014984898683687194 | 43.295775894709074 | 14.741142630664925 | 2.4814496209643969e-10 | 2.730946823123324e-25 | 4.471149443760827e-18 | PASS |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.022421524663677129 | 2.7999999999999998 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 32 | r90 | 20.8125 | 18.9375 | 0.090090090090090086 | 0.006006006006006006 | 15 | -21.957751641341996 | 8.1207854767429677e-13 | 2.5303126792687863e-18 | 2.9422487852212645e-20 | PASS |
| 32 | r99 | 24.375 | 22.0625 | 0.094871794871794868 | 0.0051282051282051282 | 18.5 | -13.136324646189436 | 1.2434982186844375e-09 | 1.5609213681963434e-16 | 1.554563354209922e-17 | PASS |
| 32 | radial_profile_L2 | null | null | 0.17487892574276687 | 0.026253620259367155 | 6.6611356458685345 | 10.397644941576297 | 2.9843462695827865e-08 | 2.6100874687186666e-21 | 3.2227504383943266e-14 | PASS |
| 36 | total_mass | 404.5 | 389.45926766739342 | 0.037183516273440319 | 0.001854140914709518 | 20.054309776808811 | -11.805436910281133 | 5.4031100015012418e-09 | 1.4630608635349803e-25 | 2.4759221327288636e-22 | PASS |
| 36 | r_K_ratio | 1.6452597661040846 | 1.893289657624615 | 0.15075424357325423 | 0.0025243036803777237 | 59.721120222229423 | 18.645946032483344 | 8.6910660132474326e-12 | 9.5829111228214293e-23 | 7.1031951951600809e-12 | PASS |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0089686098654708519 | 7 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 36 | r90 | 21.3125 | 19 | 0.10850439882697947 | 0 | null | -19.322619811082461 | 5.195040620846464e-12 | 7.7864380968352111e-16 | 1.5335258890741899e-18 | PASS |
| 36 | r99 | 24.625 | 22.5 | 0.086294416243654817 | 0.0025380710659898475 | 34 | -13.728738502483221 | 6.730952139452569e-10 | 4.5587789692221704e-17 | 1.8183189209815112e-18 | PASS |
| 36 | radial_profile_L2 | null | null | 0.18656080889003238 | 0.021200415960106734 | 8.7998654951434805 | 12.425650837851045 | 2.6806419917659458e-09 | 2.2767547083138812e-21 | 1.0718923015478952e-13 | PASS |
| 40 | total_mass | 421.125 | 406.64803311283043 | 0.034376887829432067 | 0.0014841199168892847 | 23.163147019471324 | -12.095655968352325 | 3.8779148268028664e-09 | 1.473597627386941e-26 | 8.6918298750280041e-23 | PASS |
| 40 | r_K_ratio | 1.4862818219349081 | 1.6996050002855141 | 0.14352808141923748 | 0.0052088165665525443 | 27.554835073455379 | 16.517311941195551 | 4.9447938999950111e-11 | 1.9915408596028513e-22 | 1.6403504379960735e-11 | PASS |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 13.9375 | 13.0625 | 0.062780269058295965 | 0.0044843049327354259 | 14 | -10.2469507659596 | 3.6215788445874798e-08 | 5.2040873952684304e-16 | 9.9200721974583402e-18 | PASS |
| 40 | r90 | 21.5 | 19.0625 | 0.11337209302325581 | 0 | null | -19.030051422496395 | 6.4760803247711866e-12 | 6.3605411260372649e-15 | 2.3244786420838592e-18 | PASS |
| 40 | r99 | 25 | 22.9375 | 0.082500000000000004 | 0.0074999999999999997 | 11 | -14.379574120909639 | 3.5189622604582387e-10 | 2.4558314713141745e-18 | 7.901250573376388e-19 | PASS |
| 40 | radial_profile_L2 | null | null | 0.19809400512033434 | 0.019135549966259355 | 10.352145899627782 | 11.926715870016167 | 4.7000113633628995e-09 | 2.666940545963059e-22 | 5.2419008818329114e-14 | PASS |
| 44 | total_mass | 433 | 439.23206584062973 | 0.014392761756650638 | 0.0072170900692840644 | 1.9942610690015123 | 4.3771790308701704 | 0.00054109250360977396 | 5.9224939034279197e-26 | 3.517435842868152e-21 | PASS |
| 44 | r_K_ratio | 1.3851122690599202 | 1.4831263019519321 | 0.070762518736863317 | 0.023883569861476831 | 2.9628116377610789 | 10.19887506362087 | 3.8539975473241908e-08 | 1.7147984527185799e-22 | 6.6649967260819194e-15 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 14.0625 | 13.0625 | 0.071111111111111111 | 0.0088888888888888889 | 8 | -7.7459666924148332 | 1.2783477850604916e-06 | 1.25956553381152e-13 | 5.7635238603637946e-15 | PASS |
| 44 | r90 | 21.875 | 19.1875 | 0.12285714285714286 | 0.0057142857142857143 | 21.5 | -22.456017618285021 | 5.8544703070058476e-13 | 4.6040517988789706e-14 | 9.7976516786557512e-20 | PASS |
| 44 | r99 | 25.4375 | 23 | 0.095823095823095825 | 0.0024570024570024569 | 39 | -15.497028577661013 | 1.2241974350227701e-10 | 5.016695892938304e-18 | 2.6446581536787945e-18 | PASS |
| 44 | radial_profile_L2 | null | null | 0.2002259127025463 | 0.022666905945975712 | 8.8334028993531195 | 13.096159457147015 | 1.2974403191799644e-09 | 4.0060695047064146e-23 | 1.0619373488247218e-14 | PASS |
| 48 | total_mass | 468.375 | 496.1627205009305 | 0.059327932748183647 | 0.0066720042700827327 | 8.8920705602977641 | 14.205603014782996 | 4.1744770622290889e-10 | 5.9167774317339044e-25 | 1.0508057359812269e-18 | PASS |
| 48 | r_K_ratio | 1.4233832881828432 | 1.419135061376732 | 0.0029845979233988911 | 0.023969344303723984 | 0.12451729532439444 | -0.36683445070060317 | 0.71886582056152537 | 5.4701841780063292e-20 | 1.8578973401267354e-15 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 14.4375 | 13.375 | 0.073593073593073599 | 0.030303030303030304 | 2.4285714285714284 | -17 | 3.2769117377966068e-11 | 3.3370933113405186e-17 | 1.4553451008790879e-19 | PASS |
| 48 | r90 | 22.375 | 19.8125 | 0.11452513966480447 | 0 | null | -20.005951495444929 | 3.141619345032701e-12 | 4.6042541344059545e-15 | 1.5128771287691139e-18 | PASS |
| 48 | r99 | 26.1875 | 23.375 | 0.10739856801909307 | 0.0071599045346062056 | 15 | -17.172737481873973 | 2.8355402110316446e-11 | 5.1788953735629637e-16 | 2.4582963915343147e-19 | PASS |
| 48 | radial_profile_L2 | null | null | 0.22627052907292156 | 0.017773720507827693 | 12.73062266132013 | 12.215444655901129 | 3.3882278396906931e-09 | 3.3668681730329636e-21 | 2.472721338474e-11 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2109375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1875, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

#### no_migration: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 0 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 4 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 4 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 4 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 8 | total_mass | 256 | 256 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 8 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 8 | radial_profile_L2 | null | null | 0 | 0.036932043236795692 | 0 | -13.759882903658182 | 6.5213491597846189e-10 | 0 | 0 | PASS |
| 12 | total_mass | 256 | 256.00000000014779 | 5.773159728050814e-13 | 0 | null | 39.511502884993675 | 1.4033866702015765e-16 | 1.3666133886517828e-193 | 1.366613388719367e-193 | PASS |
| 12 | r_K_ratio | 1 | 1.0000000000011546 | 1.1546458233979706e-12 | 0 | null | 39.500857025097488 | 1.4090197472796449e-16 | 4.4970670758744223e-189 | 4.4970670763194719e-189 | PASS |
| 12 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 12 | r50 | 13 | 13.0625 | 0.004807692307692308 | 0 | null | 1 | 0.33317013591547739 | 8.8277683744579181e-19 | 1.5672359974554864e-18 | PASS |
| 12 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 12 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 12 | radial_profile_L2 | null | null | 3.5152776611498511e-14 | 0.036932043236795692 | 9.5182322803292711e-13 | -13.759882903649116 | 6.5213491598445885e-10 | 2.2876041024113787e-194 | 2.2876041024182707e-194 | PASS |
| 16 | total_mass | 256 | 256.00520735219328 | 2.0341219505087826e-05 | 0 | null | 41.038918914450051 | 7.9840648437680403e-17 | 1.2377497124113274e-80 | 1.239909652423626e-80 | PASS |
| 16 | r_K_ratio | 1 | 1.0000406824384904 | 4.0682438490369233e-05 | 0 | null | 41.038918693047371 | 7.9840654845035165e-17 | 4.0523238935208809e-76 | 4.0664792811293981e-76 | PASS |
| 16 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 16 | r50 | 13 | 13.0625 | 0.004807692307692308 | 0 | null | 1 | 0.33317013591547739 | 8.8277683744579181e-19 | 1.5672359974554864e-18 | PASS |
| 16 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 16 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 16 | radial_profile_L2 | null | null | 1.1933896334477404e-06 | 0.036932043236795692 | 3.2313122396075738e-05 | -13.759647933014801 | 6.5229043169828046e-10 | 1.9559990456430911e-81 | 1.9561991360751853e-81 | PASS |
| 20 | total_mass | 256 | 259.14294079375759 | 0.012277112475615692 | 0 | null | 48.165992154089778 | 7.3532724300682568e-18 | 3.4357029280173264e-40 | 9.8450568991714422e-40 | PASS |
| 20 | r_K_ratio | 1 | 1.0245542066276185 | 0.024554206627618591 | 0 | null | 48.165974271291155 | 7.3533131338035998e-18 | 6.8280069968818597e-36 | 5.6209090932848839e-35 | PASS |
| 20 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 20 | r50 | 13 | 13.0625 | 0.004807692307692308 | 0 | null | 1 | 0.33317013591547739 | 8.8277683744579181e-19 | 1.5672359974554864e-18 | PASS |
| 20 | r90 | 17.375 | 17.375 | 0 | 0.010791366906474821 | 0 | null | 1 | 4.7721586994314007e-25 | 4.7721586994314007e-25 | PASS |
| 20 | r99 | 18.9375 | 18.9375 | 0 | 0.0033003300330033004 | 0 | null | 1 | 4.0195340979314922e-30 | 4.0195340979314922e-30 | PASS |
| 20 | radial_profile_L2 | null | null | 0.00063501902145822139 | 0.036932043236795692 | 0.017194256418110839 | -13.65494739268801 | 7.2566506696869955e-10 | 1.0635224662041332e-40 | 1.1230143065589982e-40 | PASS |
| 24 | total_mass | 283 | 297.08392175376116 | 0.049766507963820379 | 0.0079505300353356883 | 6.2595207794494074 | 19.288997068272053 | 5.3274717602656704e-12 | 1.9667852482372862e-27 | 4.1901159709178463e-22 | PASS |
| 24 | r_K_ratio | 1.2109375 | 1.3209545801508285 | 0.090852814576167992 | 0.014516129032258065 | 6.2587494485804616 | 19.286745565529355 | 5.336467248195099e-12 | 3.7684159615241649e-24 | 3.0366844378662579e-17 | PASS |
| 24 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 24 | r50 | 13.0625 | 13.0625 | 0 | 0.0095693779904306216 | 0 | null | 1 | 1.054554752596167e-27 | 1.054554752596167e-27 | PASS |
| 24 | r90 | 17.6875 | 17.5625 | 0.0070671378091872791 | 0.0070671378091872791 | 1 | -1 | 0.33317013591547739 | 1.5298594698288029e-16 | 1.4420322929484584e-15 | PASS |
| 24 | r99 | 19.0625 | 19 | 0.0032786885245901639 | 0.0065573770491803279 | 0.5 | -1 | 0.33317013591547739 | 1.4135396879852e-24 | 1.0821089527179208e-20 | PASS |
| 24 | radial_profile_L2 | null | null | 0.010972576338378889 | 0.031004829841119182 | 0.35389893750769286 | -13.996735745454139 | 5.1371459942190642e-10 | 2.1805058859519743e-19 | 5.5634252767773985e-19 | PASS |
| 28 | total_mass | 354.25 | 351.64717459303216 | 0.0073474252843128912 | 0.0061750176429075515 | 1.189863043185299 | -2.2739001532793344 | 0.038095494685434876 | 2.1973093362007688e-25 | 1.2290504817305744e-21 | PASS |
| 28 | r_K_ratio | 1.767578125 | 1.7461208658592229 | 0.012139355447599927 | 0.0096685082872928173 | 1.2555561920089069 | -2.4002894987255758 | 0.029808058154546454 | 2.2115691668041131e-22 | 8.1955444543162061e-19 | PASS |
| 28 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 28 | r50 | 13.125 | 13.0625 | 0.0047619047619047623 | 0.014285714285714285 | 0.33333333333333331 | -1 | 0.33317013591547739 | 5.4301107888314235e-20 | 2.3575928785707159e-17 | PASS |
| 28 | r90 | 18 | 17.5625 | 0.024305555555555556 | 0 | null | -3.415650255319866 | 0.0038327022878305922 | 1.8264536248980875e-15 | 1.0127157816886049e-16 | PASS |
| 28 | r99 | 19.9375 | 19 | 0.047021943573667714 | 0.006269592476489028 | 7.5 | -15 | 1.9412759304364392e-10 | 7.8745481414756495e-24 | 7.3606614810643734e-22 | PASS |
| 28 | radial_profile_L2 | null | null | 0.035917777305194278 | 0.035579270731715518 | 1.0095141515415327 | -5.4958555805698035 | 6.1529877947217449e-05 | 3.8461695488513462e-21 | 8.3710689501974908e-20 | PASS |
| 32 | total_mass | 377.125 | 373.59491934604534 | 0.0093605055457863396 | 0.0041431885979449782 | 2.2592516185309908 | -3.224394776907141 | 0.0056732368974357402 | 6.98687355826189e-26 | 1.330894245467728e-22 | PASS |
| 32 | r_K_ratio | 1.7887694692973377 | 1.8983990434750628 | 0.06128770423434704 | 0.0027342015763558028 | 22.415210628337242 | 10.53418487564759 | 2.5090134635364904e-08 | 6.2883214813276903e-23 | 1.5142085334053166e-16 | PASS |
| 32 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 32 | r50 | 13.125 | 13 | 0.0095238095238095247 | 0.0095238095238095247 | 1 | -1.4638501094227998 | 0.16387561365565406 | 2.6121019521635124e-18 | 1.6679850051069694e-15 | PASS |
| 32 | r90 | 18 | 17.5625 | 0.024305555555555556 | 0 | null | -3.415650255319866 | 0.0038327022878305922 | 1.8264536248980875e-15 | 1.0127157816886049e-16 | PASS |
| 32 | r99 | 20.0625 | 19 | 0.052959501557632398 | 0 | null | -17 | 3.2769117377966068e-11 | 1.0235469112862909e-23 | 5.1995965524563549e-22 | PASS |
| 32 | radial_profile_L2 | null | null | 0.037947102722894387 | 0.026598751298356628 | 1.4266497813090535 | -4.5328240799512107 | 0.00039651319405034364 | 5.659590274262924e-20 | 1.4614894831104819e-18 | PASS |
| 36 | total_mass | 394.1875 | 383.09120066840319 | 0.028149800111867584 | 0.0019026478515934676 | 14.795065775462405 | -9.1221676453946667 | 1.6501877724278708e-07 | 3.7980178864944126e-25 | 1.4213484823186936e-22 | PASS |
| 36 | r_K_ratio | 1.6262279579119843 | 1.8590911853699434 | 0.14319224209928538 | 0.0062477019091505366 | 22.919185995343749 | 14.184787769910242 | 4.2611946814913291e-10 | 5.8175451191989562e-21 | 8.255757762802596e-11 | PASS |
| 36 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 36 | r50 | 13.0625 | 13 | 0.0047846889952153108 | 0 | null | -1 | 0.33317013591547739 | 1.981502770015363e-20 | 2.2905328922157176e-17 | PASS |
| 36 | r90 | 18 | 17.5625 | 0.024305555555555556 | 0 | null | -3.415650255319866 | 0.0038327022878305922 | 1.8264536248980875e-15 | 1.0127157816886049e-16 | PASS |
| 36 | r99 | 20.3125 | 19 | 0.064615384615384616 | 0 | null | -10.966892325208963 | 1.464381425614242e-08 | 2.9229113868319653e-19 | 4.3802709106868439e-18 | PASS |
| 36 | radial_profile_L2 | null | null | 0.042485382589496093 | 0.022848576447345625 | 1.859432367150021 | -4.2397745933582529 | 0.0007132842061844133 | 1.5558616830400435e-19 | 5.932271884556231e-18 | PASS |
| 40 | total_mass | 408.6875 | 398.22424636036465 | 0.025602088734388337 | 0.0013763572411683743 | 18.601339803796147 | -8.1562442622931073 | 6.7802220509233907e-07 | 6.4951744170669196e-25 | 1.8396975596743288e-22 | PASS |
| 40 | r_K_ratio | 1.4717456025519788 | 1.6769805898396903 | 0.13945004281435452 | 0.0092852152085018697 | 15.018504114656295 | 12.439058335072291 | 2.6411842787440857e-09 | 2.0577233729962925e-20 | 3.5027054837761317e-10 | PASS |
| 40 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 40 | r50 | 13.0625 | 13 | 0.0047846889952153108 | 0 | null | -1 | 0.33317013591547739 | 1.981502770015363e-20 | 2.2905328922157176e-17 | PASS |
| 40 | r90 | 18.0625 | 17.625 | 0.024221453287197232 | 0.0034602076124567475 | 7 | -3.415650255319866 | 0.0038327022878305922 | 1.1216947122086982e-15 | 1.7869936330514343e-16 | PASS |
| 40 | r99 | 20.1875 | 19 | 0.058823529411764705 | 0.0030959752321981426 | 19 | -11.783299785974805 | 5.5430657455111653e-09 | 1.7075994239206279e-20 | 4.6936374131487617e-19 | PASS |
| 40 | radial_profile_L2 | null | null | 0.044995517799775314 | 0.024227567023808398 | 1.8572033153621359 | -2.8113863386119946 | 0.013153945947298998 | 1.011128869248533e-19 | 4.795637554265591e-18 | PASS |
| 44 | total_mass | 418.875 | 425.28371207215082 | 0.015299819927545975 | 0.0032826022082960309 | 4.6608815070187779 | 3.8407257093995542 | 0.0016040201076298938 | 4.9385096470234179e-24 | 3.0861586791854761e-20 | PASS |
| 44 | r_K_ratio | 1.3642916638385258 | 1.4526938276683921 | 0.06479711499602725 | 0.02354061237077195 | 2.7525670944940845 | 5.6927671337940957 | 4.2648139970186079e-05 | 8.3727463737289795e-19 | 3.9033998965655983e-12 | PASS |
| 44 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 44 | r50 | 13.0625 | 13 | 0.0047846889952153108 | 0 | null | -1 | 0.33317013591547739 | 1.981502770015363e-20 | 2.2905328922157176e-17 | PASS |
| 44 | r90 | 18.0625 | 17.6875 | 0.020761245674740483 | 0.0034602076124567475 | 6 | -3 | 0.0089727374772233335 | 5.712038532266603e-16 | 1.6369019078997911e-16 | PASS |
| 44 | r99 | 20.25 | 19.0625 | 0.058641975308641972 | 0.0030864197530864196 | 19 | -11.783299785974805 | 5.5430657455111653e-09 | 4.1986501502418116e-20 | 4.033425490932737e-19 | PASS |
| 44 | radial_profile_L2 | null | null | 0.040205628977388816 | 0.021099978364348116 | 1.9054820001769688 | -3.4652566370477 | 0.0034617453379834689 | 1.5743510783731816e-20 | 4.9524854910883851e-19 | PASS |
| 48 | total_mass | 442.9375 | 468.17728919255723 | 0.056982732761523353 | 0.0046564131508395655 | 12.237473547906546 | 13.417851554938149 | 9.2621268237133873e-10 | 2.1407427937622949e-24 | 8.4525539979633163e-19 | PASS |
| 48 | r_K_ratio | 1.3583958021746529 | 1.3363017346427704 | 0.016264823180778647 | 0.023109012801366548 | 0.70383028996447783 | -1.3294936184223325 | 0.20355146261925461 | 4.7568875420667918e-17 | 3.3936195746172337e-13 | PASS |
| 48 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 48 | r50 | 13.125 | 13.1875 | 0.0047619047619047623 | 0.0047619047619047623 | 1 | 1 | 0.33317013591547739 | 1.104347949798057e-18 | 4.7540344799007721e-18 | PASS |
| 48 | r90 | 18.3125 | 17.9375 | 0.020477815699658702 | 0.010238907849829351 | 2 | -3 | 0.0089727374772233335 | 2.8181580506755761e-17 | 1.1954952800484915e-15 | PASS |
| 48 | r99 | 20.6875 | 19.3125 | 0.066465256797583083 | 0.0030211480362537764 | 22 | -11 | 1.4062516106729139e-08 | 1.5711298046302084e-17 | 1.3231466562902307e-18 | PASS |
| 48 | radial_profile_L2 | null | null | 0.04106655321796486 | 0.023214227865607628 | 1.7690251623146094 | -2.6358345861812764 | 0.01871335585702447 | 7.38288132447601e-21 | 2.5051409077953684e-19 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1640625, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1640625, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.15625, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

### ablation-active-r20/ablation_report.json

Equivalence by intervention: `{"enlarged_boundary": true, "full": true, "no_activation": true, "no_division": true, "no_exchange": true, "no_migration": true}`.

| Intervention | Hours | Metric | ABM effect | PDE effect | Interaction/error change | SE | Paired t | Paired p |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| no_migration | 0 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 0 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_migration | 4 | total_mass | -9 | -8.0308750880559998 | 0.96912491194400019 | 0.82640988171043861 | 1.172692792513784 | 0.2592066110160216 |
| no_migration | 4 | r_K_ratio | -0.005181147506914327 | -0.0027691246447194978 | 0.0024120228621948292 | 0.0016142820262789396 | 1.4941768680623617 | 0.15586648277510518 |
| no_migration | 4 | active_fraction | -0.0059052929338261461 | -0.0045556296660311541 | 0.001349663267794992 | 0.0010825570236751428 | 1.2467364196788984 | 0.23161198605079591 |
| no_migration | 4 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_migration | 4 | r90 | -1.4375 | -2 | -0.5625 | 0.12808688457449499 | -4.3915503282683988 | 0.00052573078993947925 |
| no_migration | 4 | r99 | -23.1875 | -21.5625 | 1.625 | 0.375 | 4.333333333333333 | 0.00059085778625744573 |
| no_migration | 4 | radial_profile_L2 | null | null | -0.040240558825754238 | null | null | null |
| no_migration | 8 | total_mass | -23.125 | -19.441339972849221 | 3.6836600271507791 | 1.0247537960294744 | 3.5946780987038456 | 0.0026544740931767423 |
| no_migration | 8 | r_K_ratio | -0.019576008242981807 | -0.013402938120085299 | 0.0061730701228965082 | 0.0029270484785604831 | 2.1089743364730373 | 0.052166334198908225 |
| no_migration | 8 | active_fraction | -0.017830384381605589 | -0.01242716671741766 | 0.0054032176641879287 | 0.0014554205284030592 | 3.7124786676717672 | 0.0020851149330711779 |
| no_migration | 8 | r50 | -0.125 | -0.9375 | -0.8125 | 0.10077822185373186 | -8.0622577482985509 | 7.8258249630702749e-07 |
| no_migration | 8 | r90 | -5 | -3.5 | 1.5 | 0.57008771254956903 | 2.6311740579210876 | 0.018888226923741501 |
| no_migration | 8 | r99 | -44.5 | -41.9375 | 2.5625 | 0.59839194234771131 | 4.2823103365101671 | 0.00065469959890511616 |
| no_migration | 8 | radial_profile_L2 | null | null | -0.047337996492085838 | null | null | null |
| no_activation | 0 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 0 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 0 | active_fraction | -1 | -1 | 0 | 0 | null | 1 |
| no_activation | 0 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 0 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 0 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 0 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_activation | 4 | total_mass | -8.3125 | -6.0381114111597469 | 2.2743885888402531 | 0.94869781375227158 | 2.3973793929645884 | 0.029978064013030822 |
| no_activation | 4 | r_K_ratio | -0.0047734532920918407 | -0.002113726589565823 | 0.0026597267025260177 | 0.001565998101050399 | 1.6984226869381234 | 0.1100688975685177 |
| no_activation | 4 | active_fraction | -0.30956904900325755 | -0.30629895168952725 | 0.0032700973137303226 | 0.0012426446652927233 | 2.63156267037931 | 0.018873585616829489 |
| no_activation | 4 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_activation | 4 | r90 | -1.0625 | -1 | 0.0625 | 0.0625 | 1 | 0.33317013591547739 |
| no_activation | 4 | r99 | -22.5625 | -20.5625 | 2 | 0.40824829046386302 | 4.8989794855663558 | 0.00019272962573700628 |
| no_activation | 4 | radial_profile_L2 | null | null | -0.024843008446116119 | null | null | null |
| no_activation | 8 | total_mass | -16.875 | -13.292147078244305 | 3.5828529217556948 | 1.0462451349593138 | 3.4244870556984948 | 0.0037638282761674378 |
| no_activation | 8 | r_K_ratio | -0.014901491011833537 | -0.01161972280946192 | 0.0032817682023716169 | 0.0022881613015022077 | 1.4342381370653778 | 0.17202282033967264 |
| no_activation | 8 | active_fraction | -0.31580496549020443 | -0.30992378821835553 | 0.0058811772718488919 | 0.0018464748333120097 | 3.1850839046096624 | 0.0061487594061861101 |
| no_activation | 8 | r50 | -0.125 | -0.9375 | -0.8125 | 0.10077822185373186 | -8.0622577482985509 | 7.8258249630702749e-07 |
| no_activation | 8 | r90 | -4.9375 | -3.125 | 1.8125 | 0.58607984382107303 | 3.092582041694214 | 0.0074289371221557829 |
| no_activation | 8 | r99 | -43.125 | -40.875 | 2.25 | 0.60207972893961481 | 3.7370465934182984 | 0.0019828295345335756 |
| no_activation | 8 | radial_profile_L2 | null | null | -0.033609859466659862 | null | null | null |
| no_division | 0 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 0 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 0 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 0 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 0 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 0 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 0 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_division | 4 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 4 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 4 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 4 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 4 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 4 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 4 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_division | 8 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 8 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 8 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 8 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 8 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 8 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_division | 8 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_exchange | 0 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 0 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 0 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 0 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 0 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 0 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 0 | radial_profile_L2 | null | null | 0 | null | null | null |
| no_exchange | 4 | total_mass | -5 | -0.85144812456339736 | 4.1485518754366026 | 0.7592926509416037 | 5.4637060826177581 | 6.5356569411931671e-05 |
| no_exchange | 4 | r_K_ratio | -0.0039743733936599135 | -0.00044225149127364444 | 0.0035321219023862691 | 0.0011367863769178304 | 3.1071113923470066 | 0.0072116971887338042 |
| no_exchange | 4 | active_fraction | -0.004043151821485241 | -0.00065255551649256484 | 0.0033905963049926761 | 0.00099000329618615431 | 3.4248333496004126 | 0.0037611545821412203 |
| no_exchange | 4 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| no_exchange | 4 | r90 | -1.0625 | -0.125 | 0.9375 | 0.11063265039459795 | 8.4739902429904816 | 4.2084503891742595e-07 |
| no_exchange | 4 | r99 | -4.25 | -0.6875 | 3.5625 | 0.54748477908827142 | 6.50703021540187 | 9.90864543886477e-06 |
| no_exchange | 4 | radial_profile_L2 | null | null | 0.015419945058857236 | null | null | null |
| no_exchange | 8 | total_mass | -11.3125 | -1.9096311219107704 | 9.4028688780892296 | 0.92802296041635191 | 10.132151120345869 | 4.2030538877814111e-08 |
| no_exchange | 8 | r_K_ratio | -0.012214140035863512 | -0.0016402511914930182 | 0.010573888844370494 | 0.001704523993368783 | 6.2034262266220717 | 1.6894043130605425e-05 |
| no_exchange | 8 | active_fraction | -0.013089059315421731 | -0.0015910168804579676 | 0.011498042434963763 | 0.0014080908720231465 | 8.1656963079686484 | 6.6835312925587801e-07 |
| no_exchange | 8 | r50 | -0.125 | 0 | 0.125 | 0.085391256382996647 | 1.4638501094227998 | 0.16387561365565406 |
| no_exchange | 8 | r90 | -4.25 | -1 | 3.25 | 0.63574103742535504 | 5.1121444246574939 | 0.00012753562275307542 |
| no_exchange | 8 | r99 | -10.125 | -1.375 | 8.75 | 1.3616778865306827 | 6.4258956443020976 | 1.1412968631222935e-05 |
| no_exchange | 8 | radial_profile_L2 | null | null | 0.01735680763121087 | null | null | null |
| enlarged_boundary | 0 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 0 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 0 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 0 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 0 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 0 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 0 | radial_profile_L2 | null | null | 0 | null | null | null |
| enlarged_boundary | 4 | total_mass | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 4 | r_K_ratio | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 4 | active_fraction | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 4 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 4 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 4 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 4 | radial_profile_L2 | null | null | 0 | null | null | null |
| enlarged_boundary | 8 | total_mass | 0 | -9.2370555648813024e-14 | -9.2370555648813024e-14 | 9.2370555648813024e-14 | -1 | 0.33317013591547739 |
| enlarged_boundary | 8 | r_K_ratio | 0 | -1.8735013540549517e-16 | -1.8735013540549517e-16 | 1.8735013540549517e-16 | -1 | 0.33317013591547739 |
| enlarged_boundary | 8 | active_fraction | 0 | -1.3530843112619095e-16 | -1.3530843112619095e-16 | 1.3530843112619095e-16 | -1 | 0.33317013591547739 |
| enlarged_boundary | 8 | r50 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 8 | r90 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 8 | r99 | 0 | 0 | 0 | 0 | null | 1 |
| enlarged_boundary | 8 | radial_profile_L2 | null | null | 1.3877787807814457e-17 | null | null | null |

#### full: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 4 | total_mass | 961 | 957.90656585838929 | 0.003218974132789498 | 0.00052029136316337154 | 6.1868682832214148 | -2.8447421948193705 | 0.012296571838415133 | 6.5482316044096946e-33 | 2.1587448948415302e-28 | PASS |
| 4 | r_K_ratio | 0.95055031041461968 | 0.94601696199436702 | 0.0047691830412166618 | 0.0013209554508153041 | 3.6104041497183603 | -1.6136684927102185 | 0.1274345568774421 | 1.0685895014976849e-26 | 3.7935550722800251e-22 | PASS |
| 4 | active_fraction | 0.30956904900325755 | 0.30629895168952725 | 0.0032700973137303226 | 0.0044779620119824587 | 0.73026463935601882 | -2.63156267037931 | 0.018873585616829489 | 1.4630200768729918e-16 | 2.0861532987049903e-17 | PASS |
| 4 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r90 | 19.0625 | 19 | 0.0032786885245901639 | 0.016393442622950821 | 0.20000000000000001 | -1 | 0.33317013591547739 | 6.3106316067013764e-23 | 8.9016836752560409e-20 | PASS |
| 4 | r99 | 42.1875 | 40.5625 | 0.038518518518518521 | 0.014814814814814815 | 2.6000000000000001 | -4.333333333333333 | 0.00059085778625744573 | 7.2263681135141884e-16 | 1.4960608813144921e-15 | PASS |
| 4 | radial_profile_L2 | null | null | 0.04402952038479812 | 0.010491987617756247 | 4.1964899301143106 | 7.4422381775128956 | 2.0725615554004012e-06 | 1.8320674186031587e-26 | 8.1220604059720003e-25 | PASS |
| 8 | total_mass | 930.6875 | 930.17624675850664 | 0.00054932857859741886 | 0.0012087838291585521 | 0.454447325771898 | -0.25820631899454072 | 0.79975686500682674 | 4.7640347720947908e-29 | 3.9086266215654425e-24 | PASS |
| 8 | r_K_ratio | 0.89478658768318697 | 0.89019784338966224 | 0.005128311439497709 | 0.00095859770278864428 | 5.3498056844690991 | -1.1556781954739936 | 0.26588985593582204 | 3.4656152627587374e-24 | 1.9236040346766121e-19 | PASS |
| 8 | active_fraction | 0.31580496549020443 | 0.30992378821835553 | 0.0058811772718488919 | 0.0041874013475623299 | 1.4044933321886099 | -3.1850839046096624 | 0.0061487594061861101 | 1.1819392671347893e-13 | 3.6532628616503073e-15 | PASS |
| 8 | r50 | 13.125 | 13.9375 | 0.061904761904761907 | 0 | null | 8.0622577482985509 | 7.8258249630702749e-07 | 3.6873899291617687e-18 | 8.5245610309508474e-13 | PASS |
| 8 | r90 | 22.9375 | 21.125 | 0.07901907356948229 | 0.073569482288828342 | 1.0740740740740742 | -3.092582041694214 | 0.0074289371221557829 | 6.4679713206351386e-08 | 2.9962515735592815e-08 | PASS |
| 8 | r99 | 63.5 | 60.9375 | 0.040354330708661415 | 0.017716535433070866 | 2.2777777777777777 | -4.2823103365101671 | 0.00065469959890511616 | 2.3014461693173097e-16 | 1.4852885343350244e-14 | PASS |
| 8 | radial_profile_L2 | null | null | 0.055529674363131575 | 0.010917333453772486 | 5.0863770533585093 | 5.2489007583475669 | 9.814579750143665e-05 | 4.6965365295233549e-25 | 5.6802490780913782e-23 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5546875, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 1.1642856634157104e-11, "maximum_r99_to_half_width": 0.5, "outer_face_mass_present": true}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

#### no_migration: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 4 | total_mass | 952 | 949.87569077033334 | 0.0022314172580532639 | 0.002232142857142857 | 0.99967493160786225 | -1.6826567732804387 | 0.11313450170332376 | 4.4535366083474703e-32 | 2.8935620772195421e-27 | PASS |
| 4 | r_K_ratio | 0.94536916290770534 | 0.94324783734964746 | 0.0022439123691460091 | 0.00062042738925599896 | 3.6167203576180813 | -0.74650279126106678 | 0.46690363030552429 | 1.1234318912804935e-26 | 5.6912070311103157e-22 | PASS |
| 4 | active_fraction | 0.30366375606943141 | 0.30174332202349607 | 0.0019204340459353306 | 0.0039408988678125718 | 0.48730863448959377 | -1.1907017076283806 | 0.25227334259736889 | 4.5660121772795128e-15 | 1.4659870327791253e-15 | PASS |
| 4 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r90 | 17.625 | 17 | 0.035460992907801421 | 0.0035460992907801418 | 10 | -5 | 0.0001583695146220272 | 5.166732354224037e-17 | 1.4603044949211049e-15 | PASS |
| 4 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | radial_profile_L2 | null | null | 0.0037889615590438784 | 0.0021513449012044586 | 1.7612060051006133 | -3.1572000031730938 | 0.0065097809909166384 | 3.884111529407592e-27 | 5.3738655412258282e-27 | PASS |
| 8 | total_mass | 907.5625 | 910.73490678565736 | 0.0034955243144768476 | 0.00061979202534260722 | 5.6398342856131443 | 1.7126437090602438 | 0.10736597254485412 | 2.26791484418969e-29 | 2.3625497542197872e-24 | PASS |
| 8 | r_K_ratio | 0.87521057944020519 | 0.87679490526957693 | 0.0018102224385645039 | 0.0035342535956991851 | 0.51219370357785143 | 0.34631983118810106 | 0.73391239296415456 | 3.0836510768096614e-23 | 2.9189208794587623e-18 | PASS |
| 8 | active_fraction | 0.29797458110859887 | 0.29749662150093792 | 0.00047795960766096324 | 0.0042406827494155183 | 0.11270817363709631 | -0.20124803433589591 | 0.84320820817127218 | 8.6147299285259282e-13 | 6.5242620856423578e-13 | PASS |
| 8 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r90 | 17.9375 | 17.625 | 0.017421602787456445 | 0 | null | -2.6111648393354674 | 0.019657034541711152 | 3.6404802514401107e-16 | 8.8997838401990054e-17 | PASS |
| 8 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | radial_profile_L2 | null | null | 0.0081916778710457352 | 0.0053869623501989275 | 1.5206488069743491 | -5.2320178687917576 | 0.00010135847594577555 | 3.0785129307496989e-25 | 6.2103157867586928e-25 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1484375, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1484375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1484375, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

#### no_activation: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 4 | total_mass | 952.6875 | 951.86845444722951 | 0.00085972110767744338 | 0.001443285442498196 | 0.59566949292396676 | -0.61012194174755463 | 0.55091319382061887 | 1.0530942844662949e-31 | 7.5814778380471533e-27 | PASS |
| 4 | r_K_ratio | 0.94577685712252779 | 0.94390323540480114 | 0.0019810399288337499 | 0.0013675401410755857 | 1.4486155611312732 | -0.65413362820031062 | 0.52292895052316557 | 1.2864824010724779e-26 | 6.33873863312491e-22 | PASS |
| 4 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r90 | 18 | 18 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r99 | 19.625 | 20 | 0.019108280254777069 | 0.006369426751592357 | 3 | 3 | 0.0089727374772233335 | 9.9141794617064285e-21 | 6.7481503001323698e-16 | PASS |
| 4 | radial_profile_L2 | null | null | 0.019186511938682 | 0.0032809938768473233 | 5.8477743814377803 | 0.5702783058743085 | 0.57693204092142802 | 1.0558375790990729e-27 | 5.4736628716398694e-27 | PASS |
| 8 | total_mass | 913.8125 | 916.88409968026235 | 0.0033613018866149344 | 0.00068394774639217568 | 4.9145594884196955 | 1.6521457270470532 | 0.1192781088232462 | 1.7748297492024395e-29 | 2.4726533817824996e-24 | PASS |
| 8 | r_K_ratio | 0.87988509667135351 | 0.87857812058020024 | 0.0014853940544027194 | 0.00078428743000585487 | 1.8939409170329693 | -0.2666446731546403 | 0.79337351160724512 | 9.5847369897507056e-23 | 6.4341241194218808e-18 | PASS |
| 8 | active_fraction | 0 | 0 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r90 | 18 | 18 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 8 | r99 | 20.375 | 20.0625 | 0.015337423312883436 | 0.0061349693251533744 | 2.5 | -2.6111648393354674 | 0.019657034541711152 | 4.0101538088962186e-20 | 2.9572284872306864e-17 | PASS |
| 8 | radial_profile_L2 | null | null | 0.021919814896471712 | 0.01063483622589392 | 2.061133282250355 | -0.91363059411633851 | 0.37536315590431879 | 1.1654071036418924e-23 | 7.6292068230317077e-23 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1640625, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1640625, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.1640625, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

#### no_division: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 4 | total_mass | 961 | 957.90656585838929 | 0.003218974132789498 | 0.00052029136316337154 | 6.1868682832214148 | -2.8447421948193705 | 0.012296571838415133 | 6.5482316044096946e-33 | 2.1587448948415302e-28 | PASS |
| 4 | r_K_ratio | 0.95055031041461968 | 0.94601696199436702 | 0.0047691830412166618 | 0.0013209554508153041 | 3.6104041497183603 | -1.6136684927102185 | 0.1274345568774421 | 1.0685895014976849e-26 | 3.7935550722800251e-22 | PASS |
| 4 | active_fraction | 0.30956904900325755 | 0.30629895168952725 | 0.0032700973137303226 | 0.0044779620119824587 | 0.73026463935601882 | -2.63156267037931 | 0.018873585616829489 | 1.4630200768729918e-16 | 2.0861532987049903e-17 | PASS |
| 4 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r90 | 19.0625 | 19 | 0.0032786885245901639 | 0.016393442622950821 | 0.20000000000000001 | -1 | 0.33317013591547739 | 6.3106316067013764e-23 | 8.9016836752560409e-20 | PASS |
| 4 | r99 | 42.1875 | 40.5625 | 0.038518518518518521 | 0.014814814814814815 | 2.6000000000000001 | -4.333333333333333 | 0.00059085778625744573 | 7.2263681135141884e-16 | 1.4960608813144921e-15 | PASS |
| 4 | radial_profile_L2 | null | null | 0.04402952038479812 | 0.010491987617756247 | 4.1964899301143106 | 7.4422381775128956 | 2.0725615554004012e-06 | 1.8320674186031587e-26 | 8.1220604059720003e-25 | PASS |
| 8 | total_mass | 930.6875 | 930.17624675850664 | 0.00054932857859741886 | 0.0012087838291585521 | 0.454447325771898 | -0.25820631899454072 | 0.79975686500682674 | 4.7640347720947908e-29 | 3.9086266215654425e-24 | PASS |
| 8 | r_K_ratio | 0.89478658768318697 | 0.89019784338966224 | 0.005128311439497709 | 0.00095859770278864428 | 5.3498056844690991 | -1.1556781954739936 | 0.26588985593582204 | 3.4656152627587374e-24 | 1.9236040346766121e-19 | PASS |
| 8 | active_fraction | 0.31580496549020443 | 0.30992378821835553 | 0.0058811772718488919 | 0.0041874013475623299 | 1.4044933321886099 | -3.1850839046096624 | 0.0061487594061861101 | 1.1819392671347893e-13 | 3.6532628616503073e-15 | PASS |
| 8 | r50 | 13.125 | 13.9375 | 0.061904761904761907 | 0 | null | 8.0622577482985509 | 7.8258249630702749e-07 | 3.6873899291617687e-18 | 8.5245610309508474e-13 | PASS |
| 8 | r90 | 22.9375 | 21.125 | 0.07901907356948229 | 0.073569482288828342 | 1.0740740740740742 | -3.092582041694214 | 0.0074289371221557829 | 6.4679713206351386e-08 | 2.9962515735592815e-08 | PASS |
| 8 | r99 | 63.5 | 60.9375 | 0.040354330708661415 | 0.017716535433070866 | 2.2777777777777777 | -4.2823103365101671 | 0.00065469959890511616 | 2.3014461693173097e-16 | 1.4852885343350244e-14 | PASS |
| 8 | radial_profile_L2 | null | null | 0.055529674363131575 | 0.010917333453772486 | 5.0863770533585093 | 5.2489007583475669 | 9.814579750143665e-05 | 4.6965365295233549e-25 | 5.6802490780913782e-23 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5546875, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 1.1642856634157104e-11, "maximum_r99_to_half_width": 0.5, "outer_face_mass_present": true}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

#### no_exchange: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 4 | total_mass | 956 | 957.05511773382591 | 0.0011036796378931959 | 0.0013729079497907949 | 0.8038992257721107 | 0.99005614211195037 | 0.33784950239991707 | 4.8494654479605545e-33 | 1.9533763378451138e-28 | PASS |
| 4 | r_K_ratio | 0.94657593702095977 | 0.94557471050309339 | 0.0010577350202007078 | 0.0020059806700443946 | 0.52729073415114158 | -0.37591697201234414 | 0.71224193676460312 | 4.4157912633677212e-27 | 2.1141554501660207e-22 | PASS |
| 4 | active_fraction | 0.30552589718177231 | 0.30564639617303468 | 0.00012049899126235358 | 0.0037872688022809586 | 0.031816857358997239 | 0.087800156593427464 | 0.93119695228041333 | 2.260056809173796e-16 | 2.4276628496373229e-16 | PASS |
| 4 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r90 | 18 | 18.875 | 0.048611111111111112 | 0 | null | 10.2469507659596 | 3.6215788445874798e-08 | 6.7648146393204853e-20 | 2.4126007329805088e-17 | PASS |
| 4 | r99 | 37.9375 | 39.875 | 0.051070840197693576 | 0.018121911037891267 | 2.8181818181818183 | 4.3811422547920502 | 0.00053681085850166654 | 4.4944764487226121e-17 | 3.0928846255920417e-11 | PASS |
| 4 | radial_profile_L2 | null | null | 0.059449465443655355 | 0.0033724803913630971 | 17.627816486614748 | 12.833927211716004 | 1.7166672211365084e-09 | 3.2722074614377446e-29 | 5.6087664478374278e-27 | PASS |
| 8 | total_mass | 919.375 | 928.26661563659582 | 0.0096713698290641411 | 0.00013596193065941537 | 71.132925092766754 | 5.2908451995786665 | 9.0612275000866066e-05 | 3.4065489061979735e-30 | 6.1263493639015173e-25 | PASS |
| 8 | r_K_ratio | 0.88257244764732345 | 0.88855759219816921 | 0.0067814767691879657 | 0.00082008136138378006 | 8.2692731337596932 | 1.3785299433809588 | 0.18825359641013961 | 1.0279495071624366e-23 | 1.4525171616692975e-18 | PASS |
| 8 | active_fraction | 0.30271590617478272 | 0.3083327713378976 | 0.005616865163114871 | 0.0036498960976138368 | 1.5389109752430936 | 2.886996510534817 | 0.011288371895751954 | 8.4755872608630673e-15 | 2.3261723245254334e-13 | PASS |
| 8 | r50 | 13 | 13.9375 | 0.072115384615384609 | 0 | null | 15 | 1.9412759304364392e-10 | 2.6605239197581584e-20 | 1.8615816799604571e-16 | PASS |
| 8 | r90 | 18.6875 | 20.125 | 0.076923076923076927 | 0.0033444816053511705 | 23 | 11.22285083890813 | 1.0733435740340221e-08 | 1.007788204400195e-18 | 3.287104854231798e-13 | PASS |
| 8 | r99 | 53.375 | 59.5625 | 0.11592505854800937 | 0.01873536299765808 | 6.1875 | 6.3984255582182206 | 1.1974874035688459e-05 | 8.7057171773066584e-16 | 8.2754179948232504e-07 | PASS |
| 8 | radial_profile_L2 | null | null | 0.072886481994342445 | 0.0087813469351496597 | 8.3001483180894553 | 10.528212630791504 | 2.5280330282908069e-08 | 3.9750330424795307e-26 | 2.2416312130719471e-23 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.4765625, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.4609375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 1.1582572135901263e-11, "maximum_r99_to_half_width": 0.4921875, "outer_face_mass_present": true}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

#### enlarged_boundary: all equivalence contrasts

16 paired realizations plus 16 disjoint ABM baseline seeds. Result: PASS.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 4 | total_mass | 961 | 957.90656585838929 | 0.003218974132789498 | 0.00052029136316337154 | 6.1868682832214148 | -2.8447421948193705 | 0.012296571838415133 | 6.5482316044096946e-33 | 2.1587448948415302e-28 | PASS |
| 4 | r_K_ratio | 0.95055031041461968 | 0.94601696199436702 | 0.0047691830412166618 | 0.0013209554508153041 | 3.6104041497183603 | -1.6136684927102185 | 0.1274345568774421 | 1.0685895014976849e-26 | 3.7935550722800251e-22 | PASS |
| 4 | active_fraction | 0.30956904900325755 | 0.30629895168952725 | 0.0032700973137303226 | 0.0044779620119824587 | 0.73026463935601882 | -2.63156267037931 | 0.018873585616829489 | 1.4630200768729918e-16 | 2.0861532987049903e-17 | PASS |
| 4 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r90 | 19.0625 | 19 | 0.0032786885245901639 | 0.016393442622950821 | 0.20000000000000001 | -1 | 0.33317013591547739 | 6.3106316067013764e-23 | 8.9016836752560409e-20 | PASS |
| 4 | r99 | 42.1875 | 40.5625 | 0.038518518518518521 | 0.014814814814814815 | 2.6000000000000001 | -4.333333333333333 | 0.00059085778625744573 | 7.2263681135141884e-16 | 1.4960608813144921e-15 | PASS |
| 4 | radial_profile_L2 | null | null | 0.04402952038479812 | 0.010491987617756247 | 4.1964899301143106 | 7.4422381775128956 | 2.0725615554004012e-06 | 1.8320674186031587e-26 | 8.1220604059720003e-25 | PASS |
| 8 | total_mass | 930.6875 | 930.17624675850652 | 0.00054932857859751807 | 0.0012087838291585521 | 0.4544473257719801 | -0.25820631899458762 | 0.79975686500678955 | 4.7640347720946557e-29 | 3.9086266215654139e-24 | PASS |
| 8 | r_K_ratio | 0.89478658768318697 | 0.89019784338966201 | 0.0051283114394979181 | 0.00095859770278864428 | 5.3498056844693176 | -1.1556781954740531 | 0.26588985593579861 | 3.4656152627579491e-24 | 1.9236040346763526e-19 | PASS |
| 8 | active_fraction | 0.31580496549020443 | 0.30992378821835542 | 0.0058811772718490272 | 0.0041874013475623299 | 1.4044933321886421 | -3.1850839046097859 | 0.0061487594061845367 | 1.1819392671345702e-13 | 3.6532628616493315e-15 | PASS |
| 8 | r50 | 13.125 | 13.9375 | 0.061904761904761907 | 0 | null | 8.0622577482985509 | 7.8258249630702749e-07 | 3.6873899291617687e-18 | 8.5245610309508474e-13 | PASS |
| 8 | r90 | 22.9375 | 21.125 | 0.07901907356948229 | 0.073569482288828342 | 1.0740740740740742 | -3.092582041694214 | 0.0074289371221557829 | 6.4679713206351386e-08 | 2.9962515735592815e-08 | PASS |
| 8 | r99 | 63.5 | 60.9375 | 0.040354330708661415 | 0.017716535433070866 | 2.2777777777777777 | -4.2823103365101671 | 0.00065469959890511616 | 2.3014461693173097e-16 | 1.4852885343350244e-14 | PASS |
| 8 | radial_profile_L2 | null | null | 0.055529674363131588 | 0.010917333453772486 | 5.0863770533585102 | 5.2489007583475633 | 9.8145797501437233e-05 | 4.6965365295239573e-25 | 5.6802490780921858e-23 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.2734375, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.27734375, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.25, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

