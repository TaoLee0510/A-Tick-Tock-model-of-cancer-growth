# Phase D report

Hybrid schema 4, `hybrid_invasion_front_v4`, keeps every current r cell
as an individual. Its K core requires occupied-volume hysteresis, distance
from the moving front, a protected halo around active r, and a small nutrient
gradient. Large K cells require a complete core footprint. Regional
conversion uses cached membership/capacity and deterministic water filling;
3D integral sums replace per-agent density scans. The PDE owns the single
vascular process when the new shared vascular law is selected.

The external capacity check now respects the native grid's accepted swap
and retained self-overlap. Ordinary PDE transport reflects at the core
interface, keeping new fractional tails from diffusing into the agent front.
Both changes apply only to the new hybrid policy. PDE-only and native ABM
limits retain their own operators and checkpoints.

## Tests and exact arithmetic

| Build | Tests passed | Validation suites | Full run + retry seconds |
| --- | ---: | ---: | ---: |
| release | 60/60 | 12 | 5001.4000000000005 |
| hdf5 | 63/63 | 12 | 4974.1499999999996 |

New assert tests cover 2D/3D classification, mass-conserving whole/fractional
K conversion, unchanged front density under reflecting core transport,
native occupied-site swaps, exact density queries, active-r representation
and one-to-four-thread continuation including dynamic vasculature. Python
tests reject missing traces and failed coverage; native all-ABM/all-PDE
sampling and adaptive restart have exact endpoint/trace checks. C++ tests
use `-UNDEBUG`; checked build logs have no new warnings.

All 416 published native reports, 96 schema-14 traces and 48 schema-15
vascular reports remain exact in both builds. All 176 new invasion reports
are exactly equal between Release and HDF5 builds. Full numerical data,
including every native realization and failed pilot, are saved in
[the numeric report](phase_d_numeric_results.json). Existing phase B/C
tables retain their values in their respective numeric reports.

The initial complete run passed 59/60 Release and 62/63 HDF5 checks.
The compatibility configuration scanner then encountered the later phase
E/F prefix-guidance YAML in the shared working directory. Recompiling only
that unchanged test against immutable phase-D headers/configurations and
rerunning the failed check passed in both builds. Simulator libraries and
all completed ensembles were unchanged. The initial failure and retry are
retained explicitly in the numeric result; the table reports their combined
test-set outcome and elapsed time, rather than a fresh all-green full run.

## Prespecified invasion and refinement

The case is 256-square, 8-hour r20 early invasion, with independent ABM
groups of 16 seeds each and paired hybrid/PDE realizations. Samples are
at 0, 4 and 8 hours. Native PDE and fine-hybrid references also have
independent baseline groups. All original phase B scalar/profile margins
and the earlier refinement margins are unchanged.

# Invasion-front hybrid equivalence

Result: PASS.
This 8-hour case measures early r20 invasion, not long-time production accuracy.

## primary_hybrid_vs_abm

Column ABM denotes the declared reference; PDE denotes the candidate in this table.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 4 | total_mass | 961 | 960.60420146538536 | 0.00041186111822539645 | 0.00052029136316337154 | 0.79159706922921202 | -0.59911491992162458 | 0.55803709401791424 | 7.2339225692645515e-33 | 2.501766620876796e-32 | PASS |
| 4 | r_K_ratio | 0.95055031041461968 | 0.95032700607307152 | 0.00023492111790559915 | 0.0013209554508153041 | 0.17784181726991927 | -0.21987739066059472 | 0.8289309634296983 | 4.1970173122572531e-29 | 5.8140243264515875e-28 | PASS |
| 4 | active_fraction | 0.30956904900325755 | 0.30973743638702933 | 0.00016838738377180501 | 0.0044779620119824587 | 0.0376035757608532 | 0.3789487310739223 | 0.7100361612944639 | 1.0769482826894484e-23 | 1.1913077660463727e-23 | PASS |
| 4 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r90 | 19.0625 | 19.125 | 0.0032786885245901639 | 0.016393442622950821 | 0.20000000000000001 | 1 | 0.33317013591547739 | 3.9224844895571197e-21 | 9.2878684736242714e-21 | PASS |
| 4 | r99 | 42.1875 | 42.0625 | 0.0029629629629629628 | 0.014814814814814815 | 0.20000000000000001 | -1.4638501094227998 | 0.16387561365565406 | 7.7121077027247792e-21 | 2.5289797274610963e-25 | PASS |
| 4 | radial_profile_L2 | null | null | 0.013563834218350852 | 0.010491987617756245 | 1.2927802350239082 | -0.92436113362330252 | 0.36993041893514766 | 9.0822379811660861e-24 | 2.9024984987133047e-23 | PASS |
| 8 | total_mass | 930.6875 | 928.9731923002405 | 0.0018419799339300113 | 0.0012087838291585521 | 1.5238290664528777 | -1.9596651241171705 | 0.068886579380291885 | 3.8592581218310392e-31 | 8.2958206642171979e-29 | PASS |
| 8 | r_K_ratio | 0.89478658768318697 | 0.89190334341500987 | 0.003222270324416123 | 0.00095859770278864428 | 3.3614417341521454 | -1.589585985829981 | 0.13277931128890752 | 2.7208843840172391e-25 | 1.7675424616172034e-24 | PASS |
| 8 | active_fraction | 0.31580496549020443 | 0.31507229203460591 | 0.00073267345559853425 | 0.0041874013475623299 | 0.17497091747010485 | -0.77450021299384708 | 0.45067148277509717 | 1.1464559179817923e-18 | 7.4026323066234482e-19 | PASS |
| 8 | r50 | 13.125 | 13.8125 | 0.052380952380952382 | 0 | null | 5.7445626465380286 | 3.8761814875216441e-05 | 2.4119653173897488e-16 | 2.2306530993797326e-12 | PASS |
| 8 | r90 | 22.9375 | 20.625 | 0.1008174386920981 | 0.073569482288828342 | 1.3703703703703705 | -3.7464952863928409 | 0.0019448517147768215 | 4.9819053680939549e-07 | 2.5921042179874925e-08 | PASS |
| 8 | r99 | 63.5 | 62.8125 | 0.010826771653543307 | 0.017716535433070866 | 0.61111111111111116 | -1.4569855927715483 | 0.16573520407028275 | 4.0947084130739473e-16 | 1.4999871033774585e-16 | PASS |
| 8 | radial_profile_L2 | null | null | 0.040654172338689773 | 0.010917333453772484 | 3.7238188712319422 | 3.0941604499388822 | 0.0074050279475077335 | 1.0232500146242289e-23 | 3.3734713086728949e-22 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5546875, "outer_face_mass_present": false}, "coarse_exchange": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "coarse_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "fine_exchange": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.53125, "outer_face_mass_present": false}, "fine_exchange_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5234375, "outer_face_mass_present": false}, "fine_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "fine_step_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5625, "outer_face_mass_present": false}, "middle_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 1.1642856634157104e-11, "maximum_r99_to_half_width": 0.5, "outer_face_mass_present": true}, "pde_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 8.178909993925213e-10, "maximum_r99_to_half_width": 0.4921875, "outer_face_mass_present": true}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.
## secondary_hybrid_vs_pde

Column ABM denotes the declared reference; PDE denotes the candidate in this table.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 4 | total_mass | 957.90656585838929 | 960.60420146538536 | 0.0028161782194056936 | 8.5074110443608325e-05 | 33.102646677362635 | 2.1691663801791687 | 0.046550379617914646 | 1.7062466868890748e-29 | 3.231673025339014e-29 | PASS |
| 4 | r_K_ratio | 0.94601696199436702 | 0.95032700607307152 | 0.0045559902748659011 | 0.00011353351945977469 | 40.129032346963442 | 1.5257950649582155 | 0.14786747965403169 | 4.8784537824099061e-24 | 7.8346658415227234e-24 | PASS |
| 4 | active_fraction | 0.30629895168952725 | 0.30973743638702933 | 0.0034384846975021276 | 0.0053910646293073919 | 0.63781181156863287 | 2.6656188099945615 | 0.017632075789053698 | 3.4702551467809271e-17 | 2.6875234802061689e-16 | PASS |
| 4 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 4 | r90 | 19 | 19.125 | 0.0065789473684210523 | 0 | null | 1.4638501094227998 | 0.16387561365565406 | 2.909519278220455e-19 | 6.3856561125241839e-19 | PASS |
| 4 | r99 | 40.5625 | 42.0625 | 0.036979969183359017 | 0.0015408320493066256 | 24 | 3.7696851746252595 | 0.0018547237781325861 | 1.9742433772332595e-16 | 1.2120756276592968e-13 | PASS |
| 4 | radial_profile_L2 | null | null | 0.040356939806304455 | 0.0016640760607995906 | 24.251860090405298 | 9.9799541723681262 | 5.1305428573420625e-08 | 4.0764695175364698e-26 | 1.3130557083714211e-24 | PASS |
| 8 | total_mass | 930.17624675850664 | 928.9731923002405 | 0.001293361835951548 | 0.00022003771420046962 | 5.8779097967415064 | -0.63928892150589156 | 0.53227717939788732 | 2.1418828179601383e-26 | 1.579148470839528e-26 | PASS |
| 8 | r_K_ratio | 0.89019784338966224 | 0.89190334341500987 | 0.0019158662740111447 | 8.4126535403164092e-05 | 22.773626238494625 | 0.40661018691245382 | 0.69003667472061503 | 5.6530596755188373e-21 | 5.8452837168232881e-21 | PASS |
| 8 | active_fraction | 0.30992378821835553 | 0.31507229203460591 | 0.0051485038162503577 | 0.0053629428447906913 | 0.96001467202869228 | 2.5683633464687325 | 0.021403776135631054 | 1.4911216461572069e-14 | 3.0868564101523652e-13 | PASS |
| 8 | r50 | 13.9375 | 13.8125 | 0.0089686098654708519 | 0.0044843049327354259 | 2 | -1.4638501094227998 | 0.16387561365565406 | 7.3596621204324012e-17 | 4.2322516392401768e-17 | PASS |
| 8 | r90 | 21.125 | 20.625 | 0.023668639053254437 | 0.0059171597633136093 | 4 | -1.3693063937629153 | 0.19105762734590448 | 4.1250338146898709e-11 | 6.3343600943196488e-10 | PASS |
| 8 | r99 | 60.9375 | 62.8125 | 0.030769230769230771 | 0.0061538461538461538 | 5 | 2.7235238970096107 | 0.015699869495500712 | 1.1633714381568576e-14 | 1.2409396409401457e-13 | PASS |
| 8 | radial_profile_L2 | null | null | 0.047467027917822362 | 0.0027442327123833759 | 17.29701264168558 | 9.8907136439287484 | 5.773053540374371e-08 | 1.8383276460934294e-22 | 1.0927784274552424e-20 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5546875, "outer_face_mass_present": false}, "coarse_exchange": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "coarse_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "fine_exchange": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.53125, "outer_face_mass_present": false}, "fine_exchange_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5234375, "outer_face_mass_present": false}, "fine_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "fine_step_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5625, "outer_face_mass_present": false}, "middle_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 1.1642856634157104e-11, "maximum_r99_to_half_width": 0.5, "outer_face_mass_present": true}, "pde_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 8.178909993925213e-10, "maximum_r99_to_half_width": 0.4921875, "outer_face_mass_present": true}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.
## native_pde_vs_abm

Column ABM denotes the declared reference; PDE denotes the candidate in this table.

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

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5546875, "outer_face_mass_present": false}, "coarse_exchange": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "coarse_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "fine_exchange": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.53125, "outer_face_mass_present": false}, "fine_exchange_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5234375, "outer_face_mass_present": false}, "fine_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "fine_step_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5625, "outer_face_mass_present": false}, "middle_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 1.1642856634157104e-11, "maximum_r99_to_half_width": 0.5, "outer_face_mass_present": true}, "pde_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 8.178909993925213e-10, "maximum_r99_to_half_width": 0.4921875, "outer_face_mass_present": true}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.
## refinement_coarse_step

Column ABM denotes the declared reference; PDE denotes the candidate in this table.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 4 | total_mass | 961.09859140857247 | 960.60420146538536 | 0.00051440086127117995 | 0.00082920809651588221 | 0.62035195197991821 | -2.3517961579269691 | 0.032763565449434186 | 1.2380203479524052e-31 | 3.0530606520853109e-32 | PASS |
| 4 | r_K_ratio | 0.95036786992961819 | 0.95032700607307152 | 4.2997935683258415e-05 | 0.0015074897199713149 | 0.028522871574922976 | -0.1032890532784852 | 0.91910185787706122 | 2.7917366420345098e-28 | 2.0995014673460491e-26 | PASS |
| 4 | active_fraction | 0.30996371999862754 | 0.30973743638702933 | 0.00022628361159815827 | 0.0043339351696545028 | 0.052212043498610382 | -1.084700133310601 | 0.29518112404129682 | 1.4792058379236836e-22 | 1.0539932960408174e-22 | PASS |
| 4 | r50 | 13.0625 | 13 | 0.0047846889952153108 | 0.0047846889952153108 | 1 | -1 | 0.33317013591547739 | 3.6686536711016011e-13 | 1.6852334329393012e-12 | PASS |
| 4 | r90 | 19.0625 | 19.125 | 0.0032786885245901639 | 0.0032786885245901639 | 1 | 1 | 0.33317013591547739 | 1.972130890296023e-15 | 6.3169801777948548e-15 | PASS |
| 4 | r99 | 42.0625 | 42.0625 | 0 | 0.019316493313521546 | 0 | null | 1 | 5.0441174529444768e-23 | 5.0441174529444768e-23 | PASS |
| 4 | radial_profile_L2 | null | null | 0.0041144225883586703 | 0.0051035490015556477 | 0.80618851452283991 | -5.5998753250842386 | 5.0665710554115387e-05 | 3.0427537139168067e-19 | 6.9038649743204517e-19 | PASS |
| 8 | total_mass | 929.91646820733433 | 928.9731923002405 | 0.0010143662784166028 | 0.0013023118967302508 | 0.77889657689789926 | -1.6634451096586069 | 0.11697006322578575 | 1.0397231967515626e-25 | 7.0417433387180797e-26 | PASS |
| 8 | r_K_ratio | 0.89427426465975901 | 0.89190334341500987 | 0.0026512238341681428 | 0.003993743681274739 | 0.66384426386670714 | -1.6543054832466657 | 0.11883393390732373 | 2.5167950164340119e-19 | 7.0059902600221096e-20 | PASS |
| 8 | active_fraction | 0.31557838097293833 | 0.31507229203460591 | 0.000506088938332417 | 0.0027147641258388512 | 0.18642096140711215 | -0.97155415899204056 | 0.34667953298624743 | 1.5724172050735461e-16 | 7.4120722490762035e-17 | PASS |
| 8 | r50 | 13.9375 | 13.8125 | 0.0089686098654708519 | 0.0089686098654708519 | 1 | -1 | 0.33317013591547739 | 1.3869436397432038e-08 | 2.7589657260040533e-09 | PASS |
| 8 | r90 | 20.5 | 20.625 | 0.0060975609756097563 | 0.024390243902439025 | 0.25 | 0.69560834364025248 | 0.49731054576338629 | 1.2627692957964178e-09 | 1.8392660835542127e-08 | PASS |
| 8 | r99 | 63.125 | 62.8125 | 0.0049504950495049506 | 0.01782178217821782 | 0.27777777777777779 | -1.3206763594884356 | 0.20640537248738033 | 1.9797649689027467e-13 | 1.5491764115596632e-14 | PASS |
| 8 | radial_profile_L2 | null | null | 0.007121085843927677 | 0.0080656150377357624 | 0.8828943373333572 | -1.2423864577438846 | 0.23316710241424599 | 1.5844895188749403e-13 | 6.3364858445956056e-13 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5546875, "outer_face_mass_present": false}, "coarse_exchange": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "coarse_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "fine_exchange": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.53125, "outer_face_mass_present": false}, "fine_exchange_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5234375, "outer_face_mass_present": false}, "fine_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "fine_step_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5625, "outer_face_mass_present": false}, "middle_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 1.1642856634157104e-11, "maximum_r99_to_half_width": 0.5, "outer_face_mass_present": true}, "pde_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 8.178909993925213e-10, "maximum_r99_to_half_width": 0.4921875, "outer_face_mass_present": true}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.
## refinement_middle_step

Column ABM denotes the declared reference; PDE denotes the candidate in this table.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 4 | total_mass | 961.09859140857247 | 960.77865103107206 | 0.00033289027822996052 | 0.00082920809651588221 | 0.40145565344655859 | -1.6942330915665964 | 0.11087641155145156 | 7.6972018721969417e-32 | 3.5086034833865922e-33 | PASS |
| 4 | r_K_ratio | 0.95036786992961819 | 0.95048738605755023 | 0.00012575775309080895 | 0.0015074897199713149 | 0.083421963960856688 | 0.40459954216030308 | 0.69148259463187078 | 1.0549094623378141e-28 | 5.7992829703384536e-28 | PASS |
| 4 | active_fraction | 0.30996371999862754 | 0.30991544394758685 | 4.8276051040654216e-05 | 0.0043339351696545028 | 0.011139080108690859 | -0.47292951899019847 | 0.64307175817471962 | 2.880287433584437e-27 | 2.679157814425639e-27 | PASS |
| 4 | r50 | 13.0625 | 13 | 0.0047846889952153108 | 0.0047846889952153108 | 1 | -1 | 0.33317013591547739 | 3.6686536711016011e-13 | 1.6852334329393012e-12 | PASS |
| 4 | r90 | 19.0625 | 19.0625 | 0 | 0.0032786885245901639 | 0 | null | 1 | 3.6418671611546434e-30 | 3.6418671611546434e-30 | PASS |
| 4 | r99 | 42.0625 | 42.0625 | 0 | 0.019316493313521546 | 0 | null | 1 | 5.0441174529444768e-23 | 5.0441174529444768e-23 | PASS |
| 4 | radial_profile_L2 | null | null | 0.0039324813502922724 | 0.0051035490015556477 | 0.77053857013885552 | -8.9265243996231138 | 2.1780467269242795e-07 | 4.073120873488489e-23 | 8.935128466795645e-23 | PASS |
| 8 | total_mass | 929.91646820733433 | 929.13141548265355 | 0.00084421854168698833 | 0.0013023118967302508 | 0.6482460490506079 | -1.1235356859770231 | 0.27887062426995179 | 2.9655897309944008e-24 | 7.5004457585662478e-25 | PASS |
| 8 | r_K_ratio | 0.89427426465975901 | 0.89279820228143036 | 0.0016505701177594555 | 0.003993743681274739 | 0.41328894628325763 | -1.442822739730615 | 0.16962719443962171 | 9.6233056379915092e-22 | 1.9379908269140985e-21 | PASS |
| 8 | active_fraction | 0.31557838097293833 | 0.31541176626189538 | 0.00016661471104297693 | 0.0027147641258388512 | 0.061373549715481687 | -0.22454664515865377 | 0.82536205068415169 | 2.2791101422030838e-14 | 1.783569517558701e-14 | PASS |
| 8 | r50 | 13.9375 | 13.75 | 0.013452914798206279 | 0.0089686098654708519 | 1.5 | -1.8605210188381269 | 0.082530706346961857 | 2.0675638264748815e-09 | 5.8764970321956245e-11 | PASS |
| 8 | r90 | 20.5 | 20.4375 | 0.0030487804878048782 | 0.024390243902439025 | 0.125 | -0.43574467033059511 | 0.66922749729499076 | 2.3479500096998271e-10 | 2.2881554340863853e-10 | PASS |
| 8 | r99 | 63.125 | 62.6875 | 0.0069306930693069308 | 0.01782178217821782 | 0.3888888888888889 | -1.8154790347345071 | 0.089492103237647933 | 5.2327065955271754e-13 | 8.7531197700087064e-15 | PASS |
| 8 | radial_profile_L2 | null | null | 0.0088158856124367722 | 0.0080656150377357624 | 1.0930208758031217 | -2.0250927239903658 | 0.061035699178193167 | 1.8424283778447887e-16 | 1.0549196350743843e-15 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5546875, "outer_face_mass_present": false}, "coarse_exchange": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "coarse_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "fine_exchange": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.53125, "outer_face_mass_present": false}, "fine_exchange_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5234375, "outer_face_mass_present": false}, "fine_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "fine_step_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5625, "outer_face_mass_present": false}, "middle_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 1.1642856634157104e-11, "maximum_r99_to_half_width": 0.5, "outer_face_mass_present": true}, "pde_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 8.178909993925213e-10, "maximum_r99_to_half_width": 0.4921875, "outer_face_mass_present": true}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.
## refinement_coarse_exchange

Column ABM denotes the declared reference; PDE denotes the candidate in this table.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 4 | total_mass | 961.02174565544249 | 960.2657682235357 | 0.0007866392569412556 | 0.00048335018898518136 | 1.6274727410220844 | -1.5748302655478348 | 0.13614774055261647 | 3.7340213550499834e-27 | 2.7043313587308583e-27 | PASS |
| 4 | r_K_ratio | 0.94655770730419042 | 0.95096516786738772 | 0.0046563041314721448 | 0.00094404091210744344 | 4.9323118010612239 | 3.0770317107764233 | 0.0076685953338274685 | 3.785031561116966e-21 | 3.7538094718903231e-19 | PASS |
| 4 | active_fraction | 0.30851441942057212 | 0.30945844804771572 | 0.00094402862714358313 | 0.0056739961346165477 | 0.16637808781436211 | 1.8640281877844618 | 0.082009484156249035 | 3.5614906975091102e-17 | 1.4505412793611066e-16 | PASS |
| 4 | r50 | 13.1875 | 13 | 0.014218009478672985 | 0 | null | -1.8605210188381269 | 0.082530706346961857 | 1.2728796844468052e-09 | 3.8879211368783638e-10 | PASS |
| 4 | r90 | 19.0625 | 19.125 | 0.0032786885245901639 | 0.013114754098360656 | 0.25 | 0.43574467033059511 | 0.66922749729499076 | 1.3987630000954023e-10 | 1.9459692649303824e-09 | PASS |
| 4 | r99 | 42.25 | 42.125 | 0.0029585798816568047 | 0.0073964497041420114 | 0.40000000000000002 | -0.80757285308724824 | 0.43195708095606955 | 2.2008749221486003e-14 | 3.0726601236896891e-14 | PASS |
| 4 | radial_profile_L2 | null | null | 0.038110544529397494 | 0.014407910313865093 | 2.6451125596418215 | 1.612763410143164 | 0.12763204336850598 | 1.1342715112485703e-16 | 2.4142411119349826e-13 | PASS |
| 8 | total_mass | 928.51397835569628 | 928.20804694101662 | 0.00032948498548343249 | 0.00177604905004365 | 0.18551570153726032 | -0.27203829781178812 | 0.78930135228982234 | 4.448640789676522e-22 | 4.0934656615122109e-21 | PASS |
| 8 | r_K_ratio | 0.88434265759862962 | 0.89472064670506368 | 0.011735257840682106 | 0.0014968649513863032 | 7.8398908530884102 | 4.5213375934936213 | 0.00040568278294048478 | 2.6602963558468482e-18 | 5.0876906933659199e-15 | PASS |
| 8 | active_fraction | 0.31191934111235159 | 0.31438920500194562 | 0.0024698638895940435 | 0.0061020419805564753 | 0.40476022575131554 | 2.4513800557233356 | 0.026968155040624997 | 3.2370393453070971e-13 | 1.175416655890267e-11 | PASS |
| 8 | r50 | 14 | 13.3125 | 0.049107142857142856 | 0 | null | -5.7445626465380286 | 3.8761814875216441e-05 | 1.323207035466362e-05 | 1.1339602890069579e-11 | PASS |
| 8 | r90 | 20.3125 | 21 | 0.033846153846153845 | 0.046153846153846156 | 0.73333333333333328 | 1.7891501337389775 | 0.093798437067253118 | 1.9408301826110954e-06 | 0.0016122726158032661 | PASS |
| 8 | r99 | 62.4375 | 62.9375 | 0.0080080080080080079 | 0.023023023023023025 | 0.34782608695652173 | 1.0540925533894598 | 0.30852450532643372 | 2.990794702427626e-10 | 1.801335712350327e-09 | PASS |
| 8 | radial_profile_L2 | null | null | 0.059722567950992386 | 0.012294047784762041 | 4.8578441369827781 | 3.2550596469914201 | 0.0053277805034846531 | 3.2578014400536092e-19 | 8.746919235877852e-14 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5546875, "outer_face_mass_present": false}, "coarse_exchange": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "coarse_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "fine_exchange": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.53125, "outer_face_mass_present": false}, "fine_exchange_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5234375, "outer_face_mass_present": false}, "fine_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "fine_step_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5625, "outer_face_mass_present": false}, "middle_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 1.1642856634157104e-11, "maximum_r99_to_half_width": 0.5, "outer_face_mass_present": true}, "pde_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 8.178909993925213e-10, "maximum_r99_to_half_width": 0.4921875, "outer_face_mass_present": true}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.
## refinement_middle_exchange

Column ABM denotes the declared reference; PDE denotes the candidate in this table.

| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | total_mass | 1000 | 1000 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r_K_ratio | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | active_fraction | 1 | 1 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r50 | 13 | 13 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r90 | 17 | 17 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | r99 | 19 | 19 | 0 | 0 | null | null | 1 | 0 | 0 | PASS |
| 0 | radial_profile_L2 | null | null | 0 | 0.0013341084706137063 | 0 | -9.2592970732966027 | 1.3618988231855223e-07 | 0 | 0 | PASS |
| 4 | total_mass | 961.02174565544249 | 960.77865103107206 | 0.00025295434309304096 | 0.00048335018898518136 | 0.52333556261585756 | -0.48389212733930609 | 0.63545085497481313 | 3.4250209649226806e-27 | 1.0194027028186709e-26 | PASS |
| 4 | r_K_ratio | 0.94655770730419042 | 0.95048738605755023 | 0.0041515469400715344 | 0.00094404091210744344 | 4.3976345588707257 | 3.1966521994991681 | 0.0060048648720942869 | 1.2711955273709415e-21 | 1.8544854677319369e-20 | PASS |
| 4 | active_fraction | 0.30851441942057212 | 0.30991544394758685 | 0.0014010245270147208 | 0.0056739961346165477 | 0.24692024699615042 | 2.8578368354700578 | 0.011975087790423776 | 1.5917912866218847e-17 | 1.283279847253759e-16 | PASS |
| 4 | r50 | 13.1875 | 13 | 0.014218009478672985 | 0 | null | -1.8605210188381269 | 0.082530706346961857 | 1.2728796844468052e-09 | 3.8879211368783638e-10 | PASS |
| 4 | r90 | 19.0625 | 19.0625 | 0 | 0.013114754098360656 | 0 | 0 | 1 | 1.8943201313454503e-13 | 3.4493405137426388e-12 | PASS |
| 4 | r99 | 42.25 | 42.0625 | 0.0044378698224852072 | 0.0073964497041420114 | 0.59999999999999998 | -1.3789156793307651 | 0.18813706218300907 | 8.0903538793302976e-15 | 2.5416876981087546e-15 | PASS |
| 4 | radial_profile_L2 | null | null | 0.028538680656523104 | 0.014407910313865093 | 1.9807647351232895 | 0.32675015347026881 | 0.74837179348294935 | 8.4061264079671283e-19 | 2.6018989236188797e-16 | PASS |
| 8 | total_mass | 928.51397835569628 | 929.13141548265355 | 0.00066497343211859407 | 0.00177604905004365 | 0.37441163694338908 | 0.65769528171628688 | 0.52069994583970358 | 6.0593745709895543e-23 | 1.9080145278164413e-22 | PASS |
| 8 | r_K_ratio | 0.88434265759862962 | 0.89279820228143036 | 0.0095613895927638536 | 0.0014968649513863032 | 6.3876100405107952 | 3.7910690523591271 | 0.0017753539416338027 | 3.4191881813004347e-18 | 1.8987856581597565e-15 | PASS |
| 8 | active_fraction | 0.31191934111235159 | 0.31541176626189538 | 0.003492425149543784 | 0.0061020419805564753 | 0.57233712266681136 | 2.7070784051008907 | 0.016226770786123425 | 6.1114123894083156e-12 | 8.9465309362295577e-10 | PASS |
| 8 | r50 | 14 | 13.75 | 0.017857142857142856 | 0 | null | -2.2360679774997898 | 0.040968955955836141 | 1.7220532635300656e-08 | 1.2208494517465023e-10 | PASS |
| 8 | r90 | 20.3125 | 20.4375 | 0.0061538461538461538 | 0.046153846153846156 | 0.13333333333333333 | 0.56493268286603204 | 0.58047096825964761 | 3.3699549768413853e-08 | 1.8801386891978101e-07 | PASS |
| 8 | r99 | 62.4375 | 62.6875 | 0.004004004004004004 | 0.023023023023023025 | 0.17391304347826086 | 0.49588470368046472 | 0.62716246544855792 | 1.1479593455955577e-09 | 2.291409801050083e-09 | PASS |
| 8 | radial_profile_L2 | null | null | 0.03318565867017554 | 0.012294047784762041 | 2.6993272883897346 | 0.16009792446867904 | 0.87493997104187671 | 1.1918539787579317e-17 | 9.4128144470295766e-15 | PASS |

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5546875, "outer_face_mass_present": false}, "coarse_exchange": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "coarse_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "fine_exchange": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.53125, "outer_face_mass_present": false}, "fine_exchange_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5234375, "outer_face_mass_present": false}, "fine_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "fine_step_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5625, "outer_face_mass_present": false}, "middle_step": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 1.1642856634157104e-11, "maximum_r99_to_half_width": 0.5, "outer_face_mass_present": true}, "pde_baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 8.178909993925213e-10, "maximum_r99_to_half_width": 0.4921875, "outer_face_mass_present": true}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

## Mixed coverage

| Seed | Scenario | Min ABM fraction | Min active fraction | Max PDE fraction | To PDE | To ABM | Result |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | --- |
| 1 | coarse_step | 0.7119331966809771 | 0.33542976939203356 | 0.28806680331902301 | 720 | 455 | PASS |
| 1 | middle_step | 0.72305072228221978 | 0.33542976939203356 | 0.27694927771778027 | 744 | 495 | PASS |
| 1 | fine_step | 0.73643876990252777 | 0.33542976939203356 | 0.26356123009747223 | 770 | 530 | PASS |
| 1 | coarse_exchange | 0.69540664597018786 | 0.33600000000000002 | 0.30459335402981208 | 474 | 195 | PASS |
| 1 | fine_exchange | 0.72409885900702853 | 0.33542976939203356 | 0.27590114099297147 | 1274 | 1033 | PASS |
| 17 | fine_step_baseline | 0.68955618039502964 | 0.29059829059829062 | 0.31044381960497047 | 843 | 553 | PASS |
| 17 | fine_exchange_baseline | 0.68108501209212435 | 0.28755364806866951 | 0.31891498790787559 | 1401 | 1105 | PASS |
| 2 | coarse_step | 0.6894273761005022 | 0.31180400890868598 | 0.31057262389949786 | 820 | 526 | PASS |
| 2 | middle_step | 0.67970674761848948 | 0.31096196868008946 | 0.32029325238151063 | 789 | 486 | PASS |
| 2 | fine_step | 0.72807415225173011 | 0.31343283582089554 | 0.27192584774826994 | 802 | 551 | PASS |
| 2 | coarse_exchange | 0.69881278280156778 | 0.31469979296066253 | 0.30118721719843228 | 487 | 213 | PASS |
| 2 | fine_exchange | 0.73970676112421674 | 0.31034482758620691 | 0.26029323887578321 | 1291 | 1061 | PASS |
| 18 | fine_step_baseline | 0.7208929645010147 | 0.32048681541582152 | 0.2791070354989853 | 770 | 507 | PASS |
| 18 | fine_exchange_baseline | 0.71016468481846262 | 0.31967213114754101 | 0.28983531518153732 | 1357 | 1084 | PASS |
| 3 | coarse_step | 0.72647831606144708 | 0.30648769574944074 | 0.27352168393855292 | 777 | 521 | PASS |
| 3 | middle_step | 0.74865317895289385 | 0.30612244897959184 | 0.25134682104710615 | 757 | 527 | PASS |
| 3 | fine_step | 0.72489901086705855 | 0.30691056910569103 | 0.27510098913294145 | 803 | 565 | PASS |
| 3 | coarse_exchange | 0.72168837163215993 | 0.30549450549450552 | 0.27831162836784018 | 498 | 237 | PASS |
| 3 | fine_exchange | 0.72902215479196497 | 0.30365296803652969 | 0.27097784520803497 | 1323 | 1070 | PASS |
| 19 | fine_step_baseline | 0.66715387655707403 | 0.27235772357723576 | 0.33284612344292597 | 762 | 460 | PASS |
| 19 | fine_exchange_baseline | 0.68776255067884839 | 0.27235772357723576 | 0.31223744932115161 | 1415 | 1150 | PASS |
| 4 | coarse_step | 0.69445344545089649 | 0.28297872340425534 | 0.30554655454910351 | 815 | 553 | PASS |
| 4 | middle_step | 0.7109530628504056 | 0.27991452991452992 | 0.28904693714959429 | 838 | 582 | PASS |
| 4 | fine_step | 0.69947001021145472 | 0.27991452991452992 | 0.30052998978854528 | 834 | 557 | PASS |
| 4 | coarse_exchange | 0.66715617219356016 | 0.28389830508474578 | 0.33284382780643984 | 478 | 163 | PASS |
| 4 | fine_exchange | 0.71614787084699283 | 0.28237791932059447 | 0.28385212915300712 | 1308 | 1051 | PASS |
| 20 | fine_step_baseline | 0.68990079283106387 | 0.31910569105691056 | 0.31009920716893619 | 772 | 480 | PASS |
| 20 | fine_exchange_baseline | 0.670120301612535 | 0.31910569105691056 | 0.329879698387465 | 1284 | 1004 | PASS |
| 5 | coarse_step | 0.68445489795620706 | 0.24675324675324675 | 0.31554510204379305 | 794 | 497 | PASS |
| 5 | middle_step | 0.69664207848771198 | 0.24481327800829875 | 0.30335792151228802 | 817 | 567 | PASS |
| 5 | fine_step | 0.71389124194200382 | 0.24458874458874458 | 0.28610875805799624 | 851 | 592 | PASS |
| 5 | coarse_exchange | 0.64005395352190253 | 0.24078091106290672 | 0.35994604647809747 | 516 | 211 | PASS |
| 5 | fine_exchange | 0.71316449942478422 | 0.23747276688453159 | 0.28683550057521578 | 1410 | 1141 | PASS |
| 21 | fine_step_baseline | 0.6781581599277704 | 0.3215767634854772 | 0.32184184007222966 | 782 | 482 | PASS |
| 21 | fine_exchange_baseline | 0.72195018415270551 | 0.32193158953722334 | 0.27804981584729443 | 1326 | 1065 | PASS |
| 6 | coarse_step | 0.70293783495107021 | 0.31707317073170732 | 0.29706216504892985 | 802 | 525 | PASS |
| 6 | middle_step | 0.71454775072935262 | 0.31707317073170732 | 0.2854522492706475 | 780 | 520 | PASS |
| 6 | fine_step | 0.70945064766979971 | 0.31707317073170732 | 0.29054935233020029 | 803 | 525 | PASS |
| 6 | coarse_exchange | 0.71977086164783655 | 0.31707317073170732 | 0.2802291383521634 | 502 | 245 | PASS |
| 6 | fine_exchange | 0.71331291115785322 | 0.31568228105906315 | 0.28668708884214683 | 1288 | 1030 | PASS |
| 22 | fine_step_baseline | 0.69552126190362973 | 0.34200000000000003 | 0.30447873809637027 | 716 | 432 | PASS |
| 22 | fine_exchange_baseline | 0.73914663982046014 | 0.34200000000000003 | 0.26085336017953986 | 1309 | 1065 | PASS |
| 7 | coarse_step | 0.71305935015047117 | 0.312 | 0.28694064984952883 | 785 | 529 | PASS |
| 7 | middle_step | 0.73170930293612702 | 0.312 | 0.26829069706387304 | 764 | 530 | PASS |
| 7 | fine_step | 0.73298654317890066 | 0.312 | 0.26701345682109939 | 798 | 546 | PASS |
| 7 | coarse_exchange | 0.71694513030622165 | 0.312 | 0.28305486969377835 | 507 | 239 | PASS |
| 7 | fine_exchange | 0.72033938852917911 | 0.312 | 0.27966061147082094 | 1275 | 1012 | PASS |
| 23 | fine_step_baseline | 0.71114638818104059 | 0.31 | 0.28885361181895941 | 816 | 545 | PASS |
| 23 | fine_exchange_baseline | 0.7130908112768517 | 0.31 | 0.2869091887231483 | 1361 | 1110 | PASS |
| 8 | coarse_step | 0.7520702941829468 | 0.30443548387096775 | 0.2479297058170532 | 807 | 576 | PASS |
| 8 | middle_step | 0.73419858040393182 | 0.30443548387096775 | 0.26580141959606823 | 809 | 560 | PASS |
| 8 | fine_step | 0.72350416511837068 | 0.30443548387096775 | 0.27649583488162938 | 823 | 563 | PASS |
| 8 | coarse_exchange | 0.67783760455793551 | 0.30443548387096775 | 0.32216239544206449 | 473 | 173 | PASS |
| 8 | fine_exchange | 0.71681921553092343 | 0.30379746835443039 | 0.28318078446907652 | 1238 | 1007 | PASS |
| 24 | fine_step_baseline | 0.74276420250679909 | 0.32200000000000001 | 0.25723579749320102 | 765 | 537 | PASS |
| 24 | fine_exchange_baseline | 0.74699245386224677 | 0.32200000000000001 | 0.25300754613775323 | 1289 | 1056 | PASS |
| 9 | coarse_step | 0.70625497232520762 | 0.30346232179226068 | 0.29374502767479238 | 786 | 511 | PASS |
| 9 | middle_step | 0.68747579990184848 | 0.30346232179226068 | 0.31252420009815152 | 802 | 507 | PASS |
| 9 | fine_step | 0.71410402873832235 | 0.30346232179226068 | 0.28589597126167759 | 772 | 505 | PASS |
| 9 | coarse_exchange | 0.7075222449596511 | 0.30346232179226068 | 0.29247775504034879 | 479 | 204 | PASS |
| 9 | fine_exchange | 0.70507558692104433 | 0.30346232179226068 | 0.29492441307895573 | 1339 | 1060 | PASS |
| 25 | fine_step_baseline | 0.72638001469529812 | 0.33054393305439328 | 0.27361998530470177 | 804 | 564 | PASS |
| 25 | fine_exchange_baseline | 0.74040039647511746 | 0.32456140350877194 | 0.25959960352488254 | 1309 | 1071 | PASS |
| 10 | coarse_step | 0.73947208161486022 | 0.31364562118126271 | 0.26052791838513978 | 726 | 484 | PASS |
| 10 | middle_step | 0.71151051025902012 | 0.31364562118126271 | 0.28848948974097982 | 746 | 477 | PASS |
| 10 | fine_step | 0.71564808310256445 | 0.31364562118126271 | 0.2843519168974355 | 770 | 506 | PASS |
| 10 | coarse_exchange | 0.69593505624060081 | 0.31364562118126271 | 0.30406494375939919 | 467 | 184 | PASS |
| 10 | fine_exchange | 0.7620865928139432 | 0.3122448979591837 | 0.23791340718605683 | 1265 | 1060 | PASS |
| 26 | fine_step_baseline | 0.70858757869490807 | 0.29065040650406504 | 0.29141242130509193 | 810 | 544 | PASS |
| 26 | fine_exchange_baseline | 0.73501737597249439 | 0.29065040650406504 | 0.26498262402750555 | 1316 | 1065 | PASS |
| 11 | coarse_step | 0.71672761596497581 | 0.29799999999999999 | 0.28327238403502419 | 786 | 521 | PASS |
| 11 | middle_step | 0.69984624413858476 | 0.29662921348314608 | 0.30015375586141518 | 806 | 526 | PASS |
| 11 | fine_step | 0.70670199942187228 | 0.29782608695652174 | 0.29329800057812772 | 784 | 510 | PASS |
| 11 | coarse_exchange | 0.71825300981941431 | 0.29799999999999999 | 0.28174699018058569 | 493 | 228 | PASS |
| 11 | fine_exchange | 0.73871365978230374 | 0.29791666666666666 | 0.26128634021769626 | 1306 | 1068 | PASS |
| 27 | fine_step_baseline | 0.69287436739805452 | 0.28222222222222221 | 0.30712563260194553 | 745 | 460 | PASS |
| 27 | fine_exchange_baseline | 0.73155620837448787 | 0.28634361233480177 | 0.26844379162551218 | 1274 | 1021 | PASS |
| 12 | coarse_step | 0.73508990612619163 | 0.318 | 0.26491009387380837 | 792 | 545 | PASS |
| 12 | middle_step | 0.72957707696430951 | 0.314 | 0.27042292303569049 | 778 | 537 | PASS |
| 12 | fine_step | 0.73636321440930541 | 0.314 | 0.26363678559069459 | 787 | 560 | PASS |
| 12 | coarse_exchange | 0.68947288512740001 | 0.314 | 0.31052711487259993 | 485 | 194 | PASS |
| 12 | fine_exchange | 0.72232615694505398 | 0.314 | 0.27767384305494597 | 1337 | 1077 | PASS |
| 28 | fine_step_baseline | 0.70855198769387395 | 0.29295154185022027 | 0.29144801230612616 | 843 | 569 | PASS |
| 28 | fine_exchange_baseline | 0.68694366729791956 | 0.29424778761061948 | 0.3130563327020805 | 1350 | 1071 | PASS |
| 13 | coarse_step | 0.70403856404948273 | 0.30208333333333331 | 0.29596143595051727 | 818 | 539 | PASS |
| 13 | middle_step | 0.70721944935603964 | 0.30208333333333331 | 0.29278055064396036 | 799 | 523 | PASS |
| 13 | fine_step | 0.73238404209643426 | 0.30208333333333331 | 0.26761595790356574 | 781 | 536 | PASS |
| 13 | coarse_exchange | 0.67786451368469913 | 0.30208333333333331 | 0.32213548631530087 | 529 | 229 | PASS |
| 13 | fine_exchange | 0.74747425824586011 | 0.30208333333333331 | 0.25252574175413994 | 1374 | 1153 | PASS |
| 29 | fine_step_baseline | 0.71106646821852759 | 0.33333333333333331 | 0.28893353178147241 | 749 | 477 | PASS |
| 29 | fine_exchange_baseline | 0.7420604396954944 | 0.33541666666666664 | 0.25793956030450566 | 1308 | 1069 | PASS |
| 14 | coarse_step | 0.72023109654977524 | 0.28199999999999997 | 0.27976890345022465 | 849 | 585 | PASS |
| 14 | middle_step | 0.71224955532529977 | 0.28199999999999997 | 0.28775044467470023 | 858 | 610 | PASS |
| 14 | fine_step | 0.72514549846413401 | 0.28199999999999997 | 0.27485450153586599 | 828 | 596 | PASS |
| 14 | coarse_exchange | 0.70656113220176597 | 0.28043478260869564 | 0.29343886779823397 | 504 | 228 | PASS |
| 14 | fine_exchange | 0.71585369575278601 | 0.28043478260869564 | 0.28414630424721393 | 1423 | 1154 | PASS |
| 30 | fine_step_baseline | 0.70267594645875586 | 0.30379746835443039 | 0.29732405354124419 | 770 | 489 | PASS |
| 30 | fine_exchange_baseline | 0.72126065349981427 | 0.3037190082644628 | 0.27873934650018573 | 1287 | 1036 | PASS |
| 15 | coarse_step | 0.73127336583756219 | 0.30364372469635625 | 0.26872663416243786 | 788 | 537 | PASS |
| 15 | middle_step | 0.74816480782142425 | 0.30364372469635625 | 0.25183519217857575 | 757 | 531 | PASS |
| 15 | fine_step | 0.70421909307572861 | 0.30364372469635625 | 0.29578090692427145 | 770 | 494 | PASS |
| 15 | coarse_exchange | 0.66172321980384141 | 0.30364372469635625 | 0.33827678019615864 | 502 | 188 | PASS |
| 15 | fine_exchange | 0.71807688194332941 | 0.30364372469635625 | 0.28192311805667064 | 1324 | 1084 | PASS |
| 31 | fine_step_baseline | 0.70233314733111929 | 0.32016632016632018 | 0.29766685266888071 | 773 | 496 | PASS |
| 31 | fine_exchange_baseline | 0.7297100799767281 | 0.3165137614678899 | 0.2702899200232719 | 1198 | 951 | PASS |
| 16 | coarse_step | 0.67147062097997889 | 0.32924335378323111 | 0.328529379020021 | 778 | 487 | PASS |
| 16 | middle_step | 0.69929156142878091 | 0.32924335378323111 | 0.30070843857121904 | 821 | 537 | PASS |
| 16 | fine_step | 0.69621534058131995 | 0.32924335378323111 | 0.30378465941867999 | 804 | 519 | PASS |
| 16 | coarse_exchange | 0.65762819695171448 | 0.328125 | 0.34237180304828557 | 518 | 197 | PASS |
| 16 | fine_exchange | 0.72438800199640874 | 0.3273542600896861 | 0.27561199800359121 | 1279 | 1024 | PASS |
| 32 | fine_step_baseline | 0.67719594154962093 | 0.27676767676767677 | 0.32280405845037902 | 824 | 520 | PASS |
| 32 | fine_exchange_baseline | 0.68438814545270754 | 0.27676767676767677 | 0.3156118545472924 | 1400 | 1103 | PASS |

## Retained failure and causal correction

The first ensemble and the transport-mask-only pilot both failed the 8-hour
primary r90 TOST: ABM mean 22.9375, hybrid mean 18, relative error
0.21525885558583105; lower-side p=0.09728578357367412 at the unchanged
0.25 margin. Reflection alone did not remove this bias. Inspection found
that the external obstacle charged both the replaced agent and the incoming
agent before native collision/swap resolution. This disabled crowding swaps
and self-overlapping large-cell moves. Correcting that v4 interface permits
the native mechanism and is tested directly; the failed datasets remain
in the numeric report, rather than being relabeled as equivalence.

## Remaining limits

This experiment covers early invasion and mixed representation, not a
long-time growing vascular tumor. All current r cells remain agents, which
is conservative for representation but may limit core acceleration. K
mass below one whole cell or without a full free footprint can remain
outside the core; `front_pde_mass` quantifies it. The result does not
assert zero PDE mass throughout the front. Whole-cell/fractional conversion
and fixed 16-voxel regional redistribution are interface approximations.

Large opt-in phase B ensembles remain unexecuted. Phase E 3D/scale work
and phase F production configurations/refactoring are separate follow-ups.
