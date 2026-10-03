# Phase F report

Structured PDE activation, transport, exchange, reaction, nutrient and
vascular kernels now live in separate source files. Six immutable policy
objects select their contracts once from the schema/model configuration.
Kernel expressions and traversal order are preserved. Full native
comparisons verify the maintenance split independently of newer fixes.

Resource-limited production checks found an actual old checkpoint omission:
positive density tails outside traversal bounds, active totals and nutrient
work/front buffers were not all persisted. New structured schema 17
stores and hashes the complete state. New hybrid v5 initializes resource
fields and individual clocks through standalone ABM assembly before core
classification, and restores without recanonicalizing saved front state.
Earlier schemas retain their formats and arithmetic. The failed schema-16
cross-process measurements remain in the numeric report.

## Tests

| Build | Passed CTests | Validation suites | Full-suite seconds |
| --- | ---: | ---: | ---: |
| release | 65/65 | 13 | 4895.79 |
| hdf5 | 68/68 | 13 | 4914.2399999999998 |

New assert tests cover 2D/3D ABM/hybrid resource and clock initialization;
separate-process dynamic-vascular continuation from one to eight threads
for ABM/PDE/hybrid; recommended graph dry-runs, registry paths and distinct
output namespaces. After selecting persistence by the PDE policy, an
additional restart check covers a migrated v4 wrapper with schema 17.
Both builds pass this expanded separate-process test; the original
full suite covers the unchanged existing paths. All C++ tests use
`-UNDEBUG`. Compiler/linker logs have
no new warnings. The 416 published, 96 schema-14, 48 schema-15 and 176
hybrid-v4 native reports remain exact against preceding fixtures in both
builds. The 48 v5 native realization reports also match between builds.

## Resource-initialized hybrid equivalence

This uses the previously declared phase B margins and phase D coverage
guards without widening them. The case is 256-square, r20, eight hours
with 16 paired hybrid/ABM seeds and 16 disjoint ABM baseline seeds.

# Resource-initialized hybrid v5

# Shared ABM/PDE statistical validation

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

Boundary diagnostics: {"abm": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}, "baseline": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.5546875, "outer_face_mass_present": false}, "pde": {"annotated": true, "boundary_influenced": false, "maximum_boundary_mass": 0, "maximum_r99_to_half_width": 0.546875, "outer_face_mass_present": false}}

A failed equivalence test is not proof of inequivalence. Review the fixed margin, confidence interval and coverage. The profile test is a nonlinear jackknife approximation.
Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json.

## Production workload and raw measurements

The three configurations define 2000-square, r200, transient nutrient
and shared VEGF angiogenesis for a configured 2160 hours. The recorded
experiment is a one-hour prefix and a half-hour checkpoint, followed
by a separate process restored with eight threads. It preserves the
configured biological horizon. All restart state/field hashes and
diagnostics match exactly; the complete 2160-hour runs are unexecuted.

| Model | Scope | Threads | Time h | Seconds | Resident bytes | Checkpoint bytes | Mass | Active mass | r99 | Vascular roots | Centerline growth | State checksum | Field checksum |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| abm | reference | 1 | 1 | 0.65745583399999996 | 796393472 | 0 | 905 | 46 | 18 | 94 | 148 | 9229937819435962344 | 15254145570685839155 |
| abm | checkpoint_prefix | 1 | 0.5 | 0.877671333 | 799506432 | 228169845 | 990 | 42 | 18 | 55 | 42 | 11311091147781930892 | 3591242577413171105 |
| abm | resumed | 8 | 1 | 0.45349208299999999 | 773865472 | 0 | 905 | 46 | 18 | 94 | 148 | 9229937819435962344 | 15254145570685839155 |
| pde | reference | 1 | 1 | 1598.1177935840001 | 8387461120 | 0 | 971.02563189167199 | 153.90012883284106 | 73 | 92.840060476975253 | 155.78218978199408 | 1253543682994000429 | 1289385624660310941 |
| pde | checkpoint_prefix | 1 | 0.5 | 405.01444720799998 | 4942020608 | 2001976223 | 995.34471763408965 | 155.43032215426302 | 47 | 48.617923395244908 | 38.921995426016224 | 16911011989152890872 | 16064382902029045039 |
| pde | resumed | 8 | 1 | 1199.6755417500001 | 8939520000 | 0 | 971.02563189167199 | 153.90012883284106 | 73 | 92.840060476975253 | 155.78218978199408 | 1253543682994000429 | 1289385624660310941 |
| hybrid | reference | 1 | 1 | 2.6449263749999998 | 2724380672 | 0 | 984.47457234368994 | 44 | 18 | 99.774714975153856 | 164.18976556558115 | 1098749683618659895 | 674490213671663169 |
| hybrid | checkpoint_prefix | 1 | 0.5 | 4.1789416250000002 | 2714533888 | 1236775335 | 998.73465277777746 | 43 | 18 | 49.924952469897605 | 39.411541329225244 | 8496880339641421625 | 1718647749645970707 |
| hybrid | resumed | 8 | 1 | 2.8618167909999999 | 2689449984 | 0 | 984.47457234368994 | 44 | 18 | 99.774714975153856 | 164.18976556558115 | 1098749683618659895 | 674490213671663169 |

Timing includes construction or restore, coupled initialization, operators,
diagnostics, optional checkpoint writing and teardown. Peak resident
memory includes derived caches and temporary seeding. Different scopes
are reported separately without a speedup claim. The common all-cell
unit-bin r99 convention and unchanged interior guard are documented in
[the production protocol](production_benchmark_contract.md).

The one-seed production endpoints differ materially across models: mass
905 / 971.02563189167199 / 984.47457234368994, active mass
46 / 153.90012883284106 / 44, and r99 18 / 73 / 18 for ABM/PDE/hybrid.
The common radial convention makes this a measurable closure discrepancy
that requires a resource-limited r200 ensemble and mechanism study. Exact
continuation does not establish equivalence between models. Peak PDE
reference memory is 8387461120 bytes and restored-process memory reaches
8939520000 bytes, so this workload has not established a decimal 8-GB bound.

## Registry and delivery

Four current recommended graphs were generated by `atcg_config_migrate`,
dry-run through `atcg_sim`, and registered with all other repository
configurations as recommended dependencies, legacy, reproduction or
benchmark inputs. Recommendations are qualified by experiment and
evidence. The migration and deprecation policy is in
[the registry](recommended_configurations.md). A reviewed PR #67
description is saved locally in `pr67_description.md`; the user's
local-only instruction prevents pushing or publishing that update.

## Remaining limits

The short production prefixes do not establish 2160-hour feasibility,
boundary independence or statistical equivalence. Long ensembles remain
opt-in and unexecuted. Growing-tumor vascular equivalence and long-time
3D validation are still needed. ABM endpoint/PDE starting-cell vascular
source sampling also needs a growing-tumor splitting study. PDE clock
and space closures and fractional
hybrid interface correlations remain approximations. Native checkpoints
require a compatible ABI. No new Linux/macOS hosted execution is claimed:
new checks ran locally on macOS, and previously active cloud jobs were
cancelled after the local-only instruction. See the numeric report for
all native realization values, statistical tables and retained failures.
