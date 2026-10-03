# Phase A: versioned bug fixes

Hybrid schema 3, `hybrid_volume_coupling_v3`, corrects the large-cell
external activation denominator, uses occupied volume to decide placement,
and distributes large-agent biological consumers over the four/eight-voxel
footprint. Anchor counts remain separate from resource sinks. Native planar
large-cell z=1 footprint storage maps to the planar PDE resource voxel.
Schemas 1/2 retain their old arithmetic and native checkpoint identifiers.

Structured schema 13 records cumulative vascular deletion of normal r,
active r and K separately for small/large stages. Adaptive hybrid v3 adds
vessel-deleted agents to those counters. Metrics include all six counters and
their sum; native checkpoint format 9 preserves them and hashes them. A v13
resume refuses a metrics header lacking these columns. The operator remains
cell deletion; conservative relocation has not been implemented.

The structured executable links through the shared-rule target once. This
removes the repeated structured-core archive responsible for Apple's
`ignoring duplicate libraries` linker warning without suppressing diagnostics.
The migration tool upgrades structured rule graphs to v13 but preserves the
hybrid wrapper's explicit model selection.

## Added verification

- `atcg3d_hybrid_volume_test` checks identical large-cell activation density
  in agent and density representations in 2D and 3D, for both phenotypes;
  one-cell footprint consumption/VEGF production; volume-limited placement;
  the published v2 placement rule; single counting of vessel-deleted agents;
  and bitwise mixed restart with four threads.
- `atcg3d_vascular_removal_test` checks all six deletion counters in 2D/3D,
  loss equal to the biological mass removed, no repeated counting, metrics
  values/header continuity and rejection, checkpoint restoration with four
  threads, and the absence of counter state in schema 12.
- `atcg3d_hybrid_volume_validation` repeats the seven v2 native/mixed/refinement
  scenarios using v3 on a 256-square grid for 48 hours with 16 paired seeds.
  It uses the already declared comparison/refinement tolerances unchanged.
- CLI coverage adds loading the v3 preset and migrating structured graphs
  to v13. Published compatibility and legacy PDE bitwise fixtures remain
  required.

## Exact unit-test values

The single-large-cell comparison uses an edge-four activation window. The
same values are asserted for r and K, before and after representation exchange.

| Dimension | Window capacity | Large cell volume | Activation density | Consumer mass per footprint voxel | Total consumers |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2D | 16 | 4 | 0.25 | 0.25 | 1 |
| 3D | 64 | 8 | 0.125 | 0.125 | 1 |

The vascular-deletion tests assert the following six masses in both dimensions;
active values retain the published binary32 storage precision. Each cumulative
counter equals its corresponding initial mass, and survivors have zero mass.

| Normal r small | Normal r large | Active r small | Active r large | K small | K large | Removed total |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.050000000000000003 | 0.02 | 0.039999999105930328 | 0.0099999997764825821 | 0.029999999999999999 | 0.02 | 0.16999999888241291 |

These ensemble checks retain the prior smoke criterion. Statistical equivalence,
time-series evidence and mechanism ablation are separate phase B work. No
sampling allowance or tolerance has been enlarged for phase A.

## Full verification

Release CTest: **47/47 PASS**. HDF5 CTest: **49/49 PASS**.
Both full runs include all eight `validation` targets.
Both build logs contain no compilation or link warnings. Published
compatibility and legacy PDE bitwise regression tests pass.

The tables below retain binary64 round-trip precision. All complete
profile arrays, coverage values and refinement reports are preserved in
[phase_a_numeric_results.json](phase_a_numeric_results.json).
Raw seed realizations remain in the corresponding build validation
directories and are uploaded as CI artifacts.

### validation

Seeds: 16; endpoint: 48 hours; PASS.

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 743.3125 | 674.38546078935178 | 9.8053030358293007 | 0.092729557501923104 | 0.34999999999999998 | 0.028110788893436121 | PASS |
| r_K_ratio | 1.2542428737808808 | 1.0176043954212508 | 0.03112628731277876 | 0.1886703790042592 | 0.34999999999999998 | 0.052884588503645399 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.050000000000000003 | 0 | PASS |
| r50 | 15.8125 | 14.875 | 0.11063265039459795 | 0.059288537549407112 | 0.25 | 0.014909608094285421 | PASS |
| r90 | 24.5 | 20 | 0.18257418583505536 | 0.18367346938775511 | 0.25 | 0.015880228163857261 | PASS |
| r99 | 29.3125 | 23.5 | 0.27716947282604315 | 0.19829424307036247 | 0.29999999999999999 | 0.020150043380547478 | PASS |
| radial_profile_L2 | - | - | - | 0.32698891172682604 | 0.34999999999999998 | 0 | PASS |

Coverage (raw report values):

```json
{
  "active_fraction": {
    "abm": 0,
    "pde": 0
  }
}
```

### validation-regular

Seeds: 16; endpoint: 48 hours; PASS.

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 468.375 | 489.72804166911277 | 1.8941015470146889 | 0.045589627262583976 | 0.34999999999999998 | 0.0086177323654941049 | PASS |
| r_K_ratio | 1.4233832881828432 | 1.3285012486394163 | 0.011453123436082308 | 0.066659514925567015 | 0.34999999999999998 | 0.01714689658429951 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.050000000000000003 | 0 | PASS |
| r50 | 14.4375 | 13.1875 | 0.11180339887498948 | 0.086580086580086577 | 0.25 | 0.016502375272907537 | PASS |
| r90 | 22.375 | 19 | 0.15478479684172258 | 0.15083798882681565 | 0.25 | 0.014741738639987075 | PASS |
| r99 | 26.1875 | 22.8125 | 0.125 | 0.12887828162291171 | 0.29999999999999999 | 0.010171837708830548 | PASS |
| radial_profile_L2 | - | - | - | 0.27458671035633064 | 0.34999999999999998 | 0 | PASS |

Coverage (raw report values):

```json
{
  "active_fraction": {
    "abm": 0,
    "pde": 0
  }
}
```

### validation-active-r20

Seeds: 16; endpoint: 8 hours; PASS.

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 930.6875 | 970.37748642098722 | 2.1778784515893821 | 0.042645878902410554 | 0.34999999999999998 | 0.0049866995960910323 | PASS |
| r_K_ratio | 0.89478658768318697 | 0.91508983096094598 | 0.0037481771122353726 | 0.02269059858209194 | 0.34999999999999998 | 0.0089265591774847081 | PASS |
| active_fraction | 0.31580496549020443 | 0.27472069115986003 | 0.0024065440219043521 | 0.041084274330344395 | 0.050000000000000003 | 0.0051283453106781736 | PASS |
| r50 | 13.125 | 14 | 0.085391256382996647 | 0.066666666666666666 | 0.25 | 0.013864287036355491 | PASS |
| r90 | 22.9375 | 22.5 | 0.61892884620662714 | 0.019073569482288829 | 0.25 | 0.057501356785452748 | PASS |
| r99 | 63.5 | 57.6875 | 0.59314662324476453 | 0.091535433070866146 | 0.29999999999999999 | 0.019905440222592018 | PASS |
| radial_profile_L2 | - | - | - | 0.055546805567441175 | 0.34999999999999998 | 0 | PASS |

Coverage (raw report values):

```json
{
  "active_fraction": {
    "abm": 0.3158049654902044,
    "pde": 0.27472069115986003
  }
}
```

### validation-active-r200

Seeds: 16; endpoint: 2 hours; PASS.

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 982.375 | 989.7017335674069 | 1.4657306528003706 | 0.0074581840614906776 | 0.34999999999999998 | 0.0031795109007431883 | PASS |
| r_K_ratio | 0.97985804784240171 | 0.98185735548383302 | 0.0024945664867440129 | 0.0020404053891619075 | 0.34999999999999998 | 0.0054251952055268437 | PASS |
| active_fraction | 0.31117281779932982 | 0.29550364174646931 | 0.00085367000780809125 | 0.015669176052860501 | 0.050000000000000003 | 0.0018191707866390423 | PASS |
| r50 | 13 | 14 | 0 | 0.076923076923076927 | 0.25 | 0 | PASS |
| r90 | 33.25 | 26.6875 | 1.1760412053438718 | 0.19736842105263158 | 0.25 | 0.075372746122941078 | PASS |
| r99 | 71.5625 | 69.25 | 0.58251430597597054 | 0.032314410480349345 | 0.29999999999999999 | 0.017346207665115011 | PASS |
| radial_profile_L2 | - | - | - | 0.081326884302303726 | 0.34999999999999998 | 0 | PASS |

Coverage (raw report values):

```json
{
  "active_fraction": {
    "abm": 0.3111728177993298,
    "pde": 0.2955036417464693
  }
}
```

### validation-vascular

Seeds: 16; endpoint: 8 hours; PASS.

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| vascular_length | 26.6875 | 33.682986743131785 | 3.1983212011413391 | 0.26212596695575774 | 0.75 | 0.25538632242181519 | PASS |
| perfused_volume | 26.6875 | 33.682986743131785 | 3.1983212011413391 | 0.26212596695575774 | 0.75 | 0.25538632242181519 | PASS |
| lesion_perfused_fraction | 0.037530736013472225 | 0.090234588360324455 | 0.0039945825957257361 | 0.05270385234685223 | 0.14999999999999999 | 0.0085124555114915422 | PASS |

Coverage (raw report values):

```json
{
  "active_fraction": {
    "abm": 0,
    "pde": 0
  }
}
```

### validation-vascular-nonlinear

Seeds: 16; endpoint: 8 hours; PASS.

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| vascular_length | 26.6875 | 44.809737788726686 | 3.2142268696439849 | 0.6790534066033419 | 0.75 | 0.25665639191424189 | PASS |
| perfused_volume | 26.6875 | 44.809737788726686 | 3.2142268696439849 | 0.6790534066033419 | 0.75 | 0.25665639191424189 | PASS |
| lesion_perfused_fraction | 0.037530736013472225 | 0.12027296915844937 | 0.004029409737840264 | 0.08274223314497714 | 0.14999999999999999 | 0.0085866721513376022 | PASS |

Coverage (raw report values):

```json
{
  "active_fraction": {
    "abm": 0,
    "pde": 0
  },
  "vascular_activity": {
    "abm_roots": 58,
    "abm_anastomoses": 2,
    "pde_branching_rate": 1.5228178080433576,
    "pde_anastomosis_rate": 0.8178972356693133
  }
}
```

### validation-hybrid

Seeds: 16; endpoint: 48 hours; PASS.

Mixed representation coverage (raw report values):

```json
{
  "abm_mass": 6.0625,
  "pde_mass": 476.1682625092975,
  "to_abm": 9,
  "to_pde": 256
}
```

Native abm reference:

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 468.375 | 482.23076250929751 | 1.9805545246701453 | 0.029582626120731263 | 0.34999999999999998 | 0.009011073802128804 | PASS |
| r_K_ratio | 1.4233832881828432 | 1.2942909243521372 | 0.013278079148936769 | 0.090694028026359133 | 0.34999999999999998 | 0.019879105579852424 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.050000000000000003 | 0 | PASS |
| r50 | 14.4375 | 13.0625 | 0.125 | 0.095238095238095233 | 0.25 | 0.018450216450216449 | PASS |
| r90 | 22.375 | 19 | 0.15478479684172258 | 0.15083798882681565 | 0.25 | 0.014741738639987075 | PASS |
| r99 | 26.1875 | 25.9375 | 0.89209491273817576 | 0.0095465393794749408 | 0.29999999999999999 | 0.072593957385968591 | PASS |
| radial_profile_L2 | - | - | - | 0.29566126697015688 | 0.34999999999999998 | 0 | PASS |

Native pde reference:

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 489.72804166911277 | 482.23076250929751 | 0.66320581017523117 | 0.015309066506101438 | 0.34999999999999998 | 0.0028858702406882289 | PASS |
| r_K_ratio | 1.3285012486394163 | 1.2942909243521372 | 0.0036041601668488845 | 0.025751066716960622 | 0.34999999999999998 | 0.0057813007879525258 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.050000000000000003 | 0 | PASS |
| r50 | 13.1875 | 13.0625 | 0.085391256382996647 | 0.0094786729857819912 | 0.25 | 0.013798579514856177 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.25 | 0 | PASS |
| r99 | 22.8125 | 25.9375 | 0.88447253584645957 | 0.13698630136986301 | 0.29999999999999999 | 0.082621850910194208 | PASS |
| radial_profile_L2 | - | - | - | 0.024348409394435524 | 0.34999999999999998 | 0 | PASS |

step refinement, coarse_vs_fine:

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 486.9802930790737 | 482.23076250929751 | 0.49424241999384083 | 0.0097530241721814192 | 0.10000000000000001 | 0.0021627786831937691 | PASS |
| r_K_ratio | 1.340664755305145 | 1.2942909243521372 | 0.0025524528120272143 | 0.034590176827951825 | 0.10000000000000001 | 0.0040571492022194429 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r99 | 26.25 | 25.9375 | 0.70544755297612305 | 0.011904761904761904 | 0.10000000000000001 | 0.057268904205414015 | PASS |
| radial_profile_L2 | - | - | - | 0.0043438312838042529 | 0.14999999999999999 | 0 | PASS |

step refinement, middle_vs_fine:

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 486.9802930790737 | 485.863310751124 | 0.26169794670775864 | 0.0022936910257440927 | 0.10000000000000001 | 0.001145176370296529 | PASS |
| r_K_ratio | 1.340664755305145 | 1.3279632083777151 | 0.0012631369495516359 | 0.0094740664115831839 | 0.10000000000000001 | 0.0020077687795125748 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r99 | 26.25 | 26.375 | 0.42695628191498325 | 0.0047619047619047623 | 0.10000000000000001 | 0.034660717590888734 | PASS |
| radial_profile_L2 | - | - | - | 0.0011602356935058633 | 0.14999999999999999 | 0 | PASS |

step refinement order:

| Metric | Coarse error | Middle error | Sampling allowance | Result |
| --- | ---: | ---: | ---: | --- |
| total_mass | 0.0097530241721814192 | 0.0022936910257440927 | 0.0022898115650144776 | PASS |
| r_K_ratio | 0.034590176827951825 | 0.0094740664115831839 | 0.00414680985892969 | PASS |
| active_fraction | 0 | 0 | 0 | PASS |
| r50 | 0 | 0 | 0 | PASS |
| r90 | 0 | 0 | 0 | PASS |
| r99 | 0.011904761904761904 | 0.0047619047619047623 | 0.057388651130111393 | PASS |
| radial_profile_L2 | 0.0043438312838042529 | 0.0011602356935058633 | 0.0013372693805159147 | PASS |

exchange refinement, coarse_vs_fine:

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 485.6608560281644 | 485.85028715879059 | 0.61231255817774666 | 0.00039004817513068067 | 0.10000000000000001 | 0.002686726849159752 | PASS |
| r_K_ratio | 1.3272983106838141 | 1.3282066995729778 | 0.0028352873678673916 | 0.00068438939600222513 | 0.10000000000000001 | 0.0045521020649929248 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 18.9375 | 0.0625 | 0.0032894736842105261 | 0.10000000000000001 | 0.007009868421052631 | PASS |
| r99 | 24.6875 | 25.6875 | 0.92195444572928875 | 0.040506329113924051 | 0.10000000000000001 | 0.079582174130597025 | PASS |
| radial_profile_L2 | - | - | - | 0.0015720915166280418 | 0.14999999999999999 | 0 | PASS |

exchange refinement, middle_vs_fine:

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 485.6608560281644 | 485.863310751124 | 0.69645645116231436 | 0.00041686440331081625 | 0.10000000000000001 | 0.0030559364194276822 | PASS |
| r_K_ratio | 1.3272983106838141 | 1.3279632083777151 | 0.003351661756234164 | 0.00050094066160488032 | 0.10000000000000001 | 0.0053811499231512594 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r99 | 24.6875 | 26.375 | 1.2803604635674544 | 0.068354430379746839 | 0.10000000000000001 | 0.11051941864758461 | PASS |
| radial_profile_L2 | - | - | - | 0.0019474915248959784 | 0.14999999999999999 | 0 | PASS |

exchange refinement order:

| Metric | Coarse error | Middle error | Sampling allowance | Result |
| --- | ---: | ---: | ---: | --- |
| total_mass | 0.00039004817513068067 | 0.00041686440331081625 | 0.0028729545688062816 | PASS |
| r_K_ratio | 0.00068438939600222513 | 0.00050094066160488032 | 0.0050094686622969902 | PASS |
| active_fraction | 0 | 0 | 0 | PASS |
| r50 | 0 | 0 | 0 | PASS |
| r90 | 0.0032894736842105261 | 0 | 0.0070098684210526301 | PASS |
| r99 | 0.040506329113924051 | 0.068354430379746839 | 0.11079996974126083 | PASS |
| radial_profile_L2 | 0.0015720915166280418 | 0.0019474915248959784 | 0.0023188905601691846 | PASS |

### validation-hybrid-volume

Seeds: 16; endpoint: 48 hours; PASS.

Mixed representation coverage (raw report values):

```json
{
  "abm_mass": 5.5,
  "pde_mass": 476.8351985179054,
  "to_abm": 8.9375,
  "to_pde": 256.5
}
```

Native abm reference:

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 468.375 | 482.33519851790538 | 1.9742094932395273 | 0.029805601319253552 | 0.34999999999999998 | 0.0089822053484781041 | PASS |
| r_K_ratio | 1.4233832881828432 | 1.2947931434296158 | 0.013217386631744866 | 0.090341193282795593 | 0.34999999999999998 | 0.01978824055761301 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.050000000000000003 | 0 | PASS |
| r50 | 14.4375 | 13.0625 | 0.125 | 0.095238095238095233 | 0.25 | 0.018450216450216449 | PASS |
| r90 | 22.375 | 19 | 0.15478479684172258 | 0.15083798882681565 | 0.25 | 0.014741738639987075 | PASS |
| r99 | 26.1875 | 26.0625 | 0.89849411053532602 | 0.0047732696897374704 | 0.29999999999999999 | 0.073114690197643134 | PASS |
| radial_profile_L2 | - | - | - | 0.29532787173799768 | 0.34999999999999998 | 0 | PASS |

Native pde reference:

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 489.72804166911277 | 482.33519851790538 | 0.6657033356306552 | 0.015095813435576963 | 0.34999999999999998 | 0.0028967379596927792 | PASS |
| r_K_ratio | 1.3285012486394163 | 1.2947931434296158 | 0.0036029438539854127 | 0.025373032388432176 | 0.34999999999999998 | 0.0057793497452156726 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.050000000000000003 | 0 | PASS |
| r50 | 13.1875 | 13.0625 | 0.085391256382996647 | 0.0094786729857819912 | 0.25 | 0.013798579514856177 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.25 | 0 | PASS |
| r99 | 22.8125 | 26.0625 | 0.88741196746494244 | 0.14246575342465753 | 0.29999999999999999 | 0.082896434089547055 | PASS |
| radial_profile_L2 | - | - | - | 0.024220026542160759 | 0.34999999999999998 | 0 | PASS |

step refinement, coarse_vs_fine:

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 487.07775221479494 | 482.33519851790538 | 0.4009039965359234 | 0.0097367487538173453 | 0.10000000000000001 | 0.0017539836560658715 | PASS |
| r_K_ratio | 1.3411611329909079 | 1.2947931434296158 | 0.0020658388295670457 | 0.034573019170252385 | 0.10000000000000001 | 0.0032824561027874778 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r99 | 25.8125 | 26.0625 | 0.48733971724044822 | 0.0096852300242130755 | 0.10000000000000001 | 0.040233256656247753 | PASS |
| radial_profile_L2 | - | - | - | 0.0045289386677212573 | 0.14999999999999999 | 0 | PASS |

step refinement, middle_vs_fine:

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 487.07775221479494 | 485.85691882733988 | 0.24317279958228497 | 0.0025064445705101461 | 0.10000000000000001 | 0.0010638983890221478 | PASS |
| r_K_ratio | 1.3411611329909079 | 1.3279515421143209 | 0.0011447918097936809 | 0.0098493689920229523 | 0.10000000000000001 | 0.0018189845251703042 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r99 | 25.8125 | 26.1875 | 0.59773879468097657 | 0.014527845036319613 | 0.10000000000000001 | 0.049347462332790734 | PASS |
| radial_profile_L2 | - | - | - | 0.0013084190740619605 | 0.14999999999999999 | 0 | PASS |

step refinement order:

| Metric | Coarse error | Middle error | Sampling allowance | Result |
| --- | ---: | ---: | ---: | --- |
| total_mass | 0.0097367487538173453 | 0.0025064445705101461 | 0.0019109650735822945 | PASS |
| r_K_ratio | 0.034573019170252385 | 0.0098493689920229523 | 0.0035043289988511456 | PASS |
| active_fraction | 0 | 0 | 0 | PASS |
| r50 | 0 | 0 | 0 | PASS |
| r90 | 0 | 0 | 0 | PASS |
| r99 | 0.0096852300242130755 | 0.014527845036319613 | 0.041876101214552701 | PASS |
| radial_profile_L2 | 0.0045289386677212573 | 0.0013084190740619605 | 0.0012330365761312973 | PASS |

exchange refinement, coarse_vs_fine:

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 485.70751326384465 | 485.71918987181817 | 0.58857217199244261 | 2.4040410441796409e-05 | 0.10000000000000001 | 0.0025823098557558577 | PASS |
| r_K_ratio | 1.3275528071184484 | 1.3275877148472155 | 0.0027321546374173091 | 2.6294794888724971e-05 | 0.10000000000000001 | 0.0043856798020516022 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 18.9375 | 0.0625 | 0.0032894736842105261 | 0.10000000000000001 | 0.007009868421052631 | PASS |
| r99 | 24.625 | 25.375 | 0.7772815877574013 | 0.030456852791878174 | 0.10000000000000001 | 0.067264449279635416 | PASS |
| radial_profile_L2 | - | - | - | 0.001356710942527452 | 0.14999999999999999 | 0 | PASS |

exchange refinement, middle_vs_fine:

| Metric | ABM/reference mean | PDE/candidate mean | Paired SE | Error | Tolerance | Smoke allowance | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total_mass | 485.70751326384465 | 485.85691882733988 | 0.6742145635088862 | 0.00030760397855751588 | 0.10000000000000001 | 0.0029580584932334952 | PASS |
| r_K_ratio | 1.3275528071184484 | 1.3279515421143209 | 0.0032461253439093909 | 0.00030035339742001758 | 0.10000000000000001 | 0.0052107103165906013 | PASS |
| active_fraction | 0 | 0 | 0 | 0 | 0.02 | 0 | PASS |
| r50 | 13.0625 | 13.0625 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r90 | 19 | 19 | 0 | 0 | 0.10000000000000001 | 0 | PASS |
| r99 | 24.625 | 26.1875 | 1.0245679983941198 | 0.063451776649746189 | 0.10000000000000001 | 0.088664138256969297 | PASS |
| radial_profile_L2 | - | - | - | 0.0022897168552985488 | 0.14999999999999999 | 0 | PASS |

exchange refinement order:

| Metric | Coarse error | Middle error | Sampling allowance | Result |
| --- | ---: | ---: | ---: | --- |
| total_mass | 2.4040410441796409e-05 | 0.00030760397855751588 | 0.0026007467519915807 | PASS |
| r_K_ratio | 2.6294794888724971e-05 | 0.00030035339742001758 | 0.0045504776924638734 | PASS |
| active_fraction | 0 | 0 | 0 | PASS |
| r50 | 0 | 0 | 0 | PASS |
| r90 | 0.0032894736842105261 | 0 | 0.0070098684210526301 | PASS |
| r99 | 0.030456852791878174 | 0.063451776649746189 | 0.090837068674312224 | PASS |
| radial_profile_L2 | 0.001356710942527452 | 0.0022897168552985488 | 0.003676201764210134 | PASS |

All **416 raw realizations per build** are identical between
Release and HDF5 builds. All 304 published-model raw results also match
the preserved pre-change results from commit `7f0f45f` exactly.

## Unresolved work

Phases B-F of the current request remain outstanding. In particular, these
smoke passes do not establish statistical equivalence; vascular deletion is
a documented loss term; the existing vascular length proxy, invasive-front
classification and production/3D scale restrictions still require the later
stages. In the v3 regular-cycle case, endpoint mean ABM mass is 5.5 out of
482.3351985179054 total mass (1.1402858462123675% of ensemble mean mass); this
is not the invasive-front coverage required by phase D. The new footprint tests verify coupling arithmetic, not a claim of
joint-distribution equivalence.
