# Phase E report

The new ABM guidance model and structured schema 16 prepare 3D directional
sector means from signed summed-volume corner kernels. Tiled spectral
convolution prepares all requested means at resource refresh; migration
queries have constant cost independent of the 70-voxel cone size.
Published sums remain selected by default. Physical exterior sites do
not enter cone counts. Sparse 3D population storage accepts dynamic
vasculature, and empty division-work channels remain unallocated.

## Verification

| Build | Passed CTests | Validation suites | Full-suite seconds |
| --- | ---: | ---: | ---: |
| release | 62/62 | 12 | 5047.04 |
| hdf5 | 65/65 | 12 | 5046.4499999999998 |

The new assert tests check direct directional means and exact counts,
uniform and gradient fields, physical corners, missing-envelope rejection,
thread determinism, 3D dense/sparse biological equality and continuation.
The 416 published, 96 schema-14, 48 schema-15 vascular and 176 hybrid-v4
native reports remain exactly equal to their preceding fixtures in both
builds. Checked build logs have no new compiler or linker warnings.

The fixed 5e-4 precision budget was declared in
[the sector protocol](three_dimensional_sector_contract.md) before these
measurements; it changes no statistical equivalence margin.

| Thin | Edge | Extent | Maximum direct-mean error | Cache allocated bytes |
| --- | ---: | ---: | ---: | ---: |
| True | 8 | 16 | 1.2989609388114332e-14 | 206080 |
| False | 8 | 16 | 1.6181500583911657e-13 | 16495360 |
| False | 70 | 17 | 6.6252558994506217e-13 | 123398072 |
| False | 70 | 128 | 1.2086454059812013e-10 | 1095239680 |

## Isolated 256-cube coupled benchmark

The measured case is r20, 256 r plus 256 K cells, transient nutrient
and dynamic VEGF vasculature for one hour. Construction, coupled
initialization, cache preparation, all operators, diagnostics and
teardown are included. Each repetition is a separate process. Test
suites finished before timing. This measures performance/determinism
over four resource updates, rather than long-time model equivalence.

| Model | Threads | Median seconds | Speedup | Maximum resident bytes |
| --- | ---: | ---: | ---: | ---: |
| abm | 1 | 82.891666458000003 | 1 | 4435918848 |
| abm | 4 | 24.356374708000001 | 3.4032842511153243 | 4738252800 |
| abm | 8 | 13.043498458 | 6.3550179213736833 | 5141200896 |
| pde | 1 | 1065.0316244579999 | 1 | 13064814592 |
| pde | 4 | 677.23256537500004 | 1.5726231709904059 | 13369917440 |
| pde | 8 | 602.86270091599999 | 1.7666238479172329 | 13775110144 |

Every biological diagnostic and checksum is exact across thread counts
and repetitions within each model. The two models have different endpoint
mass (ABM 475; PDE 511.99999007859481), active mass (146;
79.319266422425045) and centerline growth (223; 278.57877937237703).
These single-trajectory differences require a separate 3D ensemble and
mechanism study before an equivalence claim. The raw measurements are:

| Model | Threads | Repetition | Seconds | Resident bytes | Mass | Active mass | Vascular roots | Centerline growth | State checksum | Field checksum |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| abm | 1 | 1 | 82.814966624999997 | 4435902464 | 475 | 146 | 92 | 223 | 17185381434236664832 | 9486803310450670894 |
| abm | 4 | 1 | 24.411438959000002 | 4738220032 | 475 | 146 | 92 | 223 | 17185381434236664832 | 9486803310450670894 |
| abm | 8 | 1 | 13.043498458 | 5141200896 | 475 | 146 | 92 | 223 | 17185381434236664832 | 9486803310450670894 |
| abm | 1 | 2 | 82.891666458000003 | 4435836928 | 475 | 146 | 92 | 223 | 17185381434236664832 | 9486803310450670894 |
| abm | 4 | 2 | 24.356374708000001 | 4738220032 | 475 | 146 | 92 | 223 | 17185381434236664832 | 9486803310450670894 |
| abm | 8 | 2 | 12.947250709 | 5141184512 | 475 | 146 | 92 | 223 | 17185381434236664832 | 9486803310450670894 |
| abm | 1 | 3 | 83.197777416999998 | 4435918848 | 475 | 146 | 92 | 223 | 17185381434236664832 | 9486803310450670894 |
| abm | 4 | 3 | 23.920976875000001 | 4738252800 | 475 | 146 | 92 | 223 | 17185381434236664832 | 9486803310450670894 |
| abm | 8 | 3 | 13.077173791 | 5141200896 | 475 | 146 | 92 | 223 | 17185381434236664832 | 9486803310450670894 |
| pde | 1 | 1 | 1066.051783333 | 13064650752 | 511.99999007859481 | 79.319266422425045 | 100 | 278.57877937237703 | 6924990288887117941 | 2700684879161137722 |
| pde | 4 | 1 | 682.16842920800002 | 13367066624 | 511.99999007859481 | 79.319266422425045 | 100 | 278.57877937237703 | 6924990288887117941 | 2700684879161137722 |
| pde | 8 | 1 | 604.28722274999996 | 13775110144 | 511.99999007859481 | 79.319266422425045 | 100 | 278.57877937237703 | 6924990288887117941 | 2700684879161137722 |
| pde | 1 | 2 | 1065.0316244579999 | 13064814592 | 511.99999007859481 | 79.319266422425045 | 100 | 278.57877937237703 | 6924990288887117941 | 2700684879161137722 |
| pde | 4 | 2 | 672.51482266699998 | 13369917440 | 511.99999007859481 | 79.319266422425045 | 100 | 278.57877937237703 | 6924990288887117941 | 2700684879161137722 |
| pde | 8 | 2 | 601.868051208 | 13767376896 | 511.99999007859481 | 79.319266422425045 | 100 | 278.57877937237703 | 6924990288887117941 | 2700684879161137722 |
| pde | 1 | 3 | 1062.2036048750001 | 13063979008 | 511.99999007859481 | 79.319266422425045 | 100 | 278.57877937237703 | 6924990288887117941 | 2700684879161137722 |
| pde | 4 | 3 | 677.23256537500004 | 13364789248 | 511.99999007859481 | 79.319266422425045 | 100 | 278.57877937237703 | 6924990288887117941 | 2700684879161137722 |
| pde | 8 | 3 | 602.86270091599999 | 13771063296 | 511.99999007859481 | 79.319266422425045 | 100 | 278.57877937237703 | 6924990288887117941 | 2700684879161137722 |

## Remaining limits

Resource and vascular arrays remain dense. The new storage reduces
population/cohort memory and empty channels; it does not promise that
a fully occupied 3D domain fits below 8 GB. FFT caches add memory and
refresh cost, and speedup is measured without a required target. The
new sums have floating-point error and therefore use new model strings.
Long-time 3D biological equivalence and production-horizon boundary
independence remain unmeasured. All raw data and source digests are in
[the numeric report](phase_e_numeric_results.json).
