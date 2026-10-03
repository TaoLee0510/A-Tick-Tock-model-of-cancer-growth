# Prepared 3D sector means

This protocol precedes the phase E precision and performance measurements.
The new sector model will use the same exact 70-voxel box and directional cone
membership as the shared rule contract. Its arithmetic is opt-in; published
models retain their original sums.

Each cone is decomposed into disjoint x-row segments. Segment sums are signed
corner differences of a 3D summed-volume prefix. Precomputed corner kernels
and tiled spectral convolution materialize every requested directional mean
at resource refresh, before transport. Each query tile (up to 32 voxels on each axis) has a halo
covering the entire cone. Padding prevents periodic convolution aliases at
query sites; physical grid-exterior sites have zero membership and do not
enter the averaging denominator. Site counts must match direct enumeration
exactly, including corners and thin-layer geometry.

A migration query is a read-only tile lookup with cost independent of cone
size. Preparation cost, kernel memory and tile memory are included in reported
whole-model timings; the work is not omitted from the benchmark. Callers must
prepare a conservative transport envelope, including newborn displacement.
A query outside that envelope fails instead of silently doing a cubic scan.
Caches are derived state and are rebuilt from the restored field and cells.

The fixed direct-mean unit-test bound is 5e-4 of the vessel nutrient value.
For the 70-voxel, 45-degree, at least 128-cube corner fixtures, every nonempty
cone contains at least 32 sites. A conservative rounding budget for a
128-cube double prefix and at most eight corner contributions per each of
70 squared rows is `512*epsilon*128^3*(8*70^2)/32`, approximately 2.92e-4.
The factor 512 covers the three 128-entry prefix scans, transform stages and division.
Rounding this budget upward to 5e-4 leaves additional rounding headroom.
Smaller-window fixtures have a smaller prefix/corner conditioning budget.
Counts, directional orientation, uniform-field behavior, linear gradients,
and thread determinism are independently checked. Actual maximum precision
errors must be reported. This unit-test arithmetic bound does not alter any
of the phase B/C/D statistical equivalence margins.

## Whole-model performance protocol

The opt-in benchmark runs shared ABM and structured PDE separately on the
256-cube r20 profile for 1 hour (four 0.25-hour resource steps), with 256 r
and 256 K cells, transient nutrient and dynamic VEGF-guided vasculature.
Initial nutrient is 0.15 and the activation threshold is 0.001, so the case
exercises hypoxia and activation. Three separate processes for each model
and each of 1, 4 and 8 threads report wall time including construction,
initialization, prefix/sector preparation and all operators; PDE initialization
uses the same resource-limited cell work and activation clocks as the ABM.
The temporary seeding environment is released before PDE time integration.
Peak resident
memory includes the derived caches. Checksums and all biological diagnostics
must agree exactly across repetitions and thread counts. Vascular roots and
centerline growth must be positive at the endpoint. Timing speedup is
measured, not required to pass a chosen threshold. The 1-hour case is a
performance and determinism benchmark, not long-time equivalence evidence.

`sparse_zero_pages_3d_v2` reserves demand-zero population/work mappings in 3D
and allows dynamic vessels; it retains full dense resource and vascular
fields. Its active-region budget remains explicit. It does not promise a
fixed memory bound for a tumor that fills the entire 3D domain.
The new storage policy also keeps empty division-work channels unallocated
and uses the existing cached conservative distribution transfers. Dense
population fields and every published storage policy retain the original
channel layout; a small 3D test compares the resulting biological fields
exactly and verifies the new layout's checkpoint restoration.

To keep timing isolated, `--test-suites MANIFEST` accepts a JSON list of
CTest logs and waits for their successful completion before spawning measured
processes. A `benchmark_dependencies.json` in the report directory is used
when no explicit manifest is supplied. This local execution metadata is not
part of the biological configuration or committed benchmark results.
