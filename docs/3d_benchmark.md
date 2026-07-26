# ATCG3D validation and scale benchmark

Last complete 10^3-through-10^7 output/checkpoint measurement: 2026-07-15.
The machine-readable report is
`benchmarks/results/atcg3d_macos_arm64_2026-07-15.json`. Current
schema-v8 core-memory, evolving-event, output, and application measurements
from 2026-07-25 are recorded first below and in
`benchmarks/results/atcg3d_schema_v8_macos_arm64_2026-07-25.json`.

Sections carrying earlier dates are retained as historical evidence and state
their own schema and limitations. They must not be cited as current
schema-v8 throughput, memory, file-size, or checksum measurements. In
particular, the complete 10^7 VTK/checkpoint round trip remains the 2026-07-15
measurement; the current 2026-07-25 run separately validates 10^7 core/event
paths and a 10^6 full I/O round trip.

## Measured platform

- Mac Studio, Apple M1 Ultra, 20 CPU cores, 128 GB RAM
- macOS 26.5.2 arm64; AppleClang 21.0.0
- Release build; HDF5 2.1.1; VTK 9.6.2
- ParaView 6.1.1 for native-loader, trame-backend, and rendering measurements

Linux and macOS x86_64 were not run. The implementation has no macOS-only
simulation path, but this report does not claim either platform as tested.

## 2026-07-25 schema-v8 and Studio acceptance

The full VTK+HDF5 Release build passed 24/24 CTest tests, and a separate
legacy-only Release build passed 2/2. Coverage includes
schema-v8 stable-slot journal reconstruction, corrupt/wrong-version failure,
continuous-versus-resume checksum, indexed event-heap behavior, one/multiple
thread checksum equality, live overwrite-only VTK-HDF, ParaView series refresh,
viewer debounce, run-control requests, migration, vascular growth, and the
legacy shared density-growth function. All five repository YAML profiles
passed strict `atcg3d --config ... --dry-run`.

The event queue now stores at most one node per actor/event kind and updates
that node in `O(log N)`. The measured runs below reported zero queue rebuilds.
The exact-time proposal cache uses rolling, bounded 8192-event windows and
spatial block versions. On the real ten-million-cell run it used all 18
workers, accepted 7,959 prefetched proposals, invalidated 183, and computed
1,858 at commit time.

| Cells | Threads | Events | Build | Event run | Events/s | Pending indexed events | Peak RSS | Checksum |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10,000,000 | 1 | 10,000 | 16.189 s | 0.987 s | 10,126.84 | 30,000,000 | 6,737,838,080 B | 10612870120539738236 |
| 10,000,000 | 18 | 10,000 | 8.050 s | 0.907 s | 11,022.60 | 30,000,000 | 6,737,756,160 B | 10612870120539738236 |

These are real typed cells, sparse occupancy/density indexes, indexed
migration/division/death schedules, and committed migration events. They are
not metadata-only counts. The equal checksum confirms deterministic output;
the 1.09x event-phase speedup is the measured result, not a claim of linear
18-core scaling.

A separate ten-million-cell core run, without VTK/checkpoint I/O, measured:

| Cells | Logical CellStore width | Tracked core | Core B/cell | Peak RSS | Build | 1M preview sample |
|---:|---:|---:|---:|---:|---:|---:|
| 10,000,000 | 99 B | 1,284,573,312 B | 128.457 | 3,776,200,704 B | 14.482 s | 0.151 s |

The tracked total is a 1,067,108,864 B CellStore, 44,982,824 B sparse grid,
and 172,481,624 B density index. The event run's 6.74 GB peak includes the
30-million-node indexed scheduler and remains below the 16 GB target.

A current one-million-cell compression-level-1 full I/O round trip wrote and
read a 7,411,635 B full VTK-HDF in 0.191/0.041 s and wrote/read/restored a
754,304 B self-contained base checkpoint in 0.516/3.570 s. The data are highly
regular and compress unusually well, so these sizes are software-validation
results, not production file-size forecasts. Current 10^7 full I/O was not
rerun; the dated complete 10^7 table below remains historical.

`ATCG3D Studio.app` built successfully with Tauri 2/Rust 1.97.1, is 14 MB,
contains a native icns, and passed
`codesign --verify --deep --strict`. It is ad-hoc signed, not notarized. The
simulator/trame viewer remain external executables selected in its right-side
panel.

The `10^8` architectural ceiling was not executed. Extrapolating bytes per cell
is insufficient to claim throughput because event count, density work, vessel
topology, and allocator overhead change with biology. The fixed-width slot
space can address the target, but 100-million-cell memory and elapsed time
remain unverified.

## 2026-07-24 scheduler, migration, and core-memory measurements

The current production profile uses
`deterministic_exact_window_v3`: proposals inside the configured look-ahead
window are calculated in parallel, but biological event times are never
rounded or overwritten. Commits retain exact time, fixed event-kind
precedence, stable seeded conflict priority, and UID ordering. Migration
direction/footprint candidates use fixed-capacity storage, local-density moves
update only the symmetric difference of old/new windows, and an active cell
has one activation-end heap event rather than appending another copy after
every migration.

The Release build passed 23/23 CTest tests plus the three-case ParaView
checkpoint-materializer test. An activated-r density-cone event benchmark
advanced 100,000 real cells for 10,000 migration events:

| Threads | Events/s | Run time | Workers observed | Pending events | Peak RSS | Checksum |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 9,007.76 | 1.1102 s | 1 | 300,000 | 125,812,736 B | 501615978897874597 |
| 18 | 13,223.42 | 0.7562 s | 18 | 300,000 | 127,205,376 B | 501615978897874597 |

This is a 1.47x measured speedup for the expensive activated-r proposal path,
with an identical biological checksum. Cheap sparse migration does not benefit
from forcing all 18 workers: the 100,000-cell/100,000-event comparison measured
25,296.61 events/s with one thread and 18,579.11 events/s with 18. The adaptive
policy therefore keeps a separate event-granularity limit; thread use is not
itself treated as a performance success. Raising the local refresh threshold
from 8 to 128 after the incremental-density change improved the activated
18-thread comparison from 8,384.45 to 13,223.42 events/s.

A fresh synthetic core run created 10,000,000 real typed cells, populated the
sparse occupancy and density indexes, and selected a stable one-million-cell
preview:

| Build | Core tracked | CellStore | Grid | Density | Core bytes/cell | Peak RSS | Build | Preview sample |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| current | 1,207,464,448 B | 990,000,000 B | 44,982,824 B | 172,481,624 B | 120.746 | 3,729,670,144 B | 13.659 s | 0.133 s |

The current logical CellStore width is 99 bytes/slot. This run deliberately
skipped VTK-HDF and checkpoint I/O because those full ten-million-cell
measurements remain the older dated evidence below. A one-million-cell
evolving benchmark also advanced 100,000 real migration events, but a complete
ten-million-cell biological trajectory has not been timed.

`10^8` had not been executed. Fixed-width slot/UID types and sparse indexes did
not impose a ten-million-cell addressing limit, but the binary heap used by
that dated build and serial deterministic commit were scale bottlenecks. The
current indexed-heap result is in the 2026-07-25 section; the 100-million
figure remains an architectural capacity target, not a verified throughput or
memory claim.

## 2026-07-20 historical incremental-storage measurements

At that date the writer supported level-1 VTK-HDF/HDF5 compression and hourly
schema-v5 row deltas between schema-v4 bases. Those formats are retained for
reading; new files use v6 bases and v7 field deltas. A 100,000-cell synthetic
comparison on the measured platform gave:

| Level | Full VTK-HDF | Full write | Checkpoint | Checkpoint write |
|---:|---:|---:|---:|---:|
| 0 | 5,563,069 B | 0.00576 s | 11,816,760 B | 0.19715 s |
| 1 | 1,513,532 B | 0.02132 s | 178,992 B | 0.20719 s |

This synthetic state is highly regular, so its checkpoint ratio is not a
production estimate. On the stopped 11,131,595-cell production run, repacking
one real uncompressed full frame with shuffle+gzip-1 reduced 416 MB to 57 MB in
2.26 s; a real one-million-point preview reduced 38 MB to 6.0 MB in 0.21 s.
The new checkpoint writer's v4-base/v5-delta restore was tested over a two-delta
chain against continuous-run checksums, including birth/death, slot reuse,
lineage suffix, missing/corrupt parent rejection, and a production OutputManager
base-then-delta integration path. These measurements do not yet constitute a
new 10^7 evolving-run completion benchmark.

## Build and test matrix

The schema-v3 implementation was rebuilt and tested on 2026-07-17:

| Build | Legacy 2D | 3D | HDF5 | VTK-HDF | Result |
|---|:---:|:---:|:---:|:---:|---:|
| legacy Release | yes | no | no | no | 2/2 passed |
| core Release | no | yes | no | no | 21/21 passed |
| HDF5 Release | no | yes | yes | no | 22/22 passed |
| VTK + HDF5 Release | no | yes | yes | yes | 23/23 passed |
| `-Wall -Wextra -Wpedantic -Werror` Release core | no | yes | no | no | 21/21 passed |

A targeted ASan + UBSan run of the new lesion-index test also passed. The
complete schema-v3 sanitizer suite was not rerun, so the older full-suite
sanitizer result is not promoted as current evidence.

Apple AddressSanitizer does not support LeakSanitizer, so leak detection was
disabled; AddressSanitizer and UndefinedBehaviorSanitizer remained enabled.
The full-feature build is reproducible with:

```sh
cmake -S . -B build-3d -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DATCG_BUILD_LEGACY_2D=OFF \
  -DATCG_BUILD_3D=ON \
  -DATCG3D_ENABLE_HDF5_CHECKPOINT=ON \
  -DATCG3D_ENABLE_VTKHDF=ON
cmake --build build-3d --parallel
ctest --test-dir build-3d --output-on-failure
```

`atcg3d --config ... --dry-run` accepted both schema-v3 YAML profiles. An
attempted `--set simulation.threads=4` override failed nonzero, confirming that
model and run values are supplied only by YAML.

## Current schema-v3 ten-million-cell lesion/angiogenesis progression

The per-lesion index, independent Poisson scheduler, bidirectional vascular
growth, displacement, and density-relief paths were measured after the final
schema-v3 changes on 2026-07-18, after enforcing an outward vessel speed of
2.0 voxel/hour above the activated-r maximum path speed. The machine-readable
result is
`benchmarks/results/atcg3d_angiogenesis_schema_v3_macos_arm64_2026-07-18.json`:

```sh
build-full/atcg3d_angiogenesis_scale_benchmark \
  --cells 10000000 --events 4096 --threads 18 --profile production
```

The run reached 720 model hours in 1,713 events. It attempted and committed
eight roots, committed 1,419 vascular growth steps, displaced 6,364 cells, and
ended with 9,993,636 live cells. The result contained 1,427 vessel nodes,
14,969 occupied vessel voxels, and 1,038,282 vascular-influence voxels. Build
and run time were 3.100 s and 21.164 s respectively (80.94 events/s), with
checksum `2698773685244278503`.

Peak RSS was 3,492,675,584 B. Explicitly tracked components used
1,156,777,916 B, or 115.68 B per initial cell; this includes a 980,032,768 B
cell store and a 5,420,380 B lesion index. The benchmark freezes the synthetic
background cells' migration, division, and death schedules so that vascular
scaling is interpretable. It is evidence for the ten-million-cell lesion and
vascular paths, not a wall-time estimate for a fully evolving tumour.

An initial run with `--events 256` stopped normally at its configured event
limit (208.53 model hours) and was rejected by the benchmark because it had not
reached the required 720-hour horizon. It is not counted as a passing result.

## Historical synthetic 10^3 through 10^7 storage benchmark

The four sizes were run in the required order with:

```sh
build-3d/atcg3d_scale_benchmark \
  --cells N --require-vtk --require-checkpoint \
  --directory RESULTS --result RESULTS/result.json
```

Each run creates N real typed `CellStore3D` records, inserts every anchor into
the sparse chunk grid and density index, samples a deterministic preview,
writes and reads a full VTK-HDF frame, writes the preview, writes an HDF5
checkpoint, restores a second simulation, and compares its checksum. It does
not modify a count-only metadata field.

| Cells | Core B/cell | CellStore | Sparse grid | Density index | Peak RSS | Build/index | Preview sample |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 239.448 | 98,000 B | 131,152 B | 10,296 B | 22.04 MB | 0.00045 s | 0.00004 s |
| 100,000 | 118.279 | 9.80 MB | 1.05 MB | 0.98 MB | 101.50 MB | 0.01705 s | 0.00016 s |
| 1,000,000 | 115.666 | 98.00 MB | 8.39 MB | 9.27 MB | 701.68 MB | 0.20598 s | 0.00161 s |
| 10,000,000 | 111.749 | 980.00 MB | 44.98 MB | 92.50 MB | 6.857 GB | 2.07838 s | 0.12779 s |

`CellStore3D` is 98 logical bytes per allocated slot, inside the 80–100 byte
target. At ten million cells the directly accounted steady simulation core is
1,117,486,844 bytes. The 6.857 GB process peak occurs during the deliberately
strict round trip: the source simulation remains alive while HDF5 read buffers
and a second restored simulation coexist. It is still below the 16 GB target.

| Cells | Full write | Full read | Full size | Preview size | Checkpoint write | Read + restore | Checkpoint size |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.00203 s | 0.00119 s | 1.83 MB | 1.83 MB | 0.00389 s | 0.00203 s | 0.16 MB |
| 100,000 | 0.00365 s | 0.00154 s | 4.16 MB | 4.16 MB | 0.10231 s | 0.11437 s | 11.74 MB |
| 1,000,000 | 0.02211 s | 0.00636 s | 32.06 MB | 32.06 MB | 1.05785 s | 1.16944 s | 117.04 MB |
| 10,000,000 | 0.23794 s | 0.04399 s | 311.16 MB | 32.06 MB | 10.71215 s | 11.88895 s | 1.170 GB |

The 10^7 full frame contains every cell and required fixed-width array. Its
311,162,142 bytes are inside the requested approximate 260–320 MB range. The
preview contains the configured maximum one million stable-hash samples.

## Historical legacy-mapped initialization and biology

```sh
build-3d/atcg3d_initialization_benchmark \
  --config configs/atcg3d_legacy_2d_mapped_v3.yaml
```

The production profile created 568,900 non-overlapping cells in 0.90792 s:
284,450 r and 284,450 K; 47,964 stage 0 and 520,936 stage 1. Biological volume
was 904,648 voxel^3, with 64 cell chunks and 67,704 exposed faces. Accounting
was 55,752,200 B for cells, 8,393,480 B for the grid, and 5,963,236 B for
density; peak RSS was 199,507,968 B. Checksum:
`7363027809122967352`.

The YAML-only smoke profile was also advanced for 24 model hours through the
real executable. It progressed 377 events, attempted 336 migrations, committed
331 migrations and 41 divisions, recorded no death or angiogenesis event, and
grew from 64 to 105 live cells. Checksum: `7994067227480867429`.

## Historical real migration-event progression

```sh
build-3d/atcg3d_event_benchmark --cells N --events 10000 --threads T
```

This benchmark creates typed cells, occupancy, density, and event queues, then
commits actual migration events. Division work is intentionally placed in the
future so migration throughput remains interpretable.

| Initial cells | Events | Threads | Run time | Events/s | Migration commits | Simulated hour | Peak RSS | Checksum |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 10,000 | 1 | 0.17958 s | 55,685 | 9,997 | 9.999673 | 17,727,488 B | 17972680546450520886 |
| 1,000 | 10,000 | 4 | 0.40420 s | 24,740 | 9,997 | 9.999673 | 18,071,552 B | 17972680546450520886 |
| 100,000 | 10,000 | 1 | 0.52622 s | 19,004 | 10,000 | 0.108109 | 82,542,592 B | 13179999894853332298 |
| 100,000 | 10,000 | 4 | 0.74767 s | 13,375 | 10,000 | 0.108109 | 82,935,808 B | 13179999894853332298 |
| 1,000,000 | 10,000 | 4 | 0.67774 s | 14,755 | 10,000 | 0.019845 | 595,492,864 B | 989825554614778225 |
| 10,000,000 | 10,000 | 4 | 0.76781 s | 13,024 | 10,000 | 0.010991 | 5,815,697,408 B | 459412592071965108 |

One and four threads produce identical state and time at both paired sizes.
These small same-time proposal batches do not amortize OpenMP overhead, so four
threads are slower here; determinism does not imply speedup.

At one migration event per cell per hour, holding population and event mix
fixed would produce 7.2 billion events over 720 hours. At the measured
10-million-cell rate that is about 6.40 days. This is a synthetic extrapolation,
not a production promise: the production model has heterogeneous rates,
division, death, density activation, vessels, and changing population. A real
ETA must use the run's observed event count and throughput.

## Historical ten-million-cell angiogenesis progression

```sh
build-3d/atcg3d_angiogenesis_scale_benchmark \
  --cells N --events 256 --threads T
```

This benchmark uses a dense true-3D tumour and a deliberately high synthetic
seed intensity so a Poisson arrival is exercised during a short run. It is a
software scale test, not the production scientific rate. Specifically, it uses
one root, two active tips, diameter 1, and inward/outward maximum length 12.
The production legacy-mapped YAML instead defaults to up to 64 roots, diameter
3, and length 128; those scientific settings were not substituted into this
short progression benchmark.

| Cells | Threads | Final cells | Events / model time | Run time | Events/s | Growth commits | In/out nodes | Deaths / displacements | Influence voxels | Peak RSS | Checksum |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 1 | 919 | 93 / 24 h | 0.02815 s | 3,304 | 14 | 3 / 11 | 78 / 3 | 14,564 | 17,137,664 B | 17584560740882085401 |
| 1,000 | 4 | 919 | 93 / 24 h | 0.02829 s | 3,287 | 14 | 3 / 11 | 78 / 3 | 14,564 | 17,203,200 B | 17584560740882085401 |
| 10,000,000 | 1 | 9,999,764 | 248 / 4.121838 h | 0.17527 s | 1,415 | 19 | 8 / 11 | 228 / 8 | 16,980 | 3,421,208,576 B | 11201788409937048902 |
| 10,000,000 | 4 | 9,999,764 | 248 / 4.121838 h | 0.17697 s | 1,401 | 19 | 8 / 11 | 228 / 8 | 16,980 | 3,423,322,112 B | 11201788409937048902 |

Every row produced exactly one attempted/committed root with no seed rejection;
both tips reached terminal state. At 10^7 cells, tracked
cell/grid/density/vessel components were 1,118,460,540 B. The tracked subtotal
intentionally excludes the surface hash, event heap, and allocator overhead,
which are included in peak RSS.

The default CTest also runs a 1,000-cell angiogenesis smoke in both thread
modes. It checks inward and outward nodes separately, permanent vessel
occupancy, displacement, influence-field activation, and equal one-/four-thread
checksums. The 256-event rows above share checksum `17584560740882085401`.
Three ten-million tests are opt-in with
`ATCG3D_ENABLE_LARGE_TESTS=ON` and `RUN_SERIAL`: the core storage round trip,
the synthetic angiogenesis progression, and the production-parameter
angiogenesis progression.

### Production-parameter angiogenesis progression

```sh
build-3d/atcg3d_angiogenesis_scale_benchmark \
  --cells 10000000 --events 100000 --threads T --profile production
```

The historical production-like integration benchmark did not parse a YAML
file; it programmatically constructed the then-current default configuration
and advanced ten million cells for 720 model hours: 10 sites per 30 days;
activation/deactivation
volumes 100,000/80,000 voxel^3; at most 64 roots and 128 active tips; diameter
3; inward/outward speeds 0.5/2.0 voxel/h; inward minimum path budget 128
voxels with lesion-scale expansion (tortuosity 1.5, exit margin 16, hard cap
4096); outward maximum length 128;
external-connection distance 64; influence cutoff 12; and maximum density
relief 0.5.

For benchmark isolation, only the synthetic `CellInit` background schedules
were neutralized. Vessel occupancy, inward cell
replacement/displacement, and vascular density-influence updates all execute
through the real model paths.

| Threads | Final cells | Events / model time | Seeds / roots / rejected | Growth | Root / in / out nodes | Tips / active | Displacements | Occupied / influenced voxels | Run time | Events/s | Peak RSS | Checksum |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 9,993,457 | 1,090 / 720 h | 7 / 7 / 0 | 1,083 | 7 / 584 / 499 | 14 / 5 | 6,543 | 11,506 / 640,700 | 3.70090 s | 294.52 | 3,450,699,776 B | 8388156886843905651 |
| 4 | 9,993,457 | 1,090 / 720 h | 7 / 7 / 0 | 1,083 | 7 / 584 / 499 | 14 / 5 | 6,543 | 11,506 / 640,700 | 3.70932 s | 293.85 | 3,450,634,240 B | 8388156886843905651 |

The tracked production subtotal is 1,143,739,664 B: 980,032,768 B cell
store, 44,982,824 B cell grid, 92,481,624 B density index, 22,122,776 B vessel
grid, 4,005,336 B influence field, and 114,336 B combined vessel node/tip
stores. The equal checksum confirms deterministic one- and four-thread results.

## Historical checkpoint, output, and visualization measurements

### July 2026 production-checkpoint hot-path measurement

A real 631 h checkpoint containing 2,771,229 live cells was resumed with
output disabled and advanced to 631.2 h on the local 18-worker macOS build.
The comparison used identical YAML dynamics and ended with the same
`completed_events=72628146`, population counters, and checksum
`8641331037086676441`.

| Implementation | Restore + 0.2 h wall time | User CPU time |
|---|---:|---:|
| block-estimated 6³ growth density | 46.67 s | 181.09 s |
| exact incremental per-slot 6³ counts | 36.50 s | 55.60 s |

The incremental build's restore-only measurement was 14.52 s, so its measured
continuous 0.2 h portion was about 21.98 s. The result is checkpoint- and
population-specific; it is not a completion-time promise for the remaining
2160 h run. A 0.005 h speculative proposal window and a disabled proposal
window both took about 46.5 s before the 6³ optimization, so proposal-window
tuning was not credited as a speedup. Raising
`min_refresh_items_per_thread` from 8 to 128 increased the same pre-optimization
measurement from 46.67 s to 49.11 s and was rejected.

The asynchronous output path was also changed to write directly from its
immutable cell/lineage/vascular snapshot. It no longer rebuilds a second cell
grid, density index, lesion index, and event queue in the writer thread. The
default async VTK-HDF/checkpoint tests cover round-trip behavior; production
frame sizes and I/O time still depend on the live cell count and filesystem.

An output-enabled run was checkpointed at 2 h, resumed to 4 h, and compared
with an uninterrupted 4 h run. Both final biological checksums were
`4062273419464254262`. Resume produced five preview frames, three full frames,
and five vessel frames, with every full frame paired to a same-time preview.
No `.tmp` or PNG file was present. The manifest used schema version 3 and the
preview carried `total_cell_count` as VTK-HDF FieldData. Continuous and resumed
frame timestamps may differ at the 2 h checkpoint boundary because of output
scheduling, while the final biological state remains checksum-identical.
Native ParaView loaded both preview and full series, and the trame app was
constructed against the live run directory with real ParaView Python.

For the final ten-million synthetic frame, ParaView's server-side off-screen
1280x720 Point Gaussian benchmark measured:

| Dataset | Points | Data load | First render | Camera-rotation FPS | Peak RSS |
|---|---:|---:|---:|---:|---:|
| preview | 1,000,000 | 0.02955 s | 0.13755 s | 525.52 | 558,858,240 B |
| full | 10,000,000 | 0.11191 s | 0.42645 s | 72.69 | 1,852,456,960 B |

The preview and even the full frame exceed 30 FPS in this local server-side
test, and full loading is inside the 1–5 s idle target. Rendering used Point
Gaussian with display radius applied through `Calculator + ScaleByArray`. The
preview frame also
reported its true `total_cell_count=10,000,000` while displaying one million
sampled points. This does not measure remote-browser latency. The
dependency-free debounce test confirms that 100 rapid slider inputs request
only the final exact-time full frame and discard stale completions.

The packaged ParaView emitted a warning about its optional OpenVKL CPU module;
VTK-HDF point/Line loading, Tube display, native loading, trame construction,
and rendering all completed successfully.
