# ATCG3D validation and scale benchmark

Last measured: 2026-07-15. The machine-readable report is
`benchmarks/results/atcg3d_macos_arm64_2026-07-15.json`.

## Measured platform

- Mac Studio, Apple M1 Ultra, 20 CPU cores, 128 GB RAM
- macOS 26.5.2 arm64; AppleClang 21.0.0
- Release build; HDF5 2.1.1; VTK 9.6.2
- ParaView 6.1.1 for native-loader, trame-backend, and rendering measurements

Linux and macOS x86_64 were not run. The implementation has no macOS-only
simulation path, but this report does not claim either platform as tested.

## Build and test matrix

Six clean builds were configured, built, and tested:

| Build | Legacy 2D | 3D | HDF5 | VTK-HDF | Result |
|---|:---:|:---:|:---:|:---:|---:|
| legacy Release | yes | no | no | no | 2/2 passed |
| core Release | no | yes | no | no | 19/19 passed |
| HDF5 Release | no | yes | yes | no | 20/20 passed |
| VTK + HDF5 Release | no | yes | yes | yes | 21/21 passed |
| ASan + UBSan Debug core | no | yes | no | no | 19/19 passed |
| `-Wall -Wextra -Wpedantic -Werror` Release core | no | yes | no | no | 19/19 passed |

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

`atcg3d --config ... --dry-run` accepted both schema-v2 YAML profiles. An
attempted `--set simulation.threads=4` override failed nonzero, confirming that
model and run values are supplied only by YAML.

## Synthetic 10^3 through 10^7 storage benchmark

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

## Legacy-mapped initialization and biology

```sh
build-3d/atcg3d_initialization_benchmark \
  --config configs/atcg3d_legacy_2d_mapped_v2.yaml
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

## Real migration-event progression

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

## Ten-million-cell angiogenesis progression

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

The production integration benchmark advances ten million cells for 720 model
hours with the YAML defaults: 10 sites per 30 days; activation/deactivation
volumes 100,000/80,000 voxel^3; at most 64 roots and 128 active tips; diameter
3; inward/outward speeds 0.5/0.25 voxel/h; inward/outward maximum lengths 128;
external-connection distance 64; influence cutoff 12; and maximum density
relief 0.5.

For benchmark isolation, only the synthetic `CellInit` background schedules
are neutralized. The benchmark does not alter the YAML, legacy-mapped density
or death rules, or any vessel parameter. Vessel occupancy, inward cell
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

## Checkpoint, output, and visualization

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
