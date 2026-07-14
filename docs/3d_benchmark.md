# ATCG3D scale benchmark

## Commands

Build used for the measurements below:

```sh
cmake -S . -B /private/tmp/ver7-3d-hdf5 \
  -DATCG_BUILD_LEGACY_2D=OFF \
  -DATCG3D_ENABLE_HDF5_CHECKPOINT=ON \
  -DATCG3D_ENABLE_VTKHDF=OFF
cmake --build /private/tmp/ver7-3d-hdf5 -j4
```

Synthetic scale, repeated in order for 10³, 10⁵, 10⁶, and 10⁷:

```sh
/private/tmp/ver7-3d-hdf5/atcg3d_scale_benchmark \
  --cells N --skip-vtk --require-checkpoint \
  --directory /private/tmp/atcg3d-scale-N \
  --result /private/tmp/atcg3d-scale-N/result.json
```

The benchmark creates N real typed records, inserts every cell into the sparse
grid and density index, hashes a deterministic preview of at most one million
UIDs, writes an HDF5 checkpoint, reads it, reconstructs another simulation, and
checks the state checksum. It does not change a count-only metadata field.

## Measured synthetic results

Platform: Mac Studio, Apple M1 Ultra (20 cores), 128 GB RAM, macOS arm64;
Apple clang 21.0.0; HDF5 1.14.6. The build did not set a CMake optimization
configuration explicitly.

| Cells | Core bytes/cell | Peak RSS | Build index | Preview hash | Checkpoint write | Checkpoint read+restore | Checkpoint size |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 227.36 | 7.00 MB | 0.00022 s | 0.00002 s | 0.00118 s | 0.00056 s | 106,464 B |
| 100,000 | 106.10 | 49.61 MB | 0.00777 s | 0.00016 s | 0.01079 s | 0.01732 s | 8,917,464 B |
| 1,000,000 | 103.66 | 430.67 MB | 0.07738 s | 0.00157 s | 0.11938 s | 0.17207 s | 89,017,464 B |
| 10,000,000 | 99.74 | 4.36 GB | 0.83595 s | 0.13051 s | 1.05891 s | 1.77238 s | 890,017,464 B |

The 10-million peak includes temporary row records, HDF5 column buffers, the
original simulation, checkpoint read buffers, and a simultaneously restored
second simulation. The steady core accounting is 900,000,000 bytes for
CellStore, 44,982,824 bytes for sparse occupancy, and 52,417,624 bytes for the
density index: 997,400,448 bytes total, below the 16 GB target. The density
index total includes compact boundary-anchor entries for exact partial-block
queries.

## Real event progression

```sh
/private/tmp/ver7-3d-hdf5/atcg3d_event_benchmark \
  --cells N --events 100000 --threads T
```

| Initial cells | Events | Threads | Wall time | Events/s | Simulated hours | Checksum agreement |
|---:|---:|---:|---:|---:|---:|:---:|
| 1,000 | 100,000 | 1 | 2.966 s | 33,719 | 400 | yes |
| 1,000 | 100,000 | 4 | 2.717 s | 36,801 | 400 | yes |
| 10,000 | 100,000 | 1 | 3.382 s | 29,571 | 40 | yes |
| 10,000 | 100,000 | 4 | 2.906 s | 34,415 | 40 | yes |
| 100,000 | 100,000 | 1 | 7.795 s | 12,828 | 4 | yes |
| 100,000 | 100,000 | 4 | 4.884 s | 20,475 | 4 | yes |

These runs advance actual scheduled migration/division/death events. They are
not a 10-million-cell biology throughput claim. Time to a requested biological
horizon must be estimated as `actual events / measured events_per_second` for
that population and parameter profile. The number of events per simulated hour
changes with population and rates, so no honest 10-million-cell completion-time
promise is made without that benchmark.

## Not measured in this environment

VTK and trame were not installed. Consequently full VTK-HDF write/read size,
preview rendering FPS, browser interaction, and idle 10-million-point load time
were not measured. The writer uses the required 31 bytes/cell of raw point and
array payload, implying about 310 MB before HDF5 overhead/compression, but this
is a schema calculation rather than a measured file result. Build with
`ATCG3D_ENABLE_VTKHDF=ON` and run:

```sh
atcg3d_scale_benchmark --cells 10000000 \
  --require-vtk --require-checkpoint --directory RESULTS --result RESULTS/result.json
```

Do not record the VTK target as passed until that command and the conditional
VTK reader round-trip test have run. Linux and macOS x86_64 are also not yet
tested; the implementation uses standard C++20/OpenMP/HDF5/VTK APIs and has no
macOS-only simulation code.
