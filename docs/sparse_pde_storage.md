# Sparse structured PDE storage

Structured schema v9 adds `storage.model: sparse_zero_pages_v1`. Published
v1-v8 schemas continue using `dense_v1`. The new backend reserves contiguous
virtual address ranges for population, direction/clock and work arrays through
POSIX anonymous mappings. Unvisited zero pages have no private physical
storage. Full zero resets remap the range; they do not depend on platform
`madvise` discard semantics. Native row arithmetic, transport, reaction and
checkpoint order are retained. Dense and sparse v9 fields compare bitwise.

```yaml
storage:
  model: sparse_zero_pages_v1
  maximum_active_voxels: 1000000
```

The v1 sparse backend supports thin-layer runs with static vasculature. Dynamic
VEGF/tip/vessel models and 3D runs remain available in dense storage; their
large-grid memory is not covered by this backend. Nutrient buffers and static
vessel fraction are dense, retaining their full-domain state and source rules.
The active bounding-region budget fails explicitly before an oversized step;
it is not a promise that arbitrary dense or fragmented tumour support fits.
Virtual address reservation is larger than resident memory and is not itself
the memory measure used below.

`ATCG3D_StructuredPDE/config/structured_sparse_v9.yaml` supplies a 10000-square,
48-hour configuration. Turn periodic output on only when needed: complete CSV
and checkpoint fields can be large even when population storage is sparse.
The dedicated verification loads this supplied configuration, runs 48 hours,
computes final diagnostics and checksum, and then checks process peak resident
memory against 8,000,000,000 bytes. On the macOS Apple-clang/libomp reference
build its peak was 4,573,052,928 bytes (4.57 GB), with final mass 668.387 and
checksum 12240940797807697984. Reporting retains full-domain nutrient/vessel
sums while avoiding empty population-page reads. This measures a sparse lesion
and is not a full-domain occupancy benchmark.

```sh
build-codex/atcg3d_sparse_pde_test --large
```

To register that benchmark in CTest, configure with
`-DATCG3D_ENABLE_LARGE_TESTS=ON`; `atcg3d_sparse_pde_10000_memory` is labelled
benchmark. Default CI uses the small equivalence/restart/budget test and the
labeled ensemble validations. Large benchmark output records measured resident
bytes so a regression cannot pass solely from a theoretical byte estimate.
