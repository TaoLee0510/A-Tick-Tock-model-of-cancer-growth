# Hosted verification

The simulator workflow builds the HDF5-enabled configuration on Linux and
macOS, runs every unit/regression test, then runs the labeled ensemble tests.
Reports and CTest logs are uploaded even when a check fails.

The published legacy PDE checksum fixture was recorded with Apple clang and
the Darwin math library. It remains the macOS regression oracle. Bitwise
floating-point checksums are not portable across math-library/compiler
implementations: the initial hosted Linux run built successfully but failed
this fixture comparison.

Linux therefore also builds the independent, pinned pre-alignment source
`e9c37820e113e2703407871459870b16ff74936b` with its current compiler and
libraries. The same regression driver dumps its six legacy structured and
continuum checksums. `ATCG_LEGACY_PDE_FIXTURE` selects that file for the current
implementation's strict bitwise comparisons. Expected results are generated
from the pinned implementation, never from the current model. The original
fixture and all published model arithmetic remain unchanged. A change to the
pinned reference requires explicit review.

The fixture override is a test-only facility. A relative override is resolved
against the source root, so CTest's working directory does not affect it.
Native restart and thread-count checks still run on both platforms.
