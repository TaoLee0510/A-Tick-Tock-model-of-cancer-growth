# CellStore Migration Plan

This document tracks the staged migration from the legacy `cell_array(row, col)`
matrix layout to a column-oriented `CellStore`.

## Current state

- Column semantics are centralized in `ATCG/cell_columns.hpp`.
- `ATCG/cell_store.hpp` provides a column-oriented container and conversion
  helpers for the current Blitz `Array<double, 2>` boundary.
- `deltah_recalculation` has a `CellStore` overload and the low-density initial
  growth path now reads the migration-rate column through `CellStore`.

The matrix layout is still the owner for most simulation state. `CellStore` is
currently a migration layer, not the single source of truth.

## Migration order

1. Keep `cell_array` at input/output and visualization boundaries while moving
   hot loop internals one function at a time.
2. Convert read-heavy calculations before mutation-heavy lifecycle code:
   `density_calculation`, then `density_growth_rate_calculation_1`, then
   `death_judgement`.
3. Move migration kernels next: `random_migration` and `migration`. Preserve the
   existing stateless RNG event indexing so every cell migration still consumes
   its own random number.
4. Move division kernels after migration: `division`, `free_living_division`,
   and `free_living_division_single_thread`.
5. Replace row sorting last. `sortRow` and `sortRowReverse` need a shared
   permutation routine that reorders every `CellStore` column consistently.
6. After the hot kernels operate on `CellStore`, remove repeated conversion in
   the growth drivers and keep conversion only at save/recovery boundaries.

## Guardrails

- One kernel family per commit, with compile verification after each code step.
- Do not mix persistent RNG state into `CellStore`; derive random values from
  `(cell id, time step, event id)` as the current stateless path does.
- Treat append, resize, erase, and sort as whole-store operations. A row index is
  only valid if every column has the same logical row count.
- Avoid converting the whole matrix when a function only needs a small set of
  columns. Use selected-column conversion during the transition period.
- Keep Blitz-compatible overloads until all call sites of a kernel family are
  migrated.

## Main risks

- Mutation-heavy functions update many columns and `Visual_range` together. They
  should be migrated only after their read/write column set is explicit.
- Sorting by a single column must reorder all columns with the same permutation;
  otherwise cells will silently mix state from different rows.
- Repeated matrix-to-store conversion inside the time loop can cancel the cache
  benefit. Full benefit comes when growth drivers own a long-lived `CellStore`.
