# Structured PDE operator composition

`StructuredOperatorStrategies3D` selects six immutable-in-use policy objects
once in the model constructor. The objects describe arithmetic contracts;
native checkpoint versions are selected separately in the same table. Runtime
stepping invokes activation, transport, exchange, reaction, nutrient and
vascular objects in the published sequence. The numerical kernels live in
six operator source files and share the model's fields, bounds and cohort
stores. No operator owns a second copy of biological state.

| Operator | Schema-selected contracts |
| --- | --- |
| Activation | Mean remaining clock; v5 grid-local refractory hysteresis; v7 transported refractory cohorts; explicit duration/rate distributions retain their selected grids. |
| Transport | v1 fixed jumps; v2 resource guidance; v3 directional sectors; v4 cached cohort substeps; v5 nutrient/footprint guidance; explicit normal transport strings select axial, fixed or feasible lattice jumps. |
| Exchange | Pre-v3 disabled closure; v3 direction-dependent conservative r/K flux with the selected stage policy. |
| Reaction | v3 local-density conversion; v4 post-transport conversion; v7 exact ABM growth edge; explicit true normal expectation and feasible small-daughter birth strings. |
| Nutrient | Legacy refreshed quasi-steady field; v5 transient shared resources; v6 moving front; v7 footprint-averaged consumers. |
| Vascular | Published v1 density field or explicit shared VEGF process; v7 shared geometry; v8 fraction-based exclusion; v9 empty-cell traversal optimization; v13 removal diagnostics. |

Kernel expressions, parenthesization, loop order, rounding, cutoff comparisons
and the macro-step sequence remain unchanged during the split. Schema
comparisons in `structured_pde_model.cpp` and the six kernels are replaced by
the constructor-selected capabilities. Model strings remain explicit choices:
upgrading a schema does not select a new migration, reaction or vascular law.
The core model retains initialization, bounds, diagnostics, native hashing and
checkpoint serialization. Derived policies are reconstructed from the saved
configuration and are not written as new checkpoint fields.

Verification compares native JSON fixtures from before the split, including
published trajectories, v14 statistical traces, v15 vascular realizations and
v4 hybrid invasion realizations. Default CTests also cover conservation,
thread counts, native checkpoint continuation and each historical schema.
Any arithmetic difference must be fixed or introduced through a new contract;
a maintenance split is not permission to change an existing result.
