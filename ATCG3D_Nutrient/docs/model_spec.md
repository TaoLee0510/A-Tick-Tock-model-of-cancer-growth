# ATCG3D effective-nutrient model specification

## Scope

The nutrient model is an extension of the event-driven spatial PDMP in
`ATCG3D`; it is not a continuum replacement for individual cells. Every cell
retains its UID, type, stage, anchor, footprint, migration state, division work,
death deadline, event sequence, and conflict behavior. The vascular network
retains its graph topology and rasterized capsule occupancy.

The only new continuum state is an effective nutrient availability field
`N(x,t)`. It represents the aggregate growth-supporting effect of oxygen,
glucose, and other unresolved diffusible resources. It must not be interpreted
as a calibrated physical oxygen concentration.

## Nutrient equation

The dynamic form is

```text
dN/dt = div(D_N grad(N))
        - lambda_N N
        - Q_cells(x, N, cell_state)
        + kappa_v I_vessel(x) (N_v - N).
```

Cell consumption is assembled from discrete agents:

```text
Q_cells(x, N, cell_state)
  = sum_i q_max(type_i, stage_i, state_i)
          * N / (K_q(type_i) + N)
          * K_epsilon(x - anchor_i).
```

`I_vessel` is evaluated from the existing rasterized perfused-vessel capsule,
not from a smoothed replacement of the vascular graph. `K_epsilon` distributes
one cell's sink over its biological footprint or a separately specified local
kernel.

The v1 solver uses the quasi-steady form obtained by setting `dN/dt` to zero.
`solver: deterministic_quasi_steady_jacobi_v1` selects fixed-count relaxed
Jacobi iterations. A future dynamic solver must use a new strategy name and
checkpoint its time-integrator history.

## Backward-compatible capacity coupling

The reference vascular rule uses

```text
effective_count = (1 - relief) * raw_count.
```

Define a nutrient-dependent carrying-capacity multiplier

```text
M(N) = 1 + (M_max - 1) * S(N),
```

where `S(N)` is bounded in `[0,1]`. The nutrient model supplies the density
growth rule with

```text
effective_r_count = raw_r_count / M(N),
effective_K_count = raw_K_count / M(N).
```

Equivalently, both density limits and carrying capacities are multiplied by
`M(N)`. The current maximum relief of `0.5` maps to

```text
M_max = 1 / (1 - 0.5) = 2.
```

Compatibility mode constrains `M(N) >= 1`: the new field represents vascular
support above the unresolved avascular baseline. A later calibrated mode may
allow nutrient deprivation below baseline to slow division work directly or
start a starvation/death process. That later mode is a biological change and
must not be introduced implicitly.

## Hybrid event coupling

A deterministic `environment_refresh` event is added to the scheduler. At each
refresh, cell sinks and perfused-vessel sources are rebuilt in stable order,
the sparse field is solved, and every living cell is refreshed in UID order.
This full refresh is the v1 correctness contract; a future dirty-region
optimization must preserve the same trajectory.

For every living cell the refresh must:

1. integrate division work using the previously stored growth rate up to the
   exact refresh time;
2. evaluate the new local nutrient value and density-growth rate;
3. reschedule division or death without redrawing the division cycle;
4. preserve migration rules unless the configured nutrient scope explicitly
   includes migration.

The v1 model uses a fixed configured refresh interval. Adaptive
refreshing may be added only with deterministic thresholding and checkpointed
state.

## Sparse spatial domain

A dense field over the full expandable ATCG3D coordinate range is forbidden.
The nutrient field must allocate sparse blocks around cells and vessels plus a
configured diffusion halo. Boundary conditions and active-domain expansion
must be deterministic. The existing vascular influence cutoff of 12 voxels can
serve as an initial spatial-range target, but it is not a calibrated nutrient
diffusion coefficient.

## Determinism and restart state

The field solver must use deterministic assembly, traversal, convergence, and
floating-point reduction rules. A checkpoint must preserve, or canonically
reconstruct without ambiguity:

- the nutrient field and active sparse blocks;
- the last and next nutrient refresh times and generation;
- numerical method, boundary-condition, and solver parameters;
- any dynamic-field history required by the chosen time integrator.

In v1, cell sinks, vessel sources, and the next active-block set are transient:
they are rebuilt in stable order from the exact cell and vascular checkpoint at
the next refresh. Field values and the refresh schedule are stored in the
nutrient sidecar.

Adding the nutrient field creates a new model trajectory. It is not expected to
match an old ATCG3D run bit for bit when consumption or PDE-based supply is
enabled. Exact restart and fixed-seed reproducibility are required within the
new model.

## Validation requirements

The registered `atcg3d_nutrient_test` covers:

- the no-source compatibility limit;
- a monotone one-vessel nutrient profile;
- one-cell consumption relative to the no-cell solution;
- capacity multipliers bounded by `M_max`;
- deterministic nutrient refresh events in the common scheduler;
- binary nutrient-state round trips;
- checkpoint/resume equality across a nutrient refresh boundary;
- strict wrapper configuration loading.

The common ATCG3D regression suite is also run with the environment hook
compiled but absent, which protects baseline behavior, thread determinism, and
existing checkpoint schemas.
