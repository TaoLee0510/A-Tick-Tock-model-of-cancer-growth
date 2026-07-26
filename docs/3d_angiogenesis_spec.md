# ATCG3D dynamic angiogenesis specification

## Spatial lesions, eligibility, and seed processes

Angiogenesis is a continuing stochastic process, not a one-shot action. The
model does not assume that the tumour nearest the coordinate origin is the only
source. It identifies every spatially disconnected solid lesion, including a
distant metastasis, and evaluates it independently.

The production detector aggregates cells and occupied sites into sparse 8³
coarse blocks. A block becomes core only after crossing the configurable
occupied-fraction, minimum-cell, and optional biological-volume thresholds; a
lower occupied-fraction threshold controls deactivation. Core blocks use
configurable 6- or 26-neighbour connectivity. Non-core occupied blocks within
the configured halo may be attributed to a nearby lesion for statistics and
surface lookup, but never connect two components. A sparse migration trail
therefore cannot falsely merge a primary and a metastatic lesion.

Stable lesion IDs follow the component with greatest prior core-block overlap.
On a split, the retained child keeps the process and new children start
independent eligibility histories. A refresh can contain a split and a merge
simultaneously, so a predecessor that remains as a retained child keeps its
process even when it also overlaps another result. Every removed predecessor
is assigned at most one successor, chosen deterministically, and its process is
transferred exactly once. On a merge, predecessor counters are combined under
the retained ID, the earliest valid pending arrival is kept, and predecessor
queue entries are invalidated. Before combining, every active predecessor
folds `refresh_time - eligibility_started_hours` into its accumulated eligible
time; the active merged process starts a new interval at the refresh time.
Thus different predecessor start times are summed without losing or
double-counting elapsed eligibility.

Existing vessel `source_lesion_id` values remain immutable historical
provenance. A sparse ownership table maps each non-current historical source to
its current lesion, or to zero after removal. Later topology changes retarget
all existing mappings directly, without alias chains. Global tip limits count
all active tips, while a per-lesion limit resolves historical tip sources
through this ownership table. A lesion therefore cannot evade its active-tip
cap merely by merging after vessel birth.

Lesions below both geometric and biological activation thresholds do not
allocate empty Poisson scheduler objects. A disappeared process with no
eligibility or event history is also removed. Processes carrying attempted
arrivals, roots, or accumulated eligibility remain checkpointed so pruning
cannot alter counters or RNG history.

Biological lesion volume is the sum of its uniquely counted cells' configured
stage volumes. Vessel voxels are excluded and co-located stage-2 cells count
separately. A lesion must have the configured minimum number of core blocks and
reach its activation volume before its process becomes eligible; it becomes
ineligible below the lower deactivation volume. Existing vessels continue to
grow after their source lesion loses seed eligibility.

Each eligible lesion has its own density-modulated Poisson process. The base
`rate_sites_per_30_days` is measured per eligible lesion per 30 days; it is not
multiplied directly by cell count, lesion surface area, or `roots_per_event`.
At each sparse lesion refresh, the model computes the relief-adjusted occupied
fraction across that lesion's core blocks and maps it to a density stress:

```
stress = clamp((effective_occupied_fraction - stress_on) /
               (stress_full - stress_on), 0, 1)
multiplier = clamp(stress^density_exponent *
                   (biological_volume / volume_reference)^volume_exponent,
                   minimum_multiplier, maximum_multiplier)
rate_per_hour = base_rate_sites_per_30_days * multiplier / 720
```

The reusable defaults use exponent 1, zero volume exponent, and multiplier
range 0–4. The long single-cell production profile uses volume exponent 2/3,
which is the initial surface-area scaling for an approximately similar-shape
lesion, and a minimum multiplier of 0.25. The base parameter remains in
sites/30 days; the multiplier only changes its instantaneous Poisson
intensity. All exponents and clamps remain subject to scientific calibration.
`homogeneous_poisson` remains available as a test/backward-comparison strategy.

One unit-exponential hazard `H=-log(U)` is sampled per arrival from
`(seed, lesion_id, event_sequence)`. Rate changes consume
`elapsed_hours * old_rate_per_hour` from the remaining hazard and recompute the
deadline with the new rate. A zero rate pauses the clock without resampling;
later density recovery resumes the same pending hazard. This makes root
formation a continuous dynamic process instead of a one-time trigger or a new
random draw at every hourly refresh.

The uniform `U` is keyed by `(seed, lesion_id, event_sequence)`, so different
lesions have independent, reproducible waiting times. Each arrival is exactly
one attempted surface site and can create at most one root. `roots_per_event`
is a schema-v3 invariant and must equal 1. Multiple starting points arise from
successive arrivals and from multiple eligible lesions. Minimum separation,
the finite surface-candidate budget, local support, capsule collisions, global
and per-lesion root limits, and global and per-lesion active-tip limits can
reject an arrival.
`attempted_events`, `committed_roots`, and `rejected_events` therefore satisfy
`attempted_events = committed_roots + rejected_events`.

A rate of one site per 30 days means an expected attempted-site count of one
and a probability `1-exp(-1)` of at least one arrival; it does not guarantee a
root because the sampled arrival may be rejected.

If eligibility is lost, its pending seed event is invalidated and a new hazard
is sampled on re-entry. Occupancy changes mark only affected coarse blocks
dirty. A deterministic lesion-refresh event at the configured interval batches
those local scans; a seed event also forces current geometry validation before
committing. Checkpoints preserve the pre-refresh dirty observations, refresh
deadline/generation, stable core identity, and every lesion process, so resume
neither applies geometry early nor resamples history.

The tumour-wide process exposed for compatibility is a derived snapshot only.
Lesion records are ordered by stable lesion ID before checked accumulation.
Each active interval is folded through the snapshot time, and an eligible
aggregate is rebased to start at that same time. Counter overflow is a hard
error. Consequently insertion or hash-table iteration order cannot change the
aggregate or the simulation checksum.

## Root and tip construction

At event time, surface candidates are restricted to faces owned by the source
lesion. The sampler excludes closed-cavity walls: a face is external when its
outward-normal axis ray is unobstructed to infinity (the outermost
positive/negative face of its transverse lattice column). Subset extrema are
computed for that lesion, so another lesion cannot hide its boundary. A stable
hash top-K sample is uniform over the eligible external-face set without
copying or sorting the complete tumour surface.

The root position is `face.inside`, not the empty outside voxel. The root
capsule must replace at least one cell owned by the source lesion and leave the
configured number of neighbouring source cells outside the capsule. Candidate
ownership, inside occupancy, outside emptiness, vascular collision, separation,
and local support are revalidated immediately before atomic commit. One root
creates two independently scheduled tips carrying the same
`source_lesion_id`:

- outward, initially aligned with the exposed-face normal;
- inward, initially aligned toward the source lesion centroid captured at root
  birth.

Directions are selected from the fixed 26-direction set using a forward cone,
turn cone, persistence probability, and stateless RNG.  A step of Euclidean
length `L` at speed `v` completes after `L/v` hours.

## Geometry and occupancy

A centreline segment is rasterized as a capsule.  A voxel belongs to the vessel
when its centre is no farther than `diameter_voxels/2` from the segment.  The
complete capsule is proposed and committed atomically.

Outward growth remains empty-space while leaving its source lesion. If it
touches another independently indexed solid lesion, the arriving tip changes
to an inward branch, retargets that lesion's current centroid, adopts the
configured inward speed/length budget, and begins whole-cell displacement.
Inward growth may intersect cells and
removes every intersected biological cell as a whole: all eight voxels of a
large cell, a complete small cell, and all intersected co-located stage-2
cells.  Removal reason is `vascular_displacement`.  Every committed vessel
voxel permanently blocks later cell placement, movement, division, and stage
recovery.

An inward tip reaching its first centroid target enters `transiting` state and
keeps its established direction through the far half of the lesion. Its path
budget is not the old fixed 128-voxel constant. At root creation it is computed
from the greater of the centre-through distance and the farthest conservative
lesion-bound distance, multiplied by `length_tortuosity_factor`, then extended
by `exit_margin_voxels`. `max_length_voxels` is the minimum and
`hard_max_length_voxels` is the safety cap.

The production `continue_to_budget` far-surface policy keeps growing after the
tip first encounters an empty capsule. This both permits emergence outside the
far tumour surface and prevents an internal void from being mistaken for that
surface. The legacy `stop_complete` policy remains available for controlled
comparisons. A temporarily blocked tip stays event-driven and retries after the
configured interval. Tips still stop at their computed budget, domain
boundaries, or true vessel collision. Vessel collision defaults to anastomosis
and termination of the arriving tip. Branching and regression are disabled in
v1 but parent-network fields are retained.

## Vascular relief

Raw density always describes real cell anchors.  Vessels never rewrite the raw
density index.  A separate sparse relief field lowers effective crowding for
growth and density-mediated death:

```
s = distance_to_nearest_rasterized_perfused_vessel_voxel_center
relief = maximum_relief * max(0, 1 - s / cutoff_radius)
effective_density = raw_density * (1 - relief)
```

The cutoff comparison is strict (`s < cutoff_radius`). Because the vessel is a
rasterized capsule, the approximate maximum centerline reach is
`diameter/2 + cutoff_radius`: with the defaults this is `1.5 + 12 = 13.5`
voxels. At a source voxel the maximum relief is 0.50, so local effective density
is halved; it falls linearly to zero relief at the cutoff. Overlapping vessels
combine by maximum relief. The default scope excludes
occupancy and migration rules. Every generated vessel is perfused immediately,
so relief begins at root creation and extends as each inward/outward segment is
committed. The stored `perfused` compatibility field is therefore always one
for newly generated vessel nodes and tips.

Perfusion activation refreshes cells through the union of density blocks
covered by the real capsule voxels plus this cutoff, not through a single
network-wide bounding box. Each candidate anchor must also pass an exact
Euclidean cutoff test and have positive stored relief. Consequently, a long or
bent vessel does not refresh unrelated cells merely because they lie inside its
coarse AABB.

## Current schema-v3 defaults

The production YAML uses 8³ lesion blocks, 26-neighbour core connectivity,
occupied-fraction activation/deactivation thresholds 0.15/0.10, at least eight
cells per core block, a one-block attribution halo, a one-hour refresh interval,
and at least four core blocks. These are configurable numerical starting values
and require scientific/sensitivity calibration.

Each lesion enables seeding at 100,000 biological voxel³ and disables it below
80,000 voxel³. Its density-modulated base rate is 10 attempted sites per
eligible lesion per 30 days, with stress onset/full occupied fractions 0.15 and
0.60. One arrival attempts one root. Defaults include an 8-voxel
minimum root separation, 64 global/per-lesion roots, 128 global/per-lesion
active tips, inside-surface roots, and at least one surviving local support
cell.

Vessel diameter is 3 voxels. Inward tips grow at 0.50 voxel/hour through the
source lesion and outward tips at 2.0 voxel/hour toward the exposed-face
normal. The inward minimum budget is 128 voxels, with tortuosity factor 1.5,
16-voxel exterior margin, and hard cap 4096; the actual root budget is scaled
from its lesion. Blocked tips retry every hour, and an outward tip touching
another lesion converts to inward growth with a newly computed budget for that
lesion. An activated-r migration
rate is measured in moves/hour, while a fixed-26 move spans at most `sqrt(3)`
voxels. Configuration validation therefore requires outward speed to exceed
both inward speed and
`sqrt(3) * activated_r_rate_upper_bound`; the default `2.0` is strictly above
the current bound `sqrt(3) * 1.0`. Direction
and turn cones are both 45°, persistence
is 0.90, forward-direction weight is `exp(cos(angle))` at the default bias 1.0,
and no extra Euclidean step-length weighting is applied. The immediately active
linear influence has maximum density relief 0.50, decay length 4 voxels,
and cutoff radius 12 voxels. These new vascular numbers are parameterized
starting values rather than values inherited from the 2D model.
The `single_r_stage0_2160h_density_vascular_v5` experiment raises the activated-r
rate upper bound to 3 moves/hour and therefore uses an outward vessel speed of
6 voxels/hour, which remains strictly above `sqrt(3) * 3`.

## Determinism and persistence

Seed hazards, surface selection, tip directions, and conflicts use independent
stateless RNG event kinds keyed by seed, actor UID/lesion ID, event sequence,
and draw index. Output never consumes this RNG. Checkpoint base schema v6 and
the default stable-slot journal schema v8 store
stable lesion identity, exact dirty-block observations, every lesion process,
historical-source ownership aliases, the lesion-refresh scheduler, centreline
nodes and tips with their immutable source lesion, pending vessel steps, event
times, sequences, generations, remaining Poisson hazard, current rate, density
stress, and hazard integration time. Sparse cell and vessel occupancy and
relief layers are rebuilt and validated on restore. Schema-v7 field deltas
remain readable for older runs.
