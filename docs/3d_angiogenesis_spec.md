# ATCG3D dynamic angiogenesis specification

## Eligibility and seed process

Angiogenesis is a continuing stochastic process, not a one-shot action.  It is
eligible while biological tumour volume is at or above an activation threshold
and becomes ineligible below a configurable deactivation threshold.  Existing
vessels continue to grow when seed generation is ineligible.

Biological tumour volume is the sum of per-cell stage volumes.  Vessel voxels
are excluded.  Stage-2 co-located cells contribute separately.

The configured `rate_sites_per_30_days` is a global Poisson intensity of
attempted surface-root sites for the whole tumour during eligible time; it is
not multiplied by surface area. With 30 days equal to 720 hours,

```
rate_per_hour = rate_sites_per_30_days / 720
waiting_hours = -log(U) / rate_per_hour
```

Each Poisson arrival is exactly one attempted surface site and can create at
most one root. `roots_per_event` remains in schema v2 for an explicit invariant,
but must equal 1; it is not a multiplier. Multiple starting points arise only
from multiple arrivals as the process continues through time. Minimum
separation, the finite surface-candidate budget, capsule collisions,
`max_total_roots`, and the two-active-tips budget can reject an arrival.
`attempted_events`, `committed_roots`, and `rejected_events` therefore satisfy
`attempted_events = committed_roots + rejected_events`.

A rate of one site per 30 days means an expected attempted-site count of one
and a probability `1-exp(-1)` of at least one arrival; it does not guarantee a
root because the sampled arrival may be rejected.

If eligibility is lost, the pending seed event is invalidated.  A new waiting
time is sampled on re-entry.  Checkpoints preserve eligibility, next seed time,
event sequence, and generation so resume does not resample history.

## Root and tip construction

At event time, the sampler first excludes closed-cavity walls: a face is
external when its outward-normal axis ray is unobstructed to infinity (the
outermost positive/negative face of its transverse lattice column). A stable
hash sample is uniform over that explicit external-face set. Configurable
minimum root separation against existing roots and a finite attempt limit
prevent duplicate roots and unbounded rejection loops. Candidate order is keyed
by the site-arrival sequence, so selection is deterministic and independent of
worker thread order. One root
creates two independently scheduled tips:

- outward, initially aligned with the exposed-face normal;
- inward, initially aligned toward the tumour centroid captured at root birth.

Directions are selected from the fixed 26-direction set using a forward cone,
turn cone, persistence probability, and stateless RNG.  A step of Euclidean
length `L` at speed `v` completes after `L/v` hours.

## Geometry and occupancy

A centreline segment is rasterized as a capsule.  A voxel belongs to the vessel
when its centre is no farther than `diameter_voxels/2` from the segment.  The
complete capsule is proposed and committed atomically.

Outward growth is empty-space only.  Inward growth may intersect cells and
removes every intersected biological cell as a whole: all eight voxels of a
large cell, a complete small cell, and all intersected co-located stage-2
cells.  Removal reason is `vascular_displacement`.  Every committed vessel
voxel permanently blocks later cell placement, movement, division, and stage
recovery.

Tips stop at configured length/node/voxel limits, domain boundaries, targets,
or when no direction is feasible.  Vessel collision defaults to anastomosis
and termination of the arriving tip.  Branching and regression are disabled in
v1 but parent-network fields are retained.

## Vascular relief

Raw density always describes real cell anchors.  Vessels never rewrite the raw
density index.  A separate sparse relief field lowers effective crowding for
growth and density-mediated death:

```
d = max(0, distance_to_centerline - vessel_radius)
relief = maximum_relief * max(0, 1 - d / influence_radius)
effective_density = raw_density * (1 - relief)
```

Overlapping vessels combine by maximum relief. The default scope excludes
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

## Current legacy-mapped-v2 defaults

The production YAML enables the process at 100,000 biological voxel³ and
disables new root seeding below 80,000 voxel³. Its global homogeneous Poisson
rate is 10 attempted surface sites per 30 days, at most one root per event, with an
8-voxel minimum root separation, at most 64 roots, and at most 128 active tips.

Vessel diameter is 3 voxels. Inward tips grow at 0.50 voxel/hour toward the
root-time tumour centroid and outward tips at 0.25 voxel/hour toward the
exposed-face normal; each is limited to 128 voxels. Direction and turn cones
are both 45°, persistence
is 0.90, forward-direction weight is `exp(cos(angle))` at the default bias 1.0,
and no extra Euclidean step-length weighting is applied. The immediately active
linear influence has maximum density relief 0.50, decay length 4 voxels,
and cutoff radius 12 voxels. These new vascular numbers are parameterized
starting values rather than values inherited from the 2D model.

## Determinism and persistence

Seed times, surface selection, tip directions, and conflicts use independent
stateless RNG event kinds keyed by seed, actor UID, event sequence, and draw
index.  Output never consumes this RNG.  Checkpoint schema v2 stores the
angiogenesis controller, centreline nodes, active tips, pending steps, event
times, sequences, and generations; sparse occupancy and relief layers are
rebuilt and validated on restore.
