# Shared VEGF-guided angiogenesis contract

This phase uses option (a): individual ABM tips follow the VEGF field and use
the same local transport, branching and connection rates as the PDE. The
published forward-persistent ABM walker and `vegf_tip_density_v1` PDE retain
their original arithmetic. The new model is `shared_vegf_lattice_v2`, selected
only by structured schema 15 or a later explicitly compatible schema.

## State and operators

VEGF production is the cell count per voxel, including footprint-distributed
large cells, multiplied by `production * clamp(1-N/(Nmax*threshold),0,1)`.
VEGF diffuses on reflecting faces and decays exponentially. Both models use
the same field solver, source assembly and operator ordering. Individual tips
are separate agents with persistent identifiers and a checkpointed random
counter. They do not use the published forward-cone walk.

For each axial face from i to j, the jump rate is
`q(i,j) = Dtip/h^2 + chi*max(T[j]-T[i],0)/h^2`.
The PDE uses the corresponding conservative face flux. A vascular substep
satisfies `dt*max(sum(q),2*d*Dvegf/h^2) <= 0.45`. ABM tips select one jump or
stay using exactly those probabilities; the PDE propagates their conditional
first moment. This is a discrete-time velocity-jump approximation to the
continuous-time process. Its substep refinement must be measured.

Movement precedes branching. The branching rate per tip is `b*T`; a Yule
offspring draw has geometric survival parameter `exp(-b*T*dt)` and the PDE
uses its expectation. Local connection survival is
`exp(-mu*(old_vessel_fraction + moved_tip_density)*dt)`. This coarse local
connection law includes tip crowding; a terminated tip is counted once. It
does not identify a unique geometric junction or distinguish self-contact.
That limitation is shared and must accompany any interpretation of the
anastomosis count. New tips are seeded after connection, with Poisson mean
`dt*seed_rate*hypoxic_cell_fraction`, on hypoxic lesion surface sites weighted
by their hypoxic consumer count. Seeding falls back to all hypoxic sites when
the surface weight is zero. Branch and connection counts are cumulative
event counts in the ABM and cumulative expected counts in the PDE.

Unlike the published model, the new law has no hard clipping of tip density:
clipping an expectation has no corresponding integer-agent transition law.
`maximum_tip_density` is a legacy parameter for this model. The code must
reject nonfinite or excessive individual populations instead of silently
clipping them. `tip_speed_voxels_per_hour` is also a legacy parameter: lattice
speed follows the declared face rates, so deposition is derived from actual
traversal rather than an independently configured speed.

## Centerlines, perfusion and exclusion

Each undirected lattice edge has a deposited occupancy. ABM traversal marks
an edge once; PDE occupancy follows the independent-arrival closure
`e_new = 1-(1-e_old)*exp(-dt*face_traversal_rate)`.
New centerline length is `h*sum(e)`. This counts unique deposited lattice
segments, including branches, rather than converting volume to length.
It is not the sum of repeated path traversals. The nonlinear independent
arrival closure can differ from correlated ABM revisits; validation measures
that difference explicitly.

New edge occupancy deposits half its tube-volume dose at each endpoint;
voxel perfusion uses the union law `1-(1-v)*exp(-dose)`. Static source voxels
retain their prescribed perfusion. Static voxel masks do not identify a
unique centerline, so the primary length metric is new centerline growth;
initial source volume is reported separately. A synthetic line's analytical
length may be reported separately but must never be substituted for growth.
Perfusion supplies nutrient through the same exchange law on both sides.
Each seeded root is treated as a perfused supply inlet from unresolved host
vasculature. This is a prescribed boundary assumption; the new lattice model
does not determine connectivity to an explicit external vascular network.
Perfusion evidence is conditional on that common assumption.
Configured vascular exclusion prevents tumor entry and removes overlapping
resident tumor mass; such removal is diagnosed and is not displacement.

The new solver reuses buffers and updates a bounding box containing consumers,
VEGF, tips and newly deposited edges, with a one-site transport halo. A
versioned activity cutoff of 1e-14 bounds numerical support; discarded VEGF
and tip mass are accumulated as diagnostics. It must not truncate static
source perfusion. Parallel loops write disjoint sites; fixed spatial
reduction order and serial identifier order preserve thread-independent
arithmetic and checkpoints.

## Prespecified verification

Before executing the new ensembles, length and perfused-volume relative
margins remain 0.75 and lesion perfusion remains 0.15 absolute, as declared in
the earlier vascular validation. New cumulative branch and connection count
margins are 0.75 relative: these use the same broad mechanism-level screening
bound as centerline growth, rather than a clinical calibration claim. No
new margin exceeds the previously declared vascular bound. Cumulative seed
count uses the same new 0.75 relative screening margin to check the common
hypoxia-dependent Poisson seeding law. No
sampling allowance is added. These broad margins are exploratory and do not
establish accurate vascular morphology.

Each ensemble uses at least 256-square voxels, 240 hours, paired seeds 1-16,
and independent ABM baseline seeds 17-32. Samples occur every 24 hours. All
positive-reference vascular contrasts must pass the phase B paired TOST
criterion at alpha 0.05. Structural zero counts at initialization are reported
and verified by exact initialization tests; they are outside relative-margin
inference. A zero reference at any post-initialization sample is undefined
and fails equivalence. Full scalar realizations, baseline differences, paired
t statistics, both TOST p values and every time table are retained. Controlled
frozen-lesion experiments must be labeled as mechanism tests; they do not
establish agreement in a growing invasive tumor.

The controlled configuration disables tumor migration, activation, division,
exchange and vascular exclusion, and uses positive growth limits to prevent
death. This fixes the tumor consumer realization while perfusion and nutrient
still interact. Exclusion is checked separately with overlapping-cell tests;
it remains enabled in growing production configurations. In a thin layer,
the extra z=1 plane of an ABM large footprint maps to the physical z=0 resource
plane for entry/exclusion checks. Large-cell removal remains a correlation
closure: ABM removal is indivisible, while PDE removal can remove only the
overlapping footprint mass.

Unit checks compare jump conditional expectations, VEGF gradient bias,
geometric offspring moments, repeated-edge length, source perfusion, event
counts, cutoff loss, two/three-dimensional nonnegativity and checkpoint
restart. One/eight-thread fields must be bitwise equal. Benchmarks use a
2000-square grid and report wall time, peak resident memory, active range,
buffer allocation and checksums for both thread counts.

The isolated offspring/connection moment check initializes 200000 tips,
uses analytic Yule or binomial variance, and accepts an absolute deviation
below six theoretical standard deviations. This fixed unit-test sampling
bound is specified before the draw and does not modify ensemble margins.

Callers may supply a validated consumer bounding box. This declares that all
consumer entries outside it are zero, and avoids full-grid source scans; it
does not restrict VEGF or tip support to that box. Source assembly computes
that declaration from every nonzero consumer, including external hybrid
agents. The full-scan fallback and supplied-box paths must have identical
arithmetic and restart results.
