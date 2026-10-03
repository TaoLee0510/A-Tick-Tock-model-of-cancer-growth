# Invasion-front hybrid v4

This contract and its acceptance margins are recorded before running the
invasion ensemble. `hybrid_invasion_front_v4` preserves the published v1-v3
models. The native all-ABM and all-PDE limits still use their own operators.

Adaptive v4 keeps every cell of phenotype r individual. This conservative policy retains
UID-specific activation, cooldown and hysteresis rather than replacing them
with a cell-mass clock. Both phenotypes at the invasion front remain agents;
only deep-core K cells are eligible for agent-to-density conversion. Ordinary core r
cells therefore remain ABM too. This is a deliberate restriction of the mixed
closure, not evidence that arbitrary r cohorts can cross the interface without
losing their joint history. Adaptive v4 rejects array initialization containing
r density and disables PDE activation internally. Native limits retain it.

Core classification uses smoothed combined occupied volume, distance from the
moving tumour boundary, an active-r exclusion halo, and nutrient gradients.
In 2D the boundary is the shared largest-connected, hole-filled mask. In 3D a
smoothed occupied-volume mask is used. Axis-neighbour graph distance is
conservative: a core location must lie beyond the configured front band;
entering the core also requires the extra distance hysteresis band. A location
enters below `gradient_off` and leaves at `gradient_on`. Occupancy uses the
existing core-on/core-off hysteresis. Every voxel of a converted large cell
must meet the core rule. Any active r protects its footprint and the configured
graph-distance halo independently of the occupancy mask.

New density queries use canonical 3D integral sums. Activation remains
quantized to the native query blocks and uses the exact edge denominator.
This changes floating-point summation only in v4. Derived prefix sums, masks
and distances are rebuilt on restart; core hysteresis flags and coverage extrema
are persisted. Conversion uses regional membership and vacant-volume caches,
with no repeated region scan per proposed individual. Mass below one cell stays
as density. Interface sites still obey the v3 biological volume constraint.
Density transport and daughter placement accept destinations only in the
classified core; large cohorts require their full footprint to lie there.
Transport reflects at that interface between exchanges. This new v4 closure
prevents diffusive fractional K tails from filling the agent front. It does not
change vessel geometry, consumers, nutrient boundary conditions or the volume
constraint. Native grid collision rules debit an existing swap partner or
retain the moving cell's overlapping footprint; the external capacity check
counts the resulting unit occupancy once. New placement still counts its
additional unit. Ultrasmall colocations are unsupported in v4, whose PDE
has only small and large compartments. Fractional residuals already outside the core can return to it or
convert at a subsequent exchange; they can still impede individual migration
and remain a stated interface approximation.
The strict representation invariant applies to every r cell, including active
r, while K cohorts below one cell or without an available full footprint can
remain as PDE interface remainders. `front_pde_mass` reports that K mass outside
the classified core. Thus this version does not claim mathematically zero PDE
mass throughout the front band. Whole-cell conversion and fractional-mass
conservation cannot in general satisfy that stronger invariant simultaneously.

The shared vascular PDE owns the single vascular process in adaptive v4.
The externally driven ABM environment has no duplicate tip process. Consumers
from both representations feed that process, and its exclusion diagnostics
record deletion in either representation.

## Prespecified invasion experiment

The first experiment uses the 256-square, 8-hour r20 activation configuration
from phase B. It exercises early invasion with an interior r99. It cannot
establish long-time invasion accuracy or production-scale suitability.
ABM A seeds 1-16, ABM B seeds 17-32, hybrid and PDE seeds 1-16 are compared at
0, 4 and 8 hours using the phase B TOST protocol and its unchanged margins:
mass and r/K 0.35 relative, active fraction 0.05 absolute, r50/r90 0.25,
r99 0.30, radial-profile L2 0.35. ABM is the primary reference; PDE is a
secondary comparison with the same bounds and the ABM baseline reported.
An additional disjoint PDE seed group supplies the secondary reference's
within-model baseline. Refinement baselines use disjoint hybrid realizations
of the fine scenario; the ABM baseline remains available in the primary table.
Each run must remain interior under the phase B r99/half-width < 0.8 guard.

Every adaptive realization must have minimum ABM mass fraction >= 0.20 and
minimum active fraction of r >= 0.10, evaluated at initialization and every
macro step, not just report times. It must also reach PDE mass fraction >= 0.10
and perform conversions in both directions, so an all-ABM run cannot satisfy
mixed coverage. Reported trajectories include mass fractions, cumulative
conversion counts and these extrema.

Splitting steps 0.25, 0.125 and 0.0625 hours and exchange intervals 2, 1 and
0.5 hours are retained. Refinement uses the previously declared 0.10 relative
scalar, 0.02 absolute active-fraction and 0.15 profile bounds, with paired TOST
instead of adding sampling allowance. Coarse/fine and middle/fine contrasts
are reported separately; monotonic drift is descriptive, not a proof of
convergence order. No failed bound is enlarged after observing the results.
