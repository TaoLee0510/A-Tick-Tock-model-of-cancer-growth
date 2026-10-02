# Nutrient-coupled angiogenesis

ABM schema v3 accepts explicit optional model choices. The default outward
`speed_policy: strict_supremum_v1` preserves the existing path-speed check.
`mean_path_speed_v1` compares the mean migration rate, including beta lower
clamping, times the distance-weighted mean lattice jump length. This is a
nominal unblocked speed; guidance and spatial correlations can change realized
speed. `disabled` omits the active-r comparison. All policies still require
positive speeds and outward speed greater than inward speed.

The new seed model `hypoxia_modulated_poisson_v2` requires a nutrient environment.
The Poisson intensity scales with the fraction of a lesion's cells whose
normalized nutrient is below `hypoxia_threshold`, then applies the configured
volume and rate multipliers. Existing lesion volume/core eligibility still
applies. Nutrient ticks update remaining hazard without resampling it. In a
thin layer this new model samples only exposed planar faces, so outward tips
have feasible directions. Published seed models retain their original sampling.
Examples for r20 and r200 explicitly select the disabled speed ordering.

Continuum v6 and structured v8 add `angiogenesis.model: vegf_tip_density_v1`.
Hypoxic cell density produces VEGF/TAF; reflecting finite-volume diffusion and
exponential decay evolve it. Tip density moves by conservative upwind VEGF
chemotaxis and ordinary diffusion, with branching and anastomosis terms.
Hypoxia-modulated tips seed at the supported tumour surface; their total rate
is configured in tips/hour. In the comparison, two tips per ABM root and the
mean of inward/outward speeds set the continuum rate and deposition speed.
Vessel volume fraction grows by saturating deposition and supplies nutrient
through perfusion exchange. Fractions at the exclusion threshold become
Dirichlet sources and displace cells, matching the ABM's replacement rule at
coarse resolution. Substeps enforce the transport CFL.

The shared ABM environment evolves VEGF for diagnostics and uses actual ABM
vessels for perfusion. The PDE evolves its continuous tip/vessel densities.
All three vascular fields participate in new-version fingerprints, state
checksums and checkpoint formats; old checkpoints keep their formats.

`angiogenesis_v8.yaml` and `angiogenesis_r200_v8.yaml` run through the shared
executable. The standalone continuum/structured executables support fresh
hypoxic-profile initialization; import of a coupled ABM checkpoint needs the
shared executable and its resource sidecar. Use the ordinary checkpoint test
paths for direct PDE restart.

The `--vascular` ensemble report compares deposited vascular volume divided by
cross-section (a common length proxy), total perfused volume, and perfused
fraction inside the smoothed, connected, hole-filled 2D lesion. ABM additionally
reports actual parent-to-child centreline length. The 3D summary uses occupied
support for lesion volume; it does not implement the 2D front's hole filling.
The early 8-hour, 16-seed smoke tolerances are 75% relative length/volume and
0.15 absolute lesion perfusion, with reported paired sampling uncertainty.
These broad tolerances are a regression check, not calibration evidence.
Root rejection, vessel rasterization, directional persistence and tip-density
closures differ. The original comparison disables branching and anastomosis.
`angiogenesis_nonlinear_v8.yaml` adds a separate 16-seed comparison with both
mechanisms enabled; its activity guard requires actual ABM roots/anastomoses
and positive PDE branch/anastomosis rates. See
[activated_vascular_validation.md](activated_vascular_validation.md) for its
measured results and scope.
