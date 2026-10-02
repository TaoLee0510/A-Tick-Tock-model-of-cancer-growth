# Phase 4 verification

ABM adds explicit strict/mean/disabled outward speed strategies and a nutrient
hypoxia-modulated Poisson seed model. Defaults omit their new identity keys
and retain published checkpoint fingerprints. Thin-layer hypoxic seeding now
uses planar surfaces; this resolves outward tips blocked by a vertical bias.
Continuum v6 and structured v8 add hypoxic VEGF production/diffusion/decay,
upwind tip chemotaxis, branching/anastomosis, saturating vessel density,
perfusion sources and thresholded cell displacement. The shared ABM environment
uses actual vessels for nutrient and evolves VEGF diagnostics. New PDE and
paired ABM checkpoints include the vascular fields.

New assert-based `atcg3d_vascular_field_test` covers speed policies (including
clamped-beta mean), hypoxic versus oxic root commitment, planar tip bias,
VEGF spatial production, nonlinear nonnegative tips, vascular serialization,
continuum/structured restart, and resource/ABM restart with four threads.
HDF5 builds additionally check checkpoint reconstruction of the hypoxic model.
The `atcg3d_vascular_validation` CTest uses 16 paired seeds for eight hours.
Default 33/33 and HDF5 34/34 CTests pass, including both validation targets and
all published PDE bitwise checksum fixtures. A 4-hour hypoxic HDF5 checkpoint
resumed with four threads reproduces the complete 8-hour seed-3 report bitwise:
ABM checksum 15161047124841193578, resource checksum 15826928194346496749.

Vascular mean length proxy and perfused volume are 26.69 ABM versus 33.68 PDE,
relative error 0.2621. Lesion perfusion is 0.0375 versus 0.0902, absolute error
0.0527. Predeclared early-smoke tolerances are 0.75 relative and 0.15 absolute,
plus reported paired sampling uncertainty. They were not enlarged after the
initial failure; instead the new thin-layer root-direction defect was fixed.
The r200 hypoxic profile passes CLI loading; the standalone continuum profile
completes with checksum 1223917129674614749.

Remaining limits: continuous tips do not retain ABM directional history or
individual root rejection. Rasterized volume-derived length is an approximation
to centreline length, which is also reported separately for ABM. This broad,
early-time ensemble check is not a vascular calibration. The comparison disables
branching/anastomosis; their nonlinear solver checks establish bounds, not ABM
statistical equivalence. The 3D lesion-perfusion summary uses occupied support
rather than the 2D connected-front mask. Hybrid coupling remains phase 5.
