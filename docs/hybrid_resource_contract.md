# Hybrid v5 resource initialization and continuation

`hybrid_resource_restart_v5` retains the v4 invasion/front policy and requires
structured schema v17. It initializes the agent environment through the same
resource assembly as standalone shared ABM before transferring its nutrient
field into the coupled density solver. Initializing an empty PDE first makes
every voxel exterior to the lesion; in a moving-front configuration that
previously replaced a requested low initial nutrient by the Dirichlet value.
Version 5 preserves the cell-limited initial field and clocks. V4 remains
available with its published initialization.

Structured v17 native checkpoints retain the complete positive population
tails and active totals, nutrient work buffer, moving-front mask, its bounds
and update box. Traversal cutoffs do not prove those values irrelevant:
vascular source assembly visits the whole grid. Earlier checkpoints could
therefore restore their hash and still diverge in a later vascular update.
Version 17 hashes the complete saved state. It preserves the uninterrupted
kernel arithmetic; old schema formats and hashing remain unchanged.

Hybrid configurations using structured v17 reconstruct guidance/activation queries without
recanonicalizing saved front bounds or active totals. Derived FFT guidance
caches remain outside checkpoints. Cross-process, one-to-eight-thread tests
exercise dynamic vascular fields, including low initial nutrient, in ABM,
PDE and hybrid. Initialization tests compare standalone ABM and hybrid
nutrient, vessels and individual r clocks exactly in 2D and 3D.
The restart test also covers a migrated v4 front wrapper with structured v17;
its published high-nutrient initialization is retained, while the selected
PDE persistence policy preserves saved workspaces during restoration.

The current v5 recommendation uses the unchanged 256-square, 8-hour r20 case
with observations at 0, 4 and 8 hours. `validate_hybrid_resources.py` uses
16 paired hybrid/ABM realizations and 16 disjoint ABM baseline realizations.
It uses the previously declared phase B TOST margins, alpha and profile
jackknife method without an additional allowance. Every realization also
requires the phase D representation minima and interior screening condition.
The native ABM baseline may reuse the preceding local invasion validation
files only after the dry-run descriptions select identical live operators and
parameters; schema metadata, file/output names and unchanged default sector
metadata are excluded from that comparison. If the baseline files are absent,
the validator runs those native seeds itself.

The four splitting/exchange refinement ensembles in phase D apply to v4.
The v5 experiment adds a primary comparison and exact resource/restart tests;
it does not claim a new low-nutrient, growing vascular multi-seed equivalence
study or a fresh four-case v5 refinement study. The production configuration
and performance prefixes stay separately qualified in their protocol/report.
