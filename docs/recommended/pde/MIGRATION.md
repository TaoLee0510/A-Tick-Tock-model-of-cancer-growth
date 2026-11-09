# Configuration migration

Source files were left unchanged. References and output directories are relative.
New schemas use exact ABM growth windows and transported cohort refractory clocks;
upgrading v5/v6 changes the old grid-local cooldown closure. Restart a new run;
old checkpoint fingerprints are intentionally not reused by changed PDE contracts.
Structured schema 17 persists full population tails, nutrient work buffers
and moving-front workspace for exact continuation. Existing live operator
choices and published sector sums are preserved. Explicit hybrid v5 migration
initializes resource fields and individual clocks before core classification;
the v4 invasion-front policy is retained.
Nutrient v1 upgrades use the explicitly supplied per-cell K uptake rate
and preserve the old r/K uptake ratio and common saturation. This is a new
calibration choice, not an exact conversion of occupied-voxel uptake.
