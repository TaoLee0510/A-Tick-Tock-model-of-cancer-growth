# v5 calibration record

This record describes the completed implementation checks and the first spatial
calibration. It is not a 2160-hour production result.

## Automated validation

The full CTest suite passes: 25/25 tests. The v5-specific coverage includes
fixed planar-edge and vessel nutrient values, transient diffusion, equal
per-cell consumption across r/K and small/large stages, nutrient-only active-r
direction selection in a uniformly crowded field, clock expiry with refractory
hysteresis, conservative r/K exchange, vessel exclusion, and deterministic v5
checkpoint/restart.

## 256 x 256, 24-hour source controls

All runs start from `N(x,0)=1` and use the same cells and vessel mask.

| Fixed nutrient sources | Mean N at 24 h | Edge N | Vessel mean N | Vessel occupancy |
|---|---:|---:|---:|---:|
| Edges and vessels | 0.836160 | 1.000000 | 1.000000 | 0 |
| Edges only | 0.831510 | 1.000000 | 0.841314 | 0 |
| Vessels only (zero-flux outer boundary) | 0.832077 | source-dependent | 1.000000 | 0 |

The simultaneous-source run has centre nutrient 0.709933 at 24 hours. Its
assembled consumption divided by cell number is exactly 0.01 per hour for both
r and K.

## 512 x 512, 360-hour calibration

Run directory:
`atcg3d_structured_pde_nutrient_chemotaxis_calibration_2d_512_360h_v5_run`.

| Time (h) | r cells | K cells | Active-r fraction | r mean radius | K mean radius | r radius 90 | Mean N |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 4,130.5 | 4,129.5 | 0.072146 | 34.7 | 35.2 | 51 | 1.0000 |
| 24 | 2,854.4 | 3,770.7 | 0.000224 | 38.5 | 37.2 | 57 | 0.8461 |
| 120 | 2,153.9 | 6,275.5 | 0 | 94.9 | 40.5 | 191 | 0.4454 |
| 240 | 21,251.9 | 7,098.1 | 0 | 152.1 | 42.9 | 230 | 0.2069 |
| 360 | 73,458.4 | 4,624.8 | 0.00000193 | 170.1 | 53.8 | 244 | 0.08715 |

At 360 hours, 88.38% of r but only 4.93% of K is outside radius 100. Inside
radius 60, only 0.01% of r remains while 74.13% of K remains. The outer region
is therefore r-dominant (local r fraction 0.9965), while the lesion core is
K-dominant. Vessel nutrient is exactly 1 and vessel occupancy is exactly zero.

The active clock is not permanent: the initial active fraction falls from
7.21% to 0.022% by 24 hours and to zero by 120 hours. A very small later episode
can re-arm only after both the 24-hour cooldown and the 0.80 off-threshold have
cleared.

## 2000 x 2000, 720-hour production calibration

The desired r/K spatial separation in the 512-grid run is reproduced without a
boundary artifact on the full plane. The PDE evolved all four million sites;
field CSV output was sampled every four sites to keep its size practical.

Run directory:
`atcg3d_structured_pde_nutrient_chemotaxis_2d_2000_720h_v5_run`.

Final exact diagnostics:

| Quantity | Value |
|---|---:|
| r cells | 15,444.9234 |
| K cells | 2,247.1643 |
| Active-r fraction | 0 |
| r mean radius | 236.382 |
| K mean radius | 64.238 |
| r radius 90 / 99 | 329 / 396 |
| Mean / maximum nutrient | 0.020492 / 1.0 |
| Maximum occupied fraction | 0.999880 |
| State checksum | 5508321473299096198 |

The stride-4 final field gives the following spatial audit:

- Outside radius 100: 93.62% of r and 13.11% of K; local r share 97.78%.
- Outside radius 200: 70.30% of r and 4.72% of K; local r share 98.92%.
- Inside radius 60: 0.005% of r and 61.83% of K; local r share 0.051%.
- Within 10 voxels of the vessel axis: local r share 78.96%.
- Outermost 20 grid layers: no r and no K; `r99=396` is far from the
  1000-voxel half-width.
- Vessel nutrient is exactly 1 and vessel occupancy is exactly zero.

Cell mass peaks earlier and declines once the common resource is depleted;
mean nutrient is 0.0205 by 720 hours. This shows that the shared nutrient field
limits total population without assigning r and K different demand or carrying
capacity. The full-plane and central-window figures are
`figures/density_720h.png` and `figures/density_720h_zoom450.png` inside the run
directory.

## 2000 x 2000 continuation to 2160 hours

The accepted 720-hour checkpoint was continued with an identical dynamics
fingerprint to 2160 hours. The restart is therefore the exact continuation of
the preceding run, not a reinitialized experiment.

Run directory:
`atcg3d_structured_pde_nutrient_chemotaxis_2d_2000_720_to_2160h_v5_run`.

| Quantity | Value |
|---|---:|
| r cells | 5,495.1841 |
| K cells | 10,757.4493 |
| Active-r fraction | 0 |
| r mean radius | 446.307 |
| K mean radius | 192.458 |
| r radius 90 / 99 | 540 / 587 |
| Mean / maximum nutrient | 0.014075 / 1.0 |
| Maximum occupied fraction | 0.999944 |
| State checksum | 5647931471677280949 |

The stride-4 spatial audit shows that r remains the peripheral phenotype:

- Outside radius 400, the local r share is approximately 92.0%.
- Outside radius 500, the sampled population is effectively all r.
- Inside radius 100, local r share is only about 0.18%.
- No cells occur in the outermost 20 grid layers; `r99=587` remains far from
  the 1000-voxel half-width.
- Vessel nutrient is exactly 1 and vessel occupancy is exactly zero.

However, the long-time geometry is no longer a broad two-dimensional r ring.
The population contracts into a narrow vessel-aligned corridor: approximately
57.8% of sampled r and 76.3% of sampled K lie within 10 voxels of the vessel
axis, where the local r share is only 26.6%. K therefore again dominates the
immediate vascular corridor, while r occupies the distal ends and larger
radii. This is not a boundary artifact. At `D_N=0.1`, shared consumption and
decay reduce mean nutrient to 0.014; over this duration the distant planar-edge
sources do not sustain the central two-dimensional interior.

Thus the run satisfies the numerical and stopping constraints, and r is
radially peripheral, but it does not satisfy a stronger requirement that r
remain broadly distributed across the two-dimensional periphery at 2160
hours. That stronger target would require another biological calibration, most
directly of nutrient transport/source geometry or the post-activation survival
mechanism.

The full-plane and central-window figures are `figures/density_2160h.png` and
`figures/density_2160h_zoom700.png` inside the continuation run directory.
