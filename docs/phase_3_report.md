# Phase 3 verification

Added `ATCG3D_ODE`, `atcg3d_ode`, ODE schema/model v1, a spatial-mean reaction
kernel with six cell compartments and active/refractory clocks, nutrient and
optional fixed vascular capacity, adaptive RK45, checkpoint/restart and a
periodic finite-volume verification mode. Existing schemas are unchanged.

The new assert-based CTest `atcg3d_ode_test` checks uniform ODE versus periodic
PDE to 2e-6, the independent published structured-v7 Euler reaction in its
small-step limit to 5e-5, nutrient exponential decay, stiff mortality,
nonnegative compartments, conservative activation/expiry/refractory changes,
nonuniform periodic mass/clock flux and bitwise restart with four threads.
The 48-hour example completes with total density 0.00998228.

Release checks: default 31/31 CTests and HDF5 32/32 CTests pass, including
the stiff mortality assertion and the legacy bitwise fixtures. The new periodic reference uses isotropic active diffusion;
the ODE's mean clocks and volumetric perfusion are closures of heterogeneous
age and Dirichlet-source distributions. These limitations are explicit in the
model documentation. Dynamic vascular capacity and angiogenesis are phase 4.
