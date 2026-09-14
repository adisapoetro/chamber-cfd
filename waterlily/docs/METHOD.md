# Airflow and CO₂ transport

WaterLily advances incompressible velocity and pressure on a staggered Cartesian
grid. Its predictor/corrector, pressure projection, immersed-boundary kernels
and Smagorinsky forcing helper are called directly. No Python momentum solver
or replacement Julia pressure algorithm is used.

$$\nabla\cdot\mathbf u=0,$$

$$\frac{\partial\mathbf u}{\partial t}+\nabla\cdot(\mathbf u\mathbf u)
=-\frac{\nabla p}{\rho}+\text{viscous and subgrid terms}
-K(\mathbf x)|\mathbf u|\mathbf u+\mathbf f_{fan}/\rho.$$

The pressure wrapper changes only tolerance and maximum iterations for the
native multigrid solver and raises an error when it does not converge.
The relevant package sources are
[Flow.jl](https://github.com/WaterLily-jl/WaterLily.jl/blob/v1.8.0/src/Flow.jl),
[Body.jl](https://github.com/WaterLily-jl/WaterLily.jl/blob/v1.8.0/src/Body.jl)
and [util.jl](https://github.com/WaterLily-jl/WaterLily.jl/blob/v1.8.0/src/util.jl).

## Physical and grid units

WaterLily uses unit-width grid cells. This adapter keeps velocity numbers in
m/s and converts distance by the physical cell width h:

$$\mathbf x_{grid}=(\mathbf x_{SI}-\mathbf x_{lower})/h,\qquad
t_{grid}=t_{SI}/h,\qquad \nu_{grid}=\nu_{SI}/h.$$

The reference velocity unit is 1 m/s. Acceleration supplied to WaterLily is
`h × acceleration_SI`. Body velocities therefore have the same numerical values
as m/s. Output timestamps convert `WaterLily.time(sim)` back to seconds; they
are not animation time or `sim_time` without its length/velocity conversion.

## Walls and opening

Two hollow shells move along x, each by half the configured gap. Distance,
normal and wall velocity are supplied through WaterLily's `AbstractBody`
interface. The nominal inside dimensions remain 4 × 6 × 4 m in the reference
case. The external ground is z=0 and the chamber floor is z=0.6 m.

At full closure the joined enclosure has a continuous signed-distance field,
including across the roof/floor seam. During movement the two open-ended halves
are measured separately. The native body data are refreshed as the walls move.
Time steps stop at operating events and requested output times.

A three-cell numerical wall provides a zero-mobility core with the one-cell
immersion kernel. This extends the solid outward and changes exterior blockage.
It does not resolve a 3 mm physical skin on a 0.5 m grid. Refinement and a
qualified thin-wall treatment are needed before comparing results quantitatively
with the zero-thickness fitted-wall reference.

Native external boundary conditions prescribe the normal wind components and
use zero tangential gradients; the special positive-x convective exit is
disabled so the same convention accommodates different wind directions.
These differ from the Python reference's exterior boundary treatment.

## Fans, plant drag and unresolved mixing

Each fan supplies a nominal momentum source:

$$Q=0.0004719474432\,Q_{CFM},\qquad A=\pi d^2/4,\qquad
F=\rho Q^2/A\;m.$$

Here m is `momentum_factor`. Gaussian forcing is sampled at velocity faces,
normalized over its chamber support and supplied through WaterLily's forcing
hook. This is an assumed thrust scale, not a measurement of installed flow.
The native immersion operation still acts on momentum near walls. No air or
CO₂ is created by the fans. Their support and housings are not resolved solids.

The prescribed crown and trunk drag envelopes retain the reference dimensions
and coefficients. A schematic display palm is scaled from the same reference
geometry. Its triangles do not represent a measured resolved canopy.

The native `sgs!` helper uses the standard Smagorinsky viscosity function
`(Cs × Δ)² × sqrt(2 S:S)`, with `Cs=0.12` and a filter width of one grid cell.
Explicit time steps also account for drag and this added diffusion. This coarse
calculation is not an established turbulent-flow solution or a resolved DNS.

## Outputs and limits

Velocity/distance arrays contain little-endian Float64 values in Julia's
column-major order on the reference platform. The summary supplies their shapes,
paths and physical grid origin. CSV columns label units explicitly.

The fixed reporting weights are exact cell intersections with the nominal
96 m³ box. That makes the reporting volume independent of grid registration.
These weights are separate from the scalar's exact fluid cut-cell volumes.
The connected closed scalar compartment is checked to contain 96 m³. Native
pressure residual, divergence and solid-core velocity remain separate checks.

## CO₂ transport

Air carries CO₂ by advection. Molecular and modelled subgrid diffusion spread it:

$$\partial_t C+\nabla\cdot(\mathbf u C)
=\nabla\cdot[(D_m+\nu_t/Sc_t)\nabla C]+s_C.$$

Fans change velocity; they do not create gas or directly homogenize CO₂. Ambient
concentration enters through inflow faces. Outflow carries local concentration
with zero diffusive gradient. Walls and ground have zero scalar flux. There is
no ambient reset when the chamber opens.

The extension integrates inventory using exact intersections with the same
rectangular shells used for airflow:

$$V_i^{n+1}C_i^{n+1}-V_i^n C_i^n
+\Delta t\sum_f[Q_f C_{upwind,f}^{n+1}-G_f(C_j^{n+1}-C_i^{n+1})]
=\Delta t S_i.$$

V is fluid volume (m³), Q signed air volume flux (m³/s), and
`G = open face area × harmonic-mean diffusivity / cell spacing` (m³/s).
Backward Euler and first-order upwinding give a sparse, monotone system.
Small cut cells are handled implicitly without merging, resetting or clipping
concentration. First-order transport introduces numerical diffusion.

Opening and closing use time-averaged apertures. Intervals split whenever a
moving rectangular edge crosses a grid face. Midpoint integration is exact
for the resulting piecewise-linear area functions. Intersection roundoff below
10⁻¹³ of a cell is treated as geometric zero.

The face flux satisfies the discrete geometric conservation law:

$$V_i^{n+1}-V_i^n+\Delta t\sum_f Q_f=0.$$

Covered and uncovered cells exchange inventory through shared fluxes. A uniform,
source-free concentration stays uniform as walls move; exposed cells do not
receive an arbitrary initial concentration.

## Mapping native BDIM velocity to scalar flux

WaterLily's immersed velocity blends fluid and body motion. The extension
removes the body contribution, divides by native face mobility and integrates
over the sharp aperture:

$$Q_{raw}=A_{open}\,[u_{BDIM}-(1-\mu_0)V_{wall}]/\mu_0.$$

This is an interface approximation: it does not invert every first-moment BDIM
term. Velocity is taken after the native step; body geometry for deblending is
sampled at its midpoint. A weighted least-squares flux correction enforces the
geometric law on the scalar fluid graph. It does not update WaterLily momentum
or replace its pressure solver. Sealed components must have compatible volume
balances before fixing a gauge; no mass defect is silently redistributed.

Each step records corrections globally, in the nominal chamber region and on
partial faces, plus the largest velocity change. A 10% chamber-region correction
screen is diagnostic, not a validation threshold. Failure remains visible in
the summary and prevents quantitative mixing claims.

## Leaflet net exchange

Source-tagged triangles are scaled to the palm envelope and clipped against
each cell. Only the 240 leaflet triangles contribute. Their approximately
0.733 m² area is a representative distribution stencil, not measured total
leaf area. Trunk and bare frond axes have no source.

$$n_{air}=P/(RT),\qquad
S_i={F_{net}\over n_{air}}\,{A_{leaflet,i}\over\sum_j A_{leaflet,j}}.$$

F is µmol/s and S is ppm m³/s. The sum equals the prescribed whole-tree rate at
every step. Uptake continues during opening, full exposure and closing. Since
the input is net exchange, adding separate respiration would double-count it.
A sink that exhausts local CO₂ raises an error instead of clipping or silently
reducing uptake.

For the sealed reference chamber, mean CO₂ declines at
`F_net / (n_air × 96) = −0.011315371 ppm/s`: approximately **3.395 ppm in five
minutes**, independent of mixing. Spatial variation depends on airflow,
diffusivity, source distribution and numerical resolution.

The geometric conservation requirement follows finite-volume principles
discussed by [Seo and Mittal (2011)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3156558/).
This project's coupling and tests are its own implementation, not a reproduction
of that paper's complete momentum method.
