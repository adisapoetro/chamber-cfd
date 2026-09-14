# Model formulation

The model represents two sliding chamber halves around a schematic oil palm.
The reference enclosure is nominally **6 × 4 × 4 m**, with a floor 0.6 m above
external ground. Solver coordinates are x × y × z = 4 × 6 × 4 m; x follows
sliding travel. Three fans sit on a fixed support 0.30 m beside the trunk,
at 0.667, 2.000 and 3.333 m above the floor. They blow toward 0°, 120° and 240°.

The implementation uses NumPy/SciPy finite-volume operators. Configuration is
compiled by `fluid_dynamic.operating`; airflow and transport are calculated by
`fluid_dynamic.recreated.fan_solver` on the fitted geometry provided by
`fluid_dynamic.scenarios.geometry`.

## Air carries CO₂; diffusion spreads gradients

Let c be CO₂ mole fraction, u air velocity, D_eff effective diffusivity, and S
the local source. At constant reference pressure and temperature:

$$\nabla\cdot\mathbf u=0,$$

$$\frac{\partial\mathbf u}{\partial t}+\nabla\cdot(\mathbf u\mathbf u)
=-\frac{\nabla p}{\rho}+\nabla\cdot[\nu_{eff}(\nabla\mathbf u+\nabla\mathbf u^T)]
-K|\mathbf u|\mathbf u+\mathbf f_{fan}/\rho,$$

$$\frac{\partial c}{\partial t}+\nabla\cdot(\mathbf u c)
=\nabla\cdot(D_{eff}\nabla c)+S,$$

$$\nu_t=(C_s\Delta)^2|S_{strain}|,\qquad D_{eff}=D_{molecular}+\nu_t/Sc_t.$$

Advection moves gas with the computed wind/fan flow. Diffusion acts on
concentration gradients. The Smagorinsky term represents unresolved mixing;
it is a modelling assumption. Porous canopy/trunk drag slows air locally.
Pressure projection enforces incompressibility. CO₂ is passive, so its uptake
strength does not feed back into airflow, temperature or buoyancy.

## Opening is a moving boundary

The two ideal zero-thickness shells move apart along x. Each travels half the
configured gap. The mesh follows their surfaces (an arbitrary Lagrangian–Eulerian,
or ALE, mesh). Cell volumes change and some cells appear/disappear as the gap
opens/closes. Conservative inventory updates use relative air/mesh face fluxes:

$$\frac{d(V_i c_i)}{dt}+\sum_f[(\mathbf u-\mathbf w)\cdot\mathbf n A]_f c_f
=\sum_f(D_{eff}\nabla c\cdot\mathbf n A)_f+V_i S_i.$$

Here w is mesh velocity. Gas enters through computed transport across the
opening and external boundaries; no fitted mixing-time formula fills the
chamber with ambient gas. Exterior wind enters faces selected by its direction;
other exterior faces permit outflow/backflow with ambient inflow concentration.
The ground is stationary. Shell surfaces impose moving-wall velocity and block
normal scalar transport. Wall motion and output times split time steps exactly.
The scalar scheme is implicit donor transport with local flux correction;
source half-steps bracket transport. Time integration remains first order.

## Fans add momentum, not air or CO₂

For each fan, the assumed force scale is

$$Q = Q_{CFM}\,0.0004719474432,\quad A=\pi d^2/4,\quad
F=\rho Q^2/A\;m,$$

where m is `momentum_factor`. The force acts along the configured unit direction.
A Gaussian, integrated over each cell and normalized inside the sealed chamber,
distributes it spatially. It adds no volume or scalar source. Three directions
120° apart can sum to zero resultant force while still producing local circulation.
The default force is about 0.513 N per fan. This is a nominal free-jet scale;
installed thrust, speed setting, mounting effects and fan pressure-flow curves
have not been measured here. Fan discs, pole, braces and housings are not resolved
solid obstacles. Fans run only while fully sealed.

## Exchange is on leaflets along the fronds

The palm illustration has a simple trunk and five fronds.
The source uses 120 representative leaflet patches, 240 triangles, along
25–90% of those fronds. Triangle/cell intersection areas allocate the total:

$$q_i=q_{tree}\frac{A_{leaf,i}}{\sum_j A_{leaf,j}}.$$

The trunk and bare frond axes have zero assigned source area. A coarse cell may
contain both leaf and trunk; its CO₂ value is an air-cell average, not a trunk
photosynthesis rate. The patch count, fine geometry and approximately 0.733 m²
one-sided proxy area are schematic. They do not establish actual organ counts,
measured leaf area, stomatal behaviour or per-leaf photosynthesis.

`exchange.net_co2_umol_s` specifies the **whole-tree net** rate. Negative values
remove CO₂; positive values release it. The reference value,
−44.40047791098563 µmol/s, represents a February 2025 daytime Q95
uptake-strength case. Daytime was defined by measured global radiation
Rg ≥ 10 W/m², with no upper radiation limit. Q95 refers to uptake strength,
equivalent to Q05 of signed net exchange.

The example supplies that value as a constant input. It is not a daily mean,
a current-year observation or gross leaf photosynthesis. It does not adjust
to local CO₂, light or temperature. The original observations are not needed
to run the example. When supplying a different rate, also update its label
and record its measurement basis or model assumption.

Exchange continues while sealed, moving and fully open. No extra respiration
is added to a net rate. The solver uses reference air molar density P/(RT) for
conversion to ppm inventory. As an independent check, the first sealed interval
must satisfy

$$c_{ppm}(t)=c_{ppm}(0)+q_{tree}t\,\frac{RT}{PV},$$

when q_tree is in µmol/s; its micro-prefix cancels the ppm factor. Fans redistribute
concentration but cannot change this closed-volume mean balance.

## What the output means

The main curve is a volume-weighted mean over the original closed chamber
region, including while the shells move. Spatial standard deviation measures
variation between air cells. A smaller value means more uniform concentration.
It is not measurement uncertainty or a confidence interval. The source view
shows depth-integrated removal per projected area, with its own labelled scale.
All plots have white canvases and fixed linear viridis concentration scales.
The scale expands to include calculated extrema; it does not clip inconvenient
values. Frames come from saved fields at printed physical times; repeated frames
slow the display during wall movement without interpolating CFD fields.

## Verification and scope

Each run checks inventory conservation, moving-cell volumes, source allocation,
saved-field statistics and the first sealed interval's analytical mean.
Results are written to `verification/checks.json`.

The reference example passes these checks and its time-step sensitivity test.
Its first-closure grid test exceeds the specified 5% criterion: the relative
RMS change in spatial standard deviation is 5.83%. The
[verification report](VERIFICATION.md) gives the metrics and test scope.

The chamber's installed clear dimensions and fan performance remain unverified.
The model omits pole and fan-housing blockage, temperature and humidity
feedback, buoyancy, leakage through imperfect seals, and dynamic plant
physiology. Leaflet geometry and the spatial allocation of net uptake are
schematic. The transfer of a historical uptake rate to this nominal chamber
geometry is a diagnostic assumption. Opening-grid convergence has not been
established, so the example cannot establish an optimal fan arrangement or a
validated chamber mixing time.
