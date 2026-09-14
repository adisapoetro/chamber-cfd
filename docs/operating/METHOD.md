# How the operating-chamber simulation works

This is the owner-selected **4 × 6 × 4 m nominal operating chamber**, with two sliding
halves, an elevated floor and a palm inside. The current configuration has three
fans on a fixed pole approximately 0.30 m beside the trunk. They sit at
0.667, 2.000 and 3.333 m above the floor, blowing toward 0°, 120° and 240°.
The support location and absolute direction are prescribed approximations.
The design documents call the plan dimensions 6 m width × 4 m depth; the
solver's x/y names are reversed because x follows sliding travel. The resulting
rectangular air volume is 96 m³. The owner adopted the documented nominal
4 m clear height on 14 September 2026; installed inside dimensions and net
air volume remain unverified. The floor datum remains the assumed 0.6 m.

The Python runtime is `fluid_dynamic.operating` →
`fluid_dynamic.recreated.fan_solver`, using NumPy/SciPy finite-volume operators
and `fluid_dynamic.scenarios.geometry`. The repository also contains historical
PhiFlow code; **the selected central-fan bundle and this reproduction do not use
PhiFlow or WaterLily**. No plant model or notebook is executed to obtain inputs.

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

The main illustration retains the accepted simple trunk and five fronds.
The source uses 120 representative leaflet patches, 240 triangles, along
25–90% of those fronds. Triangle/cell intersection areas allocate the total:

$$q_i=q_{tree}\frac{A_{leaf,i}}{\sum_j A_{leaf,j}}.$$

The trunk and bare frond axes have zero assigned source area. A coarse cell may
contain both leaf and trunk; its CO₂ value is an air-cell average, not a trunk
photosynthesis rate. The patch count, fine geometry and approximately 0.733 m²
one-sided proxy area are schematic. They do not establish actual organ counts,
measured leaf area, stomatal behaviour or per-leaf photosynthesis.

`exchange.net_co2_umol_s` is a prescribed **whole-tree net** rate. Negative
removes CO₂; positive releases it. The default −44.40047791098563 µmol/s retains
the owner-selected February 2025 daytime Q95 uptake-strength diagnostic. It is
not a 2026 observation, a typical daily mean or gross leaf photosynthesis. The
local observation provenance remains with the preserved study and is not needed
or exported by this workflow. Changing the rate does not rerun a statistical
selection. Update its label and provenance when supplying another value.

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

Each completed run records inventory, moving-volume, source-area, saved-field
and sealed-mean checks in `verification/checks.json`. Read the same file for
that run's first-closure grid and time-step sensitivity results. The previous
3.7 m run failed the grid screen (about 6.87% change in spatial standard
deviation against 5%) and passed the time-step screen (about 0.56%); those
numbers do not describe the revised 4 m run. Opening-grid convergence remains
unestablished. These engineering
checks do not validate installed fan performance, leaf physiology or actual
chamber mixing time. The original uptake-to-geometry transfer and plant identity
remain conditional. Selecting this as the current design does not remove those
limits. No calibration, qualified forecast or field validation is claimed.
