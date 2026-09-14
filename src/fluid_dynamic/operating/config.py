"""Public, unit-explicit inputs. No research table, release or older run is read."""
from dataclasses import asdict
from pathlib import Path
from typing import Literal
import math
import numpy as np
import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator
from ..scenarios.geometry import ScenarioStudy
from ..recreated.fan_forcing import FanArray


class Strict(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True, allow_inf_nan=False, validate_default=True)


class Chamber(Strict):
    width_m: float = Field(4., gt=0)
    depth_m: float = Field(6., gt=0)
    height_m: float = Field(4.0, gt=0)
    floor_m: float = Field(.6, ge=0)
    maximum_gap_m: float = Field(2., gt=0)


class Domain(Strict):
    # Null sizes retain the original exterior clearances when the chamber changes.
    x_m: float | None = Field(None, gt=0)
    y_m: float | None = Field(None, gt=0)
    z_m: float | None = Field(None, gt=0)


class Air(Strict):
    wind_speed_m_s: float = Field(1.5, ge=0)
    wind_to_deg: float = Field(0., ge=0, lt=360)
    ambient_co2_ppm: float = Field(400., gt=0)
    initial_co2_ppm: float = Field(400., gt=0)
    temperature_k: float = Field(298.15, gt=0)
    pressure_pa: float = Field(101325., gt=0)


class Cycle(Strict):
    closed_s: float = Field(300., gt=0)
    open_phase_s: float = Field(600., gt=0)
    opening_travel_s: float = Field(20., gt=0)
    closing_travel_s: float = Field(20., gt=0)
    cycles: int = Field(2, ge=1)

    @model_validator(mode='after')
    def travel_fits(self):
        if self.open_phase_s < self.opening_travel_s + self.closing_travel_s:
            raise ValueError('open_phase_s includes both travel intervals and must contain them')
        return self


class Palm(Strict):
    height_m: float = Field(2.39, gt=0)
    crown_radius_x_m: float = Field(1.52, gt=0)
    crown_radius_y_m: float = Field(1.52, gt=0)
    trunk_radius_m: float = Field(.14483099821362477, gt=0)


class Exchange(Strict):
    net_co2_umol_s: float = -44.40047791098563
    label: str = Field('Prescribed February 2025 daytime Q95 net uptake', min_length=1, max_length=100)
    # Whole-tree net rate, negative uptake / positive release; constant in every phase.
    model: Literal['constant_leaflet_net_exchange'] = 'constant_leaflet_net_exchange'


class Fans(Strict):
    count: int = Field(3, ge=0, le=100)
    pole_x_m: float = .3
    pole_y_m: float = 0.
    heights_m: list[float] | None = None  # null: centres of count equal height bands
    azimuths_deg: list[float] | None = None  # null: 360/count, starting at first_azimuth_deg
    first_azimuth_deg: float = Field(0., ge=0, lt=360)
    model: str = Field('Delta AFC1212DE', min_length=1)
    diameter_m: float = Field(.12, gt=0)
    free_air_cfm_per_fan: float = Field(148.3, gt=0)
    momentum_factor: float = Field(1., ge=0)
    kernel_sigma_m: float = Field(.3, gt=0)
    schedule: Literal['closed_only'] = 'closed_only'

    @model_validator(mode='after')
    def lengths_match(self):
        for name in ('heights_m', 'azimuths_deg'):
            values = getattr(self, name)
            if values is not None and len(values) != self.count:
                raise ValueError(f'{name} must contain exactly count entries (or use null)')
        if self.azimuths_deg and not all(0 <= a < 360 for a in self.azimuths_deg):
            raise ValueError('Fan azimuths use flow TO in [0, 360) degrees')
        return self


class Numerics(Strict):
    cell_m: float = Field(.5, gt=0)
    dt_s: float = Field(.5, gt=0)
    save_every_s: float = Field(10., gt=0)
    duration_s: float | None = Field(None, gt=0)  # null: complete configured cycles
    viscosity_m2_s: float = Field(1.5e-5, gt=0)
    molecular_diffusivity_m2_s: float = Field(1.6e-5, ge=0)
    smagorinsky: float = Field(.12, ge=0)
    turbulent_schmidt: float = Field(.7, gt=0)
    canopy_drag_per_m: float = Field(.3, ge=0)
    trunk_drag_per_m: float = Field(20., ge=0)


class Screens(Strict):
    enabled: bool = True
    duration_s: float | None = Field(None, gt=0)  # null: first closed phase
    grid_refinement: float = Field(1.5, gt=1)
    timestep_refinement: float = Field(2., gt=1)
    relative_tolerance: float = Field(.05, gt=0)


class Display(Strict):
    co2_min_ppm: float = Field(380., ge=0)
    co2_max_ppm: float = Field(405., gt=0)
    fps: int = Field(6, ge=1, le=60)

    @model_validator(mode='after')
    def ordered(self):
        if self.co2_min_ppm >= self.co2_max_ppm:
            raise ValueError('CO2 colour minimum must be below maximum')
        return self


class Inputs(Strict):
    schema_version: Literal[1] = 1
    title: str = Field('Operating chamber · central fans and leaf uptake', min_length=1, max_length=80)
    chamber: Chamber = Field(default_factory=Chamber)
    domain: Domain = Field(default_factory=Domain)
    air: Air = Field(default_factory=Air)
    cycle: Cycle = Field(default_factory=Cycle)
    palm: Palm = Field(default_factory=Palm)
    exchange: Exchange = Field(default_factory=Exchange)
    fans: Fans = Field(default_factory=Fans)
    numerics: Numerics = Field(default_factory=Numerics)
    screens: Screens = Field(default_factory=Screens)
    display: Display = Field(default_factory=Display)


def load(path, overrides=()):
    raw = yaml.safe_load(Path(path).read_text())
    if not isinstance(raw, dict):
        raise ValueError('Input must be a YAML mapping')
    # Resolve defaults before applying typed dotted overrides. Unknown keys still fail.
    raw = Inputs.model_validate(raw).model_dump()
    for item in overrides:
        key, sep, value = item.partition('=')
        if not sep:
            raise ValueError('Override syntax: group.variable=value')
        node = raw
        parts = key.split('.')
        for part in parts[:-1]:
            if part not in node or not isinstance(node[part], dict):
                raise ValueError('Unknown input: '+key)
            node = node[part]
        if parts[-1] not in node:
            raise ValueError('Unknown input: '+key)
        node[parts[-1]] = yaml.safe_load(value)
    return Inputs.model_validate(raw)


def compile_study(p):
    c, d, a, cy, n, f = p.chamber, p.domain, p.air, p.cycle, p.numerics, p.fans
    period = cy.closed_s + cy.open_phase_s
    duration = n.duration_s if n.duration_s is not None else period*cy.cycles
    if duration > period*cy.cycles:
        raise ValueError('duration_s cannot exceed cycle.cycles complete periods')
    plant = dict(height_m=p.palm.height_m,
        crown_radii_m=[p.palm.crown_radius_x_m, p.palm.crown_radius_y_m],
        trunk_radius_m=p.palm.trunk_radius_m, net_exchange_umol_s=p.exchange.net_co2_umol_s)
    if (plant['height_m'] >= c.height_m or plant['crown_radii_m'][0] >= c.width_m/2
        or plant['crown_radii_m'][1] >= c.depth_m/2):
        raise ValueError('Palm envelope must fit inside the closed chamber; change palm dimensions explicitly')
    cfg = ScenarioStudy(width=c.width_m, depth=c.depth_m, height=c.height_m, gap=c.maximum_gap_m,
        geometry=dict(kind='sliding_shells', floor_z_m=c.floor_m), plant=plant,
        outer_x=d.x_m or c.width_m+c.maximum_gap_m+8., outer_y=d.y_m or c.depth_m+4.,
        outer_z=d.z_m or c.floor_m+c.height_m+4.7,
        wind=a.wind_speed_m_s, wind_angle_deg=a.wind_to_deg, inside_ppm=a.initial_co2_ppm,
        ambient_ppm=a.ambient_co2_ppm, temperature_k=a.temperature_k, pressure_pa=a.pressure_pa,
        open_at=cy.closed_s, opening_s=cy.opening_travel_s,
        close_at=period-cy.closing_travel_s, closing_s=cy.closing_travel_s, period_s=period,
        duration_s=duration, dt_s=n.dt_s, save_every_s=n.save_every_s, target_cell_m=n.cell_m,
        viscosity=n.viscosity_m2_s, diffusivity=n.molecular_diffusivity_m2_s,
        smagorinsky=n.smagorinsky, turbulent_schmidt=n.turbulent_schmidt,
        canopy_drag_per_m=n.canopy_drag_per_m, trunk_drag_per_m=n.trunk_drag_per_m,
        tree_source_umol_s=p.exchange.net_co2_umol_s)
    heights = f.heights_m if f.heights_m is not None else [(c.height_m/(2*f.count))*(2*i+1) for i in range(f.count)]
    angles = f.azimuths_deg if f.azimuths_deg is not None else [(f.first_azimuth_deg+360*i/f.count)%360 for i in range(f.count)]
    # Stable special-angle vectors preserve the accepted 120-degree baseline.
    special = {0.: [1.,0.,0.], 120.: [-.5,math.sqrt(3)/2,0.], 240.: [-.5,-math.sqrt(3)/2,0.]}
    axes = [special.get(v, [math.cos(math.radians(v)),math.sin(math.radians(v)),0.]) for v in angles]
    fans = FanArray(model=f.model, heights_above_floor_m=heights, axes=axes if f.count else None,
        closed_x_m=f.pole_x_m, closed_y_m=f.pole_y_m, mounting_shell=0,
        diameter_m=f.diameter_m, nominal_free_air_cfm=f.free_air_cfm_per_fan,
        momentum_factor=f.momentum_factor, kernel_sigma_m=f.kernel_sigma_m,
        enabled=f.count>0, schedule=f.schedule)
    cfg.validate(); fans.validate(cfg)
    cases = {'main': dict(config=asdict(cfg), fans=asdict(fans), role='main', source_geometry='geometry/palm.npz')}
    if p.screens.enabled:
        end = p.screens.duration_s if p.screens.duration_s is not None else min(cy.closed_s, duration)
        if end > duration or end < min(n.save_every_s, duration):
            raise ValueError('Screen duration must include a saved sample and not exceed the main run')
        for name, key, value in [('grid','target_cell_m',n.cell_m/p.screens.grid_refinement),
                                 ('dt','dt_s',n.dt_s/p.screens.timestep_refinement)]:
            case = dict(config={**asdict(cfg), 'duration_s': end, key: value}, fans=asdict(fans),
                role='numerical_screen', source_geometry='geometry/palm.npz')
            cases[name] = case
    geometry = dict(pairs_per_frond=12, bare_base_fraction=.25, last_leaflet_fraction=.9,
        minimum_leaflet_length_m=.08, length_amplitude_m=.22, leaflet_half_width_m=.025, tip_forward_m=.035)
    return dict(workflow_version=1, title=p.title, inputs=p.model_dump(), cases=cases, geometry=geometry,
        display=dict(preferred_limits_ppm=[p.display.co2_min_ppm,p.display.co2_max_ppm], fps=p.display.fps))
