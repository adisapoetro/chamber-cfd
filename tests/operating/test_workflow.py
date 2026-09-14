"""Portable physical/configuration checks; no prior output or research data needed."""
from pathlib import Path
import json
import numpy as np
import pytest
from pydantic import ValidationError
from fluid_dynamic.operating.config import Inputs, compile_study, load
from fluid_dynamic.operating.workflow import prepare, run, verify
from fluid_dynamic.operating.audit import audit
from fluid_dynamic.recreated.fan_forcing import FanArray, FanForcing
from fluid_dynamic.scenarios.geometry import ScenarioStudy, ScenarioMesh
from fluid_dynamic.recreated.mesh import opening

ROOT=Path(__file__).resolve().parents[2]


def test_operating_schedule_and_three_fans():
    study=compile_study(Inputs()); case=study['cases']['main']
    cfg=ScenarioStudy(**case['config']); f=FanArray(**case['fans'])
    assert cfg.width*cfg.depth*cfg.height == pytest.approx(96.0)
    np.testing.assert_allclose(f.heights_above_floor_m, [2/3, 2, 10/3])
    assert ScenarioMesh(cfg,0).reporting_weights('fixed').sum() == pytest.approx(96.0)
    assert cfg.duration_s==1800
    for t,gap,on in [(0,0,True),(300,0,False),(310,1,False),(320,2,False),
                     (880,2,False),(890,1,False),(900,0,True),(1800,0,True)]:
        assert opening(cfg,t)==pytest.approx(gap)
        assert f.active(cfg,t)==on
    np.testing.assert_allclose(f.directions()@f.directions().T,np.eye(3)*1.5-.5,atol=1e-14)
    np.testing.assert_array_equal(f.positions(cfg,0),f.positions(cfg,2))
    assert study['cases']['grid']['config']['target_cell_m']==1/3
    assert study['cases']['dt']['config']['dt_s']==.25


@pytest.mark.parametrize('override',[
    'fans.count=-1','fans.count=1.5','air.wind_speed_m_s=-1','air.wind_to_deg=360',
    'chamber.height_m=2','chamber.width_m=3','fans.heights_m=[1]',
    'fans.azimuths_deg=[0,120,400]','fans.kernel_sigma_m=0.01',
    'cycle.open_phase_s=10','exchange.net_co2_umol_s=.nan','numerics.cell_m=0',
    'air.wind_speed=2','screens.duration_s=2000','fans.schedule=always',
    'domain.x_m=5','exchange.model=xpalm','fans.count=true'])
def test_invalid_physical_inputs_fail(override):
    with pytest.raises((ValueError,ValidationError)):
        compile_study(load(ROOT/'configs/operating/current.yaml',[override]))


@pytest.mark.parametrize('count',[0,1,4])
def test_force_count_direction_and_volume_conservation(count):
    p=load(ROOT/'configs/operating/current.yaml',[f'fans.count={count}',
        'chamber.height_m=4.2','chamber.width_m=5','air.wind_to_deg=90'])
    case=compile_study(p)['cases']['main']; cfg=ScenarioStudy(**case['config']); f=FanArray(**case['fans'])
    mesh=ScenarioMesh(cfg,0); source=FanForcing(cfg,f).source(mesh,0)
    scale=f.properties(cfg)['force_n_per_fan']/f.properties(cfg)['air_density_kg_m3']
    np.testing.assert_allclose(source.sum(axis=0),scale*f.directions().sum(axis=0),atol=3e-14)
    np.testing.assert_allclose(cfg.wind_vector,[0,1.5,0],atol=1e-14)
    assert len(f.positions(cfg,0))==count
    np.testing.assert_array_equal(source[mesh.reporting_weights('fixed')==0],0)
    if count==0: np.testing.assert_array_equal(source,0)
    else: assert np.linalg.norm(source)>0


@pytest.mark.parametrize('count,rate',[(0,0),(1,-25),(4,10)])
def test_changed_inputs_through_opening_conserve_leaf_exchange(tmp_path,count,rate):
    p=load(ROOT/'configs/operating/current.yaml',[
        'chamber.width_m=5','chamber.depth_m=5','chamber.height_m=4.2',
        'palm.height_m=2.6','palm.crown_radius_x_m=1.6','air.wind_to_deg=90',
        'air.wind_speed_m_s=0.8',f'fans.count={count}',f'exchange.net_co2_umol_s={rate}',
        'cycle.closed_s=2','cycle.open_phase_s=2','cycle.opening_travel_s=0.5',
        'cycle.closing_travel_s=0.5','cycle.cycles=1','numerics.cell_m=1',
        'numerics.dt_s=0.25','numerics.save_every_s=0.5','screens.enabled=false'])
    out=tmp_path/'run'; study=prepare(ROOT,out,p); run(out,study,1)
    result=audit(out,ROOT,media=False)
    assert result['cases']['main']['volume_m3']==105
    assert result['cases']['main']['max_budget_fraction']<1e-6
    assert result['cases']['main']['sealed_error_ppm']<1e-6
    assert result['saved_fields']==9
    if rate==0:
        for frame in json.loads((out/'data/cfd_runs/main/summary.json').read_text())['frames']:
            with np.load(out/'data/cfd_runs/main'/frame['path']) as d:
                np.testing.assert_allclose(d['co2'],400,atol=1e-8,rtol=0)
    with pytest.raises(ValueError,match='already exists'):
        prepare(ROOT,out,p)
    (out/'configs/inputs.json').write_text('{}')
    with pytest.raises(ValueError,match='Frozen input changed'):
        verify(out,ROOT)
