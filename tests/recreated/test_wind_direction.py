"""Boundary, inventory and reflection oracles for chamber-relative wind."""
from dataclasses import replace
import json
import numpy as np
import pytest
from fluid_dynamic.recreated.mesh import Mesh, Study
from fluid_dynamic.recreated.operators import Projector, face_velocity, relative_flux, transport
from fluid_dynamic.recreated.solver import run_case


@pytest.mark.parametrize('angle,expected',[(0,(1,0,0)),(45,(2**-.5,2**-.5,0)),
    (90,(0,1,0)),(135,(-2**-.5,2**-.5,0)),(180,(-1,0,0)),
    (225,(-2**-.5,-2**-.5,0)),(270,(0,-1,0)),(315,(2**-.5,-2**-.5,0))])
def test_flow_to_vector_and_upstream_faces(angle,expected):
    cfg=Study(resolution=1,wind=1.5,wind_angle_deg=angle)
    cfg.validate();np.testing.assert_allclose(cfg.wind_vector,1.5*np.array(expected),atol=1e-14)
    m=Mesh(cfg,2.)
    for f in m.boundaries:
        normal_speed=cfg.wind_vector[f['axis']]*f['side']
        assert f['kind']==('wall' if (f['axis'],f['side'])==(2,-1) else
                           'inlet' if normal_speed<0 else 'open')
        if f['kind']=='inlet':
            np.testing.assert_allclose(m.boundary_wall_velocity(f,0),
                np.broadcast_to(cfg.wind_vector,(len(f['cell']),3)),atol=0)
    # The pressure operator must be rebuilt when a different side is prescribed.
    reused=Projector();reused.matrix(Mesh(replace(cfg,wind_angle_deg=(angle+90)%360),2.))
    a,_=reused.matrix(m);b,_=Projector().matrix(m)
    assert (a!=b).nnz==0


@pytest.mark.parametrize('angle',[-1,360,float('nan')])
def test_invalid_direction_rejected(angle):
    with pytest.raises(ValueError):Study(wind_angle_deg=angle).validate()


def test_calm_has_no_artificial_fixed_horizontal_inlet():
    m=Mesh(Study(resolution=1,wind=0,wind_angle_deg=135),0)
    assert all(f['kind']=='open' for f in m.boundaries if f['axis']<2)
    assert np.array_equal(m.cfg.wind_vector,np.zeros(3))


@pytest.mark.parametrize('angle',[0,45,90,135,180,225,270,315])
@pytest.mark.parametrize('g0,g1',[(0,0),(0,.02),(.02,0)])
def test_rotated_wind_preserves_uniform_gas_and_sealed_inventory(angle,g0,g1):
    cfg=Study(resolution=1,wind_angle_deg=angle);m=Mesh(cfg,(g0+g1)/2)
    u=np.broadcast_to(cfg.wind_vector,(m.size,3)).copy()
    qi,qb=face_velocity(m,u,(g1-g0)/.2)
    qi,qb,_,div=Projector().project(m,qi,qb,.2)
    qi,qb,gcl=relative_flux(m,qi,qb,g0,g1,.2)
    x,y,z=m.xyz.reshape(-1,3).T
    old=np.where((abs(x)<2)&(abs(y)<2)&(z<6),340.,400.) if g0==g1 else np.full(m.size,400.)
    new,stats=transport(m,old,qi,qb,g0,g1,.2,.1,[400.]*6,correct=True,bounds=(340.,400.))
    np.testing.assert_allclose(new,old,atol=2e-6,rtol=0)
    assert abs(stats['balance_residual'][0])<1e-6
    assert div<2e-7 and gcl<1e-8


@pytest.mark.parametrize('a,b,axis',[(0,180,0),(45,135,0),(45,315,1),(90,270,1)])
def test_coupled_opening_reflects_wind_and_solution(tmp_path,a,b,axis):
    cfg=Study(resolution=1,outer_y=14,open_at=1,opening_s=2,duration_s=6,dt_s=.2,save_every_s=1)
    for angle in (a,b):run_case(replace(cfg,wind_angle_deg=angle),tmp_path/str(angle),progress=False)
    with np.load(tmp_path/str(a)/'snapshots/0006.npz') as da, np.load(tmp_path/str(b)/'snapshots/0006.npz') as db:
        np.testing.assert_allclose(da['co2'],np.flip(db['co2'],axis=axis),atol=2e-5,rtol=0)
        reflected=np.flip(db['velocity'],axis=axis).copy();reflected[...,axis]*=-1
        np.testing.assert_allclose(da['velocity'],reflected,atol=2e-5,rtol=0)
    for angle in (a,b):
        summary=json.loads((tmp_path/str(angle)/'summary.json').read_text())
        assert summary['max_budget_fraction']<.001
        assert summary['max_gcl_m3']<1e-8
