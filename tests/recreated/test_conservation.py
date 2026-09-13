"""Verification of the chamber's actual moving, three-dimensional mesh."""
import numpy as np
import pytest
from fluid_dynamic.recreated.mesh import Mesh,Study
from fluid_dynamic.recreated.operators import Projector,face_velocity,relative_flux,transport


def fluxes(m,g0,g1,dt):
    u=np.zeros((m.size,3));u[:,0]=m.cfg.wind
    qi,qb=face_velocity(m,u,(g1-g0)/dt)
    qi,qb,_,error=Projector().project(m,qi,qb,dt)
    ri,rb,gcl=relative_flux(m,qi,qb,g0,g1,dt)
    return ri,rb,error,gcl


@pytest.mark.parametrize('g0,g1',[(0,0),(0,.02),(.6,.62),(2,2),(.02,0)])
def test_uniform_field_and_geometric_conservation(g0,g1):
    cfg=Study(resolution=1,wind=.3);m=Mesh(cfg,(g0+g1)/2)
    ri,rb,div,gcl=fluxes(m,g0,g1,.2)
    old=np.full(m.size,400.)
    new,stats=transport(m,old,ri,rb,g0,g1,.2,cfg.diffusivity,
                        [400.]*6,correct=True,bounds=(400.,400.))
    assert np.max(abs(new-400.))<2e-6
    assert abs(stats['balance_residual'][0])<1e-6
    assert div<2e-7 and gcl<1e-8


def test_closed_panels_are_impermeable_and_interior_is_fluid():
    cfg=Study(resolution=1);m=Mesh(cfg,0.)
    x,y,z=m.xyz.reshape(-1,3).T
    inside=(abs(x)<2)&(abs(y)<2)&(z<6)
    assert np.isclose(sum(m.volume[inside]),96.)
    old=np.where(inside,340.,400.)
    ri,rb,_,_=fluxes(m,0.,0.,.2)
    new,stats=transport(m,old,ri,rb,0.,0.,.2,.1,[400.]*6,correct=True,bounds=(340.,400.))
    assert np.max(abs(new-old))<2e-7
    assert abs(stats['balance_residual'][0])<1e-6


@pytest.mark.parametrize('g0,g1',[(0,.1),(1.,1.1),(.1,0)])
def test_deficit_conservation_through_cell_birth_and_collapse(g0,g1):
    cfg=Study(resolution=1);m=Mesh(cfg,(g0+g1)/2)
    x,y,z=m.xyz.reshape(-1,3).T
    old=np.where((abs(x)<2.5)&(abs(y)<2)&(z<6),340.,400.)
    ri,rb,_,_=fluxes(m,g0,g1,.2)
    new,stats=transport(m,old,ri,rb,g0,g1,.2,.001,[400.]*6,correct=True,bounds=(340.,400.))
    assert new.min()>340.-2e-6 and new.max()<400.+2e-6
    assert abs(stats['balance_residual'][0])<1e-6


def test_closed_plant_source_changes_inventory_exactly():
    cfg=Study(resolution=1,wind=0);m=Mesh(cfg,0.)
    x,y,z=m.xyz.reshape(-1,3).T
    inside=(abs(x)<2)&(abs(y)<2)&(z<6)
    source=inside*m.volume*.7
    ri,rb,_,_=fluxes(m,0,0,.5)
    new,stats=transport(m,np.full(m.size,400.),ri,rb,0,0,.5,0,[400.]*6,source=source)
    assert np.isclose(np.sum((new-400)*m.volume),.5*source.sum(),atol=1e-7)
    assert abs(stats['balance_residual'][0])<1e-7
