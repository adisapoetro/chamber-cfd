"""Analytic transport oracles, independent of chamber output or target mixing."""
import numpy as np
from fluid_dynamic.recreated.operators import transport
from fluid_dynamic.recreated.verification import periodic_mesh


def advect_profile(resolution, profile, end=0.5, correct=True):
    m=periodic_mesh(resolution)
    x=m.xyz.reshape(-1,3)[:,0]
    c=profile(x)
    initial=c.copy()
    flux=[f['area'] if f['axis']==0 else np.zeros(len(f['left'])) for f in m.faces]
    dt=0.05/resolution
    steps=round(end/dt);dt=end/steps
    for _ in range(steps):
        c,stats=transport(m,c,flux,[],0,0,dt,0,[],correct=correct,
                          bounds=(float(initial.min()),float(initial.max())))
        assert abs(stats['balance_residual'][0])<1e-8
    return m,initial,c,profile(x-end)


def test_smooth_scalar_translation_converges_and_conserves_inventory(record_property):
    errors=[]
    profile=lambda x: 0.5+0.3*np.sin(2*np.pi*x/6)
    for resolution in (1,2,4):
        m,initial,c,exact=advect_profile(resolution,profile)
        errors.append(float(np.mean(abs(c-exact))))
        assert abs(np.dot(m.volume,c-initial))<1e-8
    record_property('scalar_translation_l1',str(errors))
    assert errors[0]>errors[1]>errors[2]
    assert errors[1]/errors[2]>1.5


def test_discontinuous_translation_is_bounded_without_mass_repair():
    profile=lambda x: ((x>-1)&(x<1)).astype(float)
    m,initial,c,_=advect_profile(4,profile,end=1.0)
    assert c.min()>=-1e-9 and c.max()<=1+1e-9
    assert abs(np.dot(m.volume,c-initial))<1e-8
    line=c.reshape(m.shape)[:,0,0]
    variation=np.sum(abs(line-np.roll(line,1)))
    assert variation<=2+1e-8
