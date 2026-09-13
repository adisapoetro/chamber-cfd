"""Analytic transient Taylor–Green vortex check of the coupled operators.

The analytic velocity is the oracle; no stored chamber output is used as truth.
This periodic, smooth laminar check cannot validate chamber turbulence or walls.
"""
import numpy as np
from .mesh import Mesh,Study
from .operators import Projector,face_velocity,interpolate,transport
from .solver import gradients,cross_stress


def periodic_mesh(resolution):
    m=Mesh(Study(resolution=resolution,width=2.,depth=2.,height=2.,outer_x=6.,outer_y=6.,outer_z=6.,wind=0),0.)
    ids=np.arange(m.size).reshape(m.shape)
    for a,f in enumerate(m.faces):
        low=[slice(None)]*3;high=low.copy();low[a]=0;high[a]=-1
        l=ids[tuple(high)].ravel();r=ids[tuple(low)].ravel()
        position=m.xyz[tuple(high)].reshape(-1,3).copy();position[:,a]=m.nodes[a][-1]
        count=len(l)
        for key,extra in [('left',l),('right',r),('area',np.full(count,1/resolution**2)),
                          ('distance',np.full(count,1/resolution)),('position',position),
                          ('wall',np.zeros(count,bool)),('normal',np.ones(count))]:
            f[key]=np.concatenate([f[key],extra])
        f['wall'][:]=False
    m.boundaries=[]
    return m


def taylor_green(resolution,dt_scale=.05,end=.5):
    m=periodic_mesh(resolution);nu=.05;k=2*np.pi/6
    x,y,_=m.xyz.reshape(-1,3).T
    initial=np.stack([np.sin(k*x)*np.cos(k*y),-np.cos(k*x)*np.sin(k*y),np.zeros_like(x)],axis=-1)
    u=initial.copy();qi,qb=face_velocity(m,u,0);projector=Projector()
    dt=dt_scale/resolution;steps=round(end/dt);dt=end/steps
    qi,qb,_,_=projector.project(m,qi,qb,dt)
    for _ in range(steps):
        # Same predictor, symmetric stress, native-face update and projection as
        # the chamber run; no canopy, turbulence or moving geometry in this oracle.
        grad=gradients(m,u,0.)
        star,_=transport(m,u,qi,[],0.,0.,dt,nu,[],source=cross_stress(m,np.full(m.size,nu),grad),
                        wall_values=[np.zeros((len(f['left']),3)) for f in m.faces])
        delta=star-u
        predicted=[q+interpolate(m,delta,f)[:,f['axis']]*f['area'] for f,q in zip(m.faces,qi)]
        qi,_,pgrad,div=projector.project(m,predicted,[],dt)
        u=star-dt*pgrad
    exact=initial*np.exp(-2*nu*k*k*end)
    relative_l2=float(np.sqrt(np.sum((u-exact)**2)/np.sum(exact**2)))
    energy=float(np.sum(u*u)/np.sum(exact*exact))
    return dict(resolution=resolution,cells=m.size,dt_s=dt,end_s=end,relative_velocity_l2=relative_l2,
                energy_ratio_to_exact=energy,max_final_divergence_s_inv=div)


def verification_suite():
    spatial=[taylor_green(r) for r in [1,2,4,8]]
    time_refined=taylor_green(8,.025)
    errors=[s['relative_velocity_l2'] for s in spatial]
    orders=[float(np.log(errors[i]/errors[i+1])/np.log(2)) for i in range(len(errors)-1)]
    passed=all(p>.5 for p in orders) and errors[-1]<.05 and time_refined['relative_velocity_l2']<.05
    return dict(status='PASS' if passed else 'FAIL',spatial=spatial,observed_orders=orders,
                time_refined=time_refined,criterion='Velocity L2 <5% on finest grid and observed refinement order >0.5; fine half timestep also <5%.',
                limitation='Smooth periodic laminar flow only; this does not qualify wall shear, turbulence or chamber flushing.')
