"""Independent analytic checks of the viscous momentum operator."""
import numpy as np
from fluid_dynamic.recreated.mesh import Study,Mesh
from fluid_dynamic.recreated.operators import transport,Projector


def test_poiseuille_viscous_solution_converges_under_refinement():
    errors=[]
    for resolution in (1,2,4):
        m=Mesh(Study(resolution=resolution,width=2.,depth=2.,height=2.,outer_x=6.,outer_y=6.,outer_z=4.),0)
        for f in m.faces:f['wall'][:]=False
        for f in m.boundaries:f['kind']='wall'
        z=m.xyz.reshape(-1,3)[:,2];height=4.;nu=.1
        exact=z*(height-z)
        bc=[]
        for f in m.boundaries:
            zz=f['position'][:,2]
            values=np.zeros((len(zz),3));values[:,0]=zz*(height-zz);bc.append(values)
        source=np.zeros((m.size,3));source[:,0]=2*nu*m.volume
        result,_=transport(m,np.zeros((m.size,3)),[np.zeros(len(f['left'])) for f in m.faces],
            [np.zeros(len(f['cell'])) for f in m.boundaries],0,0,1e7,nu,bc,source=source,
            wall_values=[np.zeros((len(f['left']),3)) for f in m.faces])
        errors.append(float(np.sqrt(np.mean((result[:,0]-exact)**2))))
    assert errors[0]/errors[1]>3.5
    assert errors[1]/errors[2]>3.5


def test_pressure_projection_annuls_an_irrotational_flux():
    m=Mesh(Study(resolution=1,wind=0),1.)
    p=np.sin(m.xyz.reshape(-1,3)[:,0])+.2*m.xyz.reshape(-1,3)[:,2]**2
    qi=[];qb=[]
    for f in m.faces:
        qi.append(f['area']*(p[f['right']]-p[f['left']])/f['distance']*(~f['wall']))
    for f in m.boundaries:
        qb.append(-p[f['cell']]*f['area']/f['distance'] if f['kind']=='open' else np.zeros(len(f['cell'])))
    qi,qb,_,residual=Projector().project(m,qi,qb,1.)
    assert max(np.max(abs(q),initial=0) for q in qi+qb)<1e-8
    assert residual<1e-8
