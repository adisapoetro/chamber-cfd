"""Compatible pressure projection and conservative ALE inventory transport."""
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import LinearOperator, cg, splu, bicgstab, gmres


def interpolate(mesh,values,f):
    l,r=f['left'],f['right'];a=f['axis']
    x=mesh.xyz.reshape(-1,3)[:,a];d=f['distance']
    weight=(f['position'][:,a]-x[l])/d
    return values[l]*(1-weight[:,None])+values[r]*weight[:,None] if values.ndim==2 else values[l]*(1-weight)+values[r]*weight


def face_velocity(mesh,u,rate):
    qi=[];qb=[]
    for f in mesh.faces:
        v=interpolate(mesh,u,f)[:,f['axis']]
        v[f['wall']]=mesh.wall_velocity(f['position'][f['wall']],rate)[:,f['axis']]
        qi.append(v*f['area'])
    for f in mesh.boundaries:
        v=u[f['cell'],f['axis']] if f['kind']=='open' else mesh.boundary_wall_velocity(f,rate)[:,f['axis']]
        qb.append(v*f['side']*f['area'])
    return qi,qb


class Projector:
    """Pressure correction on open faces only, with component-wise gauges.

    An earlier factorization is used solely as a preconditioner when the mesh
    moves. CG still solves the CURRENT matrix to the stated residual tolerance.
    """
    def __init__(self):
        self.signature=None;self.preconditioner=None;self.last_gap=None
        self.last_pressure=None;self.factorizations=0;self.refactorizations_after_cg_failure=0
        self.matrix_key=None;self.cached_matrix=None

    def matrix(self,m):
        key=(m.shape,m.gap,tuple(f['kind'] for f in m.boundaries))
        if key==self.matrix_key:return self.cached_matrix
        rows=[];cols=[];vals=[];diag=np.zeros(m.size)
        for f in m.faces:
            k=f['area']/f['distance']*(~f['wall']);l,r=f['left'],f['right']
            np.add.at(diag,l,k);np.add.at(diag,r,k)
            rows.extend([l,r]);cols.extend([r,l]);vals.extend([-k,-k])
        outside=np.zeros(m.size,bool)
        for f in m.boundaries:
            if f['kind']=='open':
                k=f['area']/f['distance'];np.add.at(diag,f['cell'],k);outside[f['cell']]=True
        rows.append(np.arange(m.size));cols.append(np.arange(m.size));vals.append(diag)
        mat=coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(m.size,m.size)).tocsr()
        mat.eliminate_zeros()
        _,labels=connected_components(mat,directed=False)
        anchored=np.unique(labels[outside]);gauges=[]
        for label in np.unique(labels):
            if label not in anchored:gauges.append(np.flatnonzero(labels==label)[0])
        if gauges:
            mat=mat.tolil()
            for g in gauges:mat[g,:]=0;mat[:,g]=0;mat[g,g]=max(diag[g],1.)
            mat=mat.tocsr()
        self.matrix_key=key;self.cached_matrix=(mat,gauges)
        return self.cached_matrix

    def project(self,m,qi,qb,dt):
        mat,gauges=self.matrix(m)
        rhs=-m.net_outflow(qi,qb)/dt
        # Never repair incompatible enclosed-volume forcing by subtracting its mean.
        rhs[gauges]=0.
        signature=(m.shape,tuple(gauges),tuple(f['kind'] for f in m.boundaries))
        ratio=(max(m.gap,1e-12)/max(self.last_gap or 0.,1e-12))
        stationary_endpoint=m.gap in (0.,m.cfg.gap) and self.last_gap!=m.gap
        if signature!=self.signature or not .25<ratio<4 or stationary_endpoint:
            factor=splu(mat.tocsc())
            self.preconditioner=LinearOperator(mat.shape,matvec=factor.solve)
            self.signature=signature;self.last_gap=m.gap;self.factorizations+=1
            self.last_pressure=None
        p,info=cg(mat,rhs,x0=self.last_pressure,M=self.preconditioner,rtol=2e-11,atol=1e-12,maxiter=300)
        if info:
            # Crossing fixed/moving planes can invalidate a previously useful
            # preconditioner without changing matrix size or pressure gauges.
            # Refactor the CURRENT matrix; retain the identical RHS/tolerances.
            factor=splu(mat.tocsc())
            self.preconditioner=LinearOperator(mat.shape,matvec=factor.solve)
            self.signature=signature;self.last_gap=m.gap;self.factorizations+=1
            self.refactorizations_after_cg_failure+=1
            p,info=cg(mat,rhs,M=self.preconditioner,rtol=2e-11,atol=1e-12,maxiter=300)
        if info:raise RuntimeError(f'Pressure solve did not converge after current-matrix refactorization: {info}')
        self.last_pressure=p
        outi=[];outb=[];grad=np.zeros((m.size,3));weight=np.zeros_like(grad)
        for f,q in zip(m.faces,qi):
            l,r=f['left'],f['right'];axis=f['axis']
            g=(p[r]-p[l])/f['distance']*(~f['wall'])
            outi.append(q-dt*f['area']*g)
            # Cell correction reconstructed from the same corrected face gradients.
            for cells in (l,r):
                np.add.at(grad[:,axis],cells,g*f['area'])
                np.add.at(weight[:,axis],cells,f['area'])
        for f,q in zip(m.boundaries,qb):
            if f['kind']=='open':
                g=-p[f['cell']]/f['distance']
                outb.append(q-dt*f['area']*g)
                np.add.at(grad[:,f['axis']],f['cell'],g*f['side']*f['area'])
                np.add.at(weight[:,f['axis']],f['cell'],f['area'])
            else:outb.append(q.copy())
        grad/=np.maximum(weight,1e-30)
        divergence=m.net_outflow(outi,outb)
        residual=float(np.max(abs(divergence)/m.volume))
        if residual>2e-7:raise RuntimeError(f'Incompressibility failed: {residual:.3g} /s')
        return outi,outb,grad,residual


def relative_flux(m,qi,qb,g0,g1,dt):
    wi,wb=m.sweep(g0,g1,dt)
    ri=[q-w for q,w in zip(qi,wi)];rb=[q-w for q,w in zip(qb,wb)]
    for f,q in zip(m.faces,ri):
        if np.max(abs(q[f['wall']]),initial=0)>1e-10:raise RuntimeError('Moving wall leaks')
        q[f['wall']]=0.
    gcl=m.volumes_at(g1)-m.volumes_at(g0)+dt*m.net_outflow(ri,rb)
    error=float(np.max(abs(gcl)))
    if error>1e-8:raise RuntimeError(f'Geometric conservation failed: {error}')
    return ri,rb,error


def transport(m,old,qi,qb,g0,g1,dt,diffusivity,boundary_values,
              source=None,wall_values=None,reaction=None,correct=False,bounds=None,
              local_bounds=True):
    """Implicit donor-cell finite volume with optional conservative limited correction.

    Each internal face adds equal and opposite inventory. No clipping, ambient
    fill, mixing relaxation, or total-mass renormalization is performed.
    The sparse solve includes zero-capacity newborn/collapsing gap cells.
    """
    vector=old.ndim==2
    previous=old if vector else old[:,None]
    ncomp=previous.shape[1]
    v0=m.volumes_at(g0);v1=m.volumes_at(g1)
    diag=v1.copy();rhs=v0[:,None]*previous
    if source is not None:rhs+=dt*(source if vector else source[:,None])
    if reaction is not None:diag+=dt*v1*reaction
    rows=[];cols=[];vals=[];conductances=[];wall_exchange=np.zeros(ncomp)
    diffusivity=np.broadcast_to(np.asarray(diffusivity,dtype=float),(m.size,))
    for index,(f,q) in enumerate(zip(m.faces,qi)):
        l,r=f['left'],f['right'];wall=f['wall']
        # Harmonic mean for discontinuous material/eddy diffusivity.
        dl=diffusivity[l];dr=diffusivity[r]
        d=np.divide(2*dl*dr,dl+dr,out=np.zeros_like(dl),where=(dl+dr)>0)
        k=d*f['area']/f['distance'];k[wall]=0.
        conductances.append(k)
        outgoing_l=dt*(np.maximum(q,0)+k);outgoing_r=dt*(np.maximum(-q,0)+k)
        np.add.at(diag,l,outgoing_l);np.add.at(diag,r,outgoing_r)
        rows.extend([l,r]);cols.extend([r,l]);vals.extend([-outgoing_r,-outgoing_l])
        if wall_values is not None:
            wv=wall_values[index]
            # Separate one-sided viscous wall fluxes; never gas diffusion through a panel.
            xyz=m.xyz.reshape(-1,3)
            for cells in (l,r):
                distance=abs(xyz[cells,f['axis']]-f['position'][:,f['axis']])
                wk=dt*diffusivity[cells]*f['area']/distance*wall
                np.add.at(diag,cells,wk)
                np.add.at(rhs,cells,wk[:,None]*wv)
    bk=[]
    for index,(f,q) in enumerate(zip(m.boundaries,qb)):
        cells=f['cell'];bc=np.broadcast_to(boundary_values[index],(len(cells),ncomp))
        k=diffusivity[cells]*f['area']/f['distance']
        if wall_values is None:k=np.where(q<0,k,0.) # scalar: inflow concentration, outward zero gradient
        elif f['kind']=='open':k=np.zeros_like(k) # momentum: zero normal gradient at open faces
        bk.append(k)
        np.add.at(diag,cells,dt*(np.maximum(q,0)+k))
        np.add.at(rhs,cells,dt*(np.maximum(-q,0)+k)[:,None]*bc)
    rows.append(np.arange(m.size));cols.append(np.arange(m.size));vals.append(diag)
    a=coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(m.size,m.size)).tocsr()
    pre=LinearOperator(a.shape,matvec=lambda x:x/diag)
    result=np.empty_like(previous);iterations=[];fallbacks=[]
    for j in range(ncomp):
        counter=[0]
        def count(_):counter[0]+=1
        out,info=bicgstab(a,rhs[:,j],x0=previous[:,j],M=pre,rtol=3e-12,atol=2e-13,maxiter=350,callback=count)
        if info<0:
            # BiCGSTAB can break down for a symmetric/zero-front Krylov vector.
            # Solve the identical inventory system with GMRES; change no flux,
            # concentration, right-hand side, or conservation constraint.
            out,info=gmres(a,rhs[:,j],x0=previous[:,j],M=pre,rtol=3e-12,atol=2e-13,
                           restart=50,maxiter=350,callback=count,callback_type='legacy')
            fallbacks.append(j)
        if info:raise RuntimeError(f'Transport solve did not converge: {info}')
        result[:,j]=out;iterations.append(counter[0])
    low=result.copy()
    if correct:
        if vector:raise ValueError('Limiter is scalar-only')
        c=result[:,0];lo,hi=bounds
        if local_bounds:
            # Zalesak-style local extrema: only face-connected AIR neighbours
            # contribute. Include both old and low-order states so an implicit
            # long step is admissible, and retain the declared source bounds.
            lo_local=np.minimum(previous[:,0],c).copy()
            hi_local=np.maximum(previous[:,0],c).copy()
            for f in m.faces:
                l=f['left'][~f['wall']];r=f['right'][~f['wall']]
                np.minimum.at(lo_local,l,np.minimum(previous[r,0],c[r]))
                np.minimum.at(lo_local,r,np.minimum(previous[l,0],c[l]))
                np.maximum.at(hi_local,l,np.maximum(previous[r,0],c[r]))
                np.maximum.at(hi_local,r,np.maximum(previous[l,0],c[l]))
            for f,q,bc in zip(m.boundaries,qb,boundary_values):
                inflow=q<0
                values=np.broadcast_to(bc,(len(q),))
                np.minimum.at(lo_local,f['cell'][inflow],values[inflow])
                np.maximum.at(hi_local,f['cell'][inflow],values[inflow])
            lo=np.maximum(lo,lo_local);hi=np.minimum(hi,hi_local)
        plus=np.zeros(m.size);minus=np.zeros(m.size);antiflux=[]
        for f,q in zip(m.faces,qi):
            l,r=f['left'],f['right']
            # Correct donor-cell flux toward linear interpolation at the actual
            # face. A fixed 1/2 weight was wrong across unequal moving cells.
            face=interpolate(m,c,f)
            donor=np.where(q>=0,c[l],c[r])
            delta=dt*q*(donor-face);delta[f['wall']]=0.
            antiflux.append(delta)
            np.add.at(plus,l,np.maximum(delta,0));np.add.at(minus,l,np.maximum(-delta,0))
            np.add.at(plus,r,np.maximum(-delta,0));np.add.at(minus,r,np.maximum(delta,0))
        rp=np.minimum(1,np.maximum(0,v1*(hi-c))/np.maximum(plus,1e-300))
        rm=np.minimum(1,np.maximum(0,v1*(c-lo))/np.maximum(minus,1e-300))
        change=np.zeros(m.size)
        for f,delta in zip(m.faces,antiflux):
            l,r=f['left'],f['right']
            factor=np.where(delta>=0,np.minimum(rp[l],rm[r]),np.minimum(rm[l],rp[r]))
            limited=delta*factor
            np.add.at(change,l,limited);np.add.at(change,r,-limited)
        c+=np.divide(change,v1,out=np.zeros_like(change),where=v1>0)
    external=np.zeros(ncomp)
    for i,(f,q,k) in enumerate(zip(m.boundaries,qb,bk)):
        bc=np.broadcast_to(boundary_values[i],(len(f['cell']),ncomp));c=low[f['cell']]
        external+=dt*np.sum(np.maximum(q,0)[:,None]*c+np.minimum(q,0)[:,None]*bc+k[:,None]*(c-bc),axis=0)
    added=0 if source is None else dt*np.sum(source,axis=0)
    residual=np.sum(v1[:,None]*result-v0[:,None]*previous,axis=0)+external-added
    return (result if vector else result[:,0]),dict(balance_residual=residual.tolist(),boundary_out=external.tolist(),iterations=iterations,gmres_fallback_components=fallbacks)
