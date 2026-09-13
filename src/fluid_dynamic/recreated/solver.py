"""Unsteady three-dimensional incompressible flow and CO2 on sliding shells.

This is a verification-stage research solver. A run is not qualified by merely
finishing: refinement, domain and closure sensitivity must be examined separately.
"""
from dataclasses import asdict
import json
import time
import numpy as np
from .mesh import Mesh,opening
from .operators import Projector,face_velocity,relative_flux,transport,interpolate


def gradients(m,u,rate):
    grad=np.zeros((m.size,3,3));weights=np.zeros((m.size,3))
    xyz=m.xyz.reshape(-1,3)
    for f in m.faces:
        l,r=f['left'],f['right'];axis=f['axis'];wall=f['wall']
        du=(u[r]-u[l])/f['distance'][:,None]
        wv=m.wall_velocity(f['position'],rate)
        for cells in (l,r):
            values=du.copy()
            distance=f['position'][:,axis]-xyz[cells,axis]
            values[wall]=(wv[wall]-u[cells[wall]])/distance[wall,None]
            np.add.at(grad[:,:,axis],cells,values*f['area'][:,None])
            np.add.at(weights[:,axis],cells,f['area'])
    for f in m.boundaries:
        cells=f['cell'];axis=f['axis']
        du=np.zeros((len(cells),3))
        if f['kind']!='open':
            du=(m.boundary_wall_velocity(f,rate)-u[cells])/(f['side']*f['distance'][:,None])
        np.add.at(grad[:,:,axis],cells,du*f['area'][:,None])
        np.add.at(weights[:,axis],cells,f['area'])
    return grad/np.maximum(weights[:,None,:],1e-30)


def eddy_viscosity(m,u,rate):
    grad=gradients(m,u,rate)
    strain=.5*(grad+grad.transpose(0,2,1))
    nu=(m.cfg.smagorinsky*np.cbrt(m.volume))**2*np.sqrt(2*np.sum(strain**2,axis=(1,2)))
    return nu,grad


def cross_stress(m,nu,grad):
    """Explicit transpose-gradient part of div[nu (grad u + grad u.T)]."""
    source=np.zeros((m.size,3))
    for f in m.faces:
        l,r=f['left'],f['right'];a=f['axis'];wall=f['wall']
        stress=interpolate(m,nu[:,None]*grad[:,a,:],f)*f['area'][:,None]
        sl=stress.copy();sr=stress.copy()
        sl[wall]=(nu[l,None]*grad[l,a,:]*f['area'][:,None])[wall]
        sr[wall]=(nu[r,None]*grad[r,a,:]*f['area'][:,None])[wall]
        np.add.at(source,l,sl);np.add.at(source,r,-sr)
    for f in m.boundaries:
        cells=f['cell'];a=f['axis']
        # Homogeneous traction approximation at open outer boundaries.
        if f['kind']!='open':
            np.add.at(source,cells,nu[cells,None]*grad[cells,a,:]*f['side']*f['area'][:,None])
    return source


def run_case(cfg,destination,progress=True,mesh_factory=Mesh):
    cfg.validate();destination.mkdir(parents=True,exist_ok=False)
    (destination/'config.json').write_text(json.dumps(asdict(cfg),indent=2)+'\n')
    snapshot_dir=destination/'snapshots';snapshot_dir.mkdir()
    m=mesh_factory(cfg,opening(cfg,0.));x,y,z=m.xyz.reshape(-1,3).T
    original=(abs(x)<cfg.width/2)&(abs(y)<cfg.depth/2)&(z<cfg.height)
    if hasattr(m,'reporting_weights'):original=m.reporting_weights('fixed')>0
    c=np.where(original,cfg.inside_ppm,cfg.ambient_ppm)
    u=np.zeros((m.size,3));u[~original]=cfg.wind_vector
    cm=m.store(c);um=m.store(u)
    initial_inventory=float(np.dot(m.volume,c))
    initial_deficit=float(np.dot(m.volume,cfg.ambient_ppm-c))
    air_moles=cfg.pressure_pa/(8.31446261815324*cfg.temperature_k)
    scale=max(abs(initial_deficit),abs(cfg.tree_source_umol_s/air_moles)*cfg.duration_s,1.)
    projector=Projector();time_s=0.;step=0;records=[];frames=[]
    external=0.;source_total=0.;max_div=0.;max_gcl=0.;worst_budget=0.;native=None
    scalar_fallbacks=0
    started=time.monotonic();next_save=0.;last_print=-100.

    def save(time_s,mesh,concentration,velocity,row):
        frame=len(frames);path=snapshot_dir/f'{frame:04d}.npz'
        np.savez_compressed(path,time_s=time_s,gap=mesh.gap,x=mesh.centers[0],y=mesh.centers[1],z=mesh.centers[2],
                            x_faces=mesh.xf,y_faces=mesh.yf,z_faces=mesh.zf,
                            cell_volume_m3=mesh.volume.reshape(mesh.shape),
                            co2=concentration.reshape(mesh.shape),
                            velocity=velocity.reshape(*mesh.shape,3))
        frames.append(dict(time_s=time_s,path=str(path.relative_to(destination))))

    def diagnostics(t,mesh,cc,uu,extra):
        from scipy.interpolate import RegularGridInterpolator
        xx,yy,zz=mesh.xyz.reshape(-1,3).T
        roi=(abs(xx)<cfg.width/2)&(abs(yy)<cfg.depth/2)&(zz<cfg.height)
        # Intersection volumes avoid changing the diagnostic ROI when the mesh moves.
        wx=np.maximum(0,np.minimum(mesh.xf[1:],cfg.width/2)-np.maximum(mesh.xf[:-1],-cfg.width/2))
        wy=np.maximum(0,np.minimum(mesh.yf[1:],cfg.depth/2)-np.maximum(mesh.yf[:-1],-cfg.depth/2))
        wz=np.maximum(0,np.minimum(mesh.zf[1:],cfg.height)-np.maximum(mesh.zf[:-1],0.))
        weights=(wx[:,None,None]*wy[None,:,None]*wz[None,None,:]).ravel()
        if hasattr(mesh,'reporting_weights'):weights=mesh.reporting_weights('fixed')
        mean=float(np.dot(weights,cc)/weights.sum());std=float(np.sqrt(np.dot(weights,(cc-mean)**2)/weights.sum()))
        row=dict(time_s=t,gap_m=mesh.gap,roi_mean_ppm=mean,roi_std_ppm=std,
            roi_volume_m3=float(weights.sum()),
            minimum_ppm=float(cc.min()),maximum_ppm=float(cc.max()),
            inventory_ppm_m3=float(np.dot(mesh.volume,cc)),
            deficit_ppm_m3=float(np.dot(mesh.volume,cfg.ambient_ppm-cc)),
            max_speed_m_s=float(np.linalg.norm(uu,axis=1).max()),**extra)
        probe=RegularGridInterpolator(mesh.centers,cc.reshape(mesh.shape))
        probes=mesh.probes() if hasattr(mesh,'probes') else [('centre_1m',(0,0,1)),('centre_3m',(0,0,3)),('centre_5m',(0,0,5)),
                               ('left_1m',(-1.5,0,1)),('right_1m',(1.5,0,1))]
        for label,position in probes:
            row[f'{label}_ppm']=float(probe([position])[0])
        halfmask=(abs(xx)>mesh.gap/2)&(abs(xx)<cfg.width/2+mesh.gap/2)&(abs(yy)<cfg.depth/2)&(zz<cfg.height)
        shell_weights=mesh.volume*halfmask
        if hasattr(mesh,'reporting_weights'):shell_weights=mesh.reporting_weights('sheltered')
        for side,label in [(-1,'left'),(1,'right')]:
            w=shell_weights*(xx*side>0)
            row[f'{label}_half_mean_ppm']=float(np.dot(w,cc)/w.sum())
        shell_volume=float(shell_weights.sum())
        shell_mean=float(np.dot(shell_weights,cc)/shell_volume)
        row['shell_volume_m3']=shell_volume
        row['shell_mean_ppm']=shell_mean
        row['shell_std_ppm']=float(np.sqrt(np.dot(shell_weights,(cc-shell_mean)**2)/shell_volume))
        return row

    records.append(diagnostics(0,m,c,u,dict(budget_residual_ppm_m3=0.,budget_fraction=0.,divergence_s_inv=0.,gcl_m3=0.,max_courant=0.,eddy_diffusivity_m2_s=0.,source_umol_s=cfg.tree_source_umol_s,mouth_inflow_m3_s=0.,mouth_outflow_m3_s=0.)))
    save(0,m,c,u,records[-1]);next_save=cfg.save_every_s
    while time_s<cfg.duration_s-1e-10:
        dt=min(cfg.dt_s,cfg.duration_s-time_s)
        # Split travel and output events exactly; no step straddles a kink in wall motion.
        events=[next_save,cfg.duration_s]
        cycle=np.floor(time_s/cfg.period_s)*cfg.period_s if cfg.period_s else 0.
        for offset in ([cycle,cycle+cfg.period_s] if cfg.period_s else [0.]):
            events.extend(offset+v for v in (cfg.open_at,cfg.open_at+cfg.opening_s,cfg.close_at,cfg.close_at+cfg.closing_s))
        future=[e for e in events if e>time_s+1e-9]
        if hasattr(m,'motion_events'):future.extend(e for e in m.motion_events() if e>time_s+1e-9)
        if future:dt=min(dt,min(future)-time_s)
        g0=opening(cfg,time_s);g1=opening(cfg,time_s+dt);mid=(g0+g1)/2;rate=(g1-g0)/dt
        if m.gap!=mid:m=mesh_factory(cfg,mid)
        c=m.collect(cm);u=m.collect(um)
        # Newborn cells carry zero inventory. Their momentum guess only helps the linear solve.
        newborn=m.volumes_at(g0)==0
        u[newborn]=cfg.wind_vector;c[newborn]=cfg.ambient_ppm
        qi,qb=face_velocity(m,u,rate)
        if native is not None and native[0]==m.shape:
            qi=[v*f['area'] for v,f in zip(native[1],m.faces)]
            qb=[v*f['area'] for v,f in zip(native[2],m.boundaries)]
            for f,q in zip(m.faces,qi):q[f['wall']]=m.wall_velocity(f['position'][f['wall']],rate)[:,f['axis']]*f['area'][f['wall']]
            for f,q in zip(m.boundaries,qb):
                if f['kind']!='open':q[:]=m.boundary_wall_velocity(f,rate)[:,f['axis']]*f['side']*f['area']
        qi,qb,_,d0=projector.project(m,qi,qb,dt)
        ri,rb,gcl0=relative_flux(m,qi,qb,g0,g1,dt)
        eddy,grad=eddy_viscosity(m,u,rate);nu=cfg.viscosity+eddy
        crown,trunk=m.canopy();drag=(cfg.canopy_drag_per_m*crown+cfg.trunk_drag_per_m*trunk)*np.linalg.norm(u,axis=1)
        wall=[m.wall_velocity(f['position'],rate) for f in m.faces]
        bc=[np.broadcast_to(cfg.wind_vector,(len(f['cell']),3)) if f['kind']=='open' else m.boundary_wall_velocity(f,rate) for f in m.boundaries]
        star,_=transport(m,u,ri,rb,g0,g1,dt,nu,bc,source=cross_stress(m,nu,grad),wall_values=wall,reaction=drag)
        # Retain the pressure-corrected native face field, adding only the predictor increment.
        delta=star-u
        predicted_i=[q+interpolate(m,delta,f)[:,f['axis']]*f['area']*(~f['wall']) for f,q in zip(m.faces,qi)]
        predicted_b=[q+(delta[f['cell'],f['axis']]*f['side']*f['area'] if f['kind']=='open' else 0.) for f,q in zip(m.boundaries,qb)]
        qi,qb,pgrad,d1=projector.project(m,predicted_i,predicted_b,dt)
        u=star-dt*pgrad
        ri,rb,gcl1=relative_flux(m,qi,qb,g0,g1,dt)
        native=(m.shape,[q/f['area'] for q,f in zip(qi,m.faces)],[q/f['area'] for q,f in zip(qb,m.boundaries)])
        weights=crown*m.volume
        if cfg.tree_source_umol_s and weights.sum()<=0:raise ValueError('Nonzero plant exchange has no resolved canopy volume')
        source=weights/weights.sum()*cfg.tree_source_umol_s/air_moles if weights.sum()>0 else np.zeros(m.size)
        scalar_diffusion=cfg.diffusivity+eddy/cfg.turbulent_schmidt
        # With sources the admissible interval includes the local explicit source increment.
        v1=m.volumes_at(g1)
        split_source=bool(cfg.tree_source_umol_s and hasattr(m,'source_weights_at'))
        post_increment=np.zeros(m.size)
        if split_source:
            # Source / source-free ALE transport / source. Endpoint source
            # inventories avoid putting a finite sink in a zero-capacity cell.
            # The transport step remains first order in time; symmetric source
            # placement does not establish second-order accuracy of the solver.
            v0=m.volumes_at(g0);endpoint_sources=[]
            for gap in (g0,g1):
                w=m.source_weights_at(gap)
                if w.sum()<=0:raise ValueError('Nonzero plant exchange has no resolved endpoint canopy')
                endpoint_sources.append(w/w.sum()*cfg.tree_source_umol_s/air_moles)
            before,after=endpoint_sources
            c=c+np.divide(.5*dt*before,v0,out=np.zeros(m.size),where=v0>0)
            post_increment=np.divide(.5*dt*after,v1,out=np.zeros(m.size),where=v1>0)
            if c.min()<0:raise RuntimeError('Prescribed plant uptake exhausted local CO2; reduce timestep or use a responsive source model')
            source=np.zeros(m.size)
        increment=np.divide(dt*source,v1,out=np.zeros_like(source),where=v1>0)
        lower=min(float(c.min()),cfg.ambient_ppm)+min(0.,float(increment.min()))
        upper=max(float(c.max()),cfg.ambient_ppm)+max(0.,float(increment.max()))
        c,stats=transport(m,c,ri,rb,g0,g1,dt,scalar_diffusion,[cfg.ambient_ppm]*6,
                          source=source,correct=cfg.flux_correction,bounds=(lower,upper),
                          local_bounds=cfg.local_scalar_bounds)
        scalar_fallbacks+=len(stats['gmres_fallback_components'])
        if c.min()<lower-2e-5 or c.max()>upper+2e-5:raise RuntimeError(f'Scalar violated its admissible bounds at t={time_s:g}, dt={dt:g}, gap={g0:g}->{g1:g}: range {c.min():.10g}..{c.max():.10g}, bounds {lower:.10g}..{upper:.10g}, active range {c[v1>0].min():.10g}..{c[v1>0].max():.10g}')
        c+=post_increment
        if c.min()<0:raise RuntimeError('Prescribed plant uptake exhausted local CO2; reduce timestep or use a responsive source model')
        source_total+=dt*(cfg.tree_source_umol_s/air_moles if split_source else source.sum());external+=stats['boundary_out'][0]
        budget=float(np.dot(v1,c)-initial_inventory+external-source_total)
        worst_budget=max(worst_budget,abs(budget)/scale);max_div=max(max_div,d0,d1);max_gcl=max(max_gcl,gcl0,gcl1)
        if worst_budget>.001:raise RuntimeError(f'CO2 budget exceeds 0.1% of deficit/source scale: {worst_budget}')
        cm=m.store(c);um=m.store(u);time_s+=dt;step+=1
        outgoing=np.zeros(m.size)
        for f,q in zip(m.faces,ri):
            np.add.at(outgoing,f['left'],np.maximum(q,0));np.add.at(outgoing,f['right'],np.maximum(-q,0))
        courant=float(np.max(dt*outgoing/m.volume))
        mouth_in=0.;mouth_out=0.
        if mid>0:
            f=m.faces[0];xx,yy,zz=f['position'].T
            mouths=np.isclose(abs(xx),mid/2,atol=1e-10)&(abs(yy)<cfg.depth/2)&(zz<cfg.height)
            if hasattr(m,'mouth_selector'):mouths=m.mouth_selector(f['position'])
            flux=ri[0][mouths]*(-np.sign(xx[mouths]))
            mouth_in=float(np.maximum(-flux,0).sum());mouth_out=float(np.maximum(flux,0).sum())
        endmesh=m if m.gap==g1 else mesh_factory(cfg,g1)
        row=diagnostics(time_s,endmesh,endmesh.collect(cm),endmesh.collect(um),dict(
            budget_residual_ppm_m3=budget,budget_fraction=abs(budget)/scale,divergence_s_inv=max(d0,d1),
            gcl_m3=max(gcl0,gcl1),max_courant=courant,eddy_diffusivity_m2_s=float(scalar_diffusion.max()),source_umol_s=cfg.tree_source_umol_s,
            mouth_inflow_m3_s=mouth_in,mouth_outflow_m3_s=mouth_out))
        records.append(row)
        if time_s>=next_save-1e-8 or time_s>=cfg.duration_s-1e-8:
            save(time_s,endmesh,endmesh.collect(cm),endmesh.collect(um),row)
            next_save+=cfg.save_every_s
        if progress and time.monotonic()-last_print>15:
            print(f'{destination.name}: t={time_s:.1f}/{cfg.duration_s:g}s, mean={row["roi_mean_ppm"]:.2f} ppm, budget={worst_budget:.2g}, speed={row["max_speed_m_s"]:.2f} m/s',flush=True)
            last_print=time.monotonic()
    import csv
    with (destination/'timeseries.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=records[0]);writer.writeheader();writer.writerows(records)
    summary=dict(status='EXECUTED_UNQUALIFIED',wall_seconds=time.monotonic()-started,steps=step,cells=m.size,
        initial_deficit_ppm_m3=initial_deficit,budget_scale_ppm_m3=scale,max_budget_fraction=worst_budget,
        max_divergence_s_inv=max_div,max_gcl_m3=max_gcl,pressure_factorizations=projector.factorizations,
        pressure_refactorizations_after_cg_failure=projector.refactorizations_after_cg_failure,
        source_time_integration='endpoint half steps around ALE transport' if hasattr(m,'source_weights_at') else 'coupled implicit inventory RHS',
        scalar_gmres_fallback_count=scalar_fallbacks,
        scientific_status='DIAGNOSTIC_ONLY',
        wind_vector_m_s=cfg.wind_vector.tolist(),wind_angle_deg=cfg.wind_angle_deg,
        wind_direction_convention='flow TO; counterclockwise from +x toward +y; not meteorological wind-from',
        scalar_scheme='implicit donor + local FCT' if cfg.flux_correction and cfg.local_scalar_bounds else ('implicit donor + global FCT' if cfg.flux_correction else 'implicit donor'),
        fixed_reporting_volume_m3=records[0]['roi_volume_m3'],
        recovery_dwell_s=cfg.recovery_dwell_s,
        final_roi_mean_ppm=records[-1]['roi_mean_ppm'],final_roi_std_ppm=records[-1]['roi_std_ppm'],
        max_courant=max(r['max_courant'] for r in records),minimum_ppm=min(r['minimum_ppm'] for r in records),
        maximum_ppm=max(r['maximum_ppm'] for r in records),frames=frames)
    for fraction,label in [(1-np.exp(-1),'t63_s'),(.95,'t95_s'),(.99,'t99_s')]:
        target=cfg.inside_ppm+fraction*(cfg.ambient_ppm-cfg.inside_ppm)
        recovered=np.array([(r['roi_mean_ppm']>=target if cfg.inside_ppm<cfg.ambient_ppm else r['roi_mean_ppm']<=target) for r in records])
        persistent=np.logical_and.accumulate(recovered[::-1])[::-1]
        hit=next((r['time_s'] for r,ok in zip(records,persistent) if ok and r['time_s']<=cfg.duration_s-cfg.recovery_dwell_s),None)
        summary[label]=hit if cfg.inside_ppm!=cfg.ambient_ppm and cfg.tree_source_umol_s==0 else None
        summary[label.replace('_s','_after_open_s')]=max(0.,hit-cfg.open_at) if summary[label] is not None and not cfg.start_open else summary[label]
    (destination/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
