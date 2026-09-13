"""Previous chamber visual language applied to unchanged saved scenario fields.

No solver is invoked. Cell colors use saved values and true ALE face coordinates.
Camera, color range and playback are fixed within a comparison group. The
existing illustrative palm asset is fitted to declared display envelopes only.
"""
import json,subprocess
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from .geometry import ScenarioStudy,panels
from ..recreated.render import BG,INK,ACCENT,CO2,style,palm,series

COLORS=['#18665e','#4877bd','#b06c29','#8a59a5','#338c79','#4d91c7','#cd9760','#b485bb','#8a897e']
FPS=6


def load_case(path):
    path=Path(path)
    return dict(path=path,scenario=json.loads((path/'scenario.json').read_text()),
                cfg=ScenarioStudy(**json.loads((path/'config.json').read_text())),
                summary=json.loads((path/'summary.json').read_text()),ts=series(path/'timeseries.csv'))


def color_limits(cases):
    low=min(c['summary']['minimum_ppm'] for c in cases)
    high=max(max(c['summary']['maximum_ppm'],c['cfg'].ambient_ppm) for c in cases)
    span=max(high-low,.25);step=.25 if span<=5 else (1. if span<=20 else 5.)
    return float(np.floor(low/step+1e-8)*step),float(np.ceil(high/step-1e-8)*step)


def repeats(gap,maximum,index):
    return 8 if index==0 else (6 if 1e-8<gap<maximum-1e-8 else 2)


def synchronized_frames(cases):
    frames=[c['summary']['frames'] for c in cases]
    times=[[f['time_s'] for f in row] for row in frames]
    if not all(len(t)==len(times[0]) and np.allclose(t,times[0],rtol=0,atol=1e-9) for t in times):
        raise ValueError('Comparison requires identical saved physical timestamps; no field interpolation is allowed')
    return frames


def extents(cfg):
    g=cfg.geometry
    extent=g.get('fixed_end_m',cfg.width/2+cfg.gap/2)
    depth=g.get('fixed_depth_m',cfg.depth)/2
    roof=g.get('fixed_roof_z_m',g['floor_z_m']+cfg.height)
    pad=.20*max(2*extent,2*depth,roof)
    return extent+pad,depth+pad,roof+.55*pad


def panel_polygons(p):
    """Quads leave the fixed-to-moving aperture open in the rendering too."""
    lo=np.array(p['lo'],float);hi=np.array(p['hi'],float)
    axes=[i for i in range(3) if i!=p['axis']];a,b=axes
    rectangles=[(lo,hi)]
    if p['hole']:
        hl,hh=map(np.array,p['hole']);rectangles=[]
        for lower,upper in [(lo,np.where(np.arange(3)==a,hl,hi)),
                            (np.where(np.arange(3)==a,hh,lo),hi),
                            (np.where(np.arange(3)==a,hl,lo),np.where(np.arange(3)==a,hh,np.where(np.arange(3)==b,hl,hi))),
                            (np.where(np.arange(3)==a,hl,np.where(np.arange(3)==b,hh,lo)),np.where(np.arange(3)==a,hh,hi))]:
            if upper[a]>lower[a] and upper[b]>lower[b]:rectangles.append((lower,upper))
    polygons=[]
    for lower,upper in rectangles:
        points=[]
        for x,y in [(lower[a],lower[b]),(upper[a],lower[b]),(upper[a],upper[b]),(lower[a],upper[b])]:
            point=np.zeros(3);point[p['axis']]=p['position'];point[a]=x;point[b]=y;points.append(point)
        polygons.append(points)
    return polygons


def shells(ax,cfg,gap):
    for p in panels(cfg,gap):
        color=(.18,.35,.35,.8) if p['speed'] else (.47,.40,.32,.8)
        ax.add_collection3d(Poly3DCollection(panel_polygons(p),facecolors=(.45,.62,.62,.025),edgecolors=color,linewidths=.9,zorder=3))


def plan_panels(ax,cfg,gap,z):
    for p in panels(cfg,gap):
        if p['axis']==2:continue
        for poly in panel_polygons(p):
            poly=np.array(poly)
            if poly[:,2].min()-1e-9<=z<=poly[:,2].max()+1e-9:
                ax.plot([poly[:,0].min(),poly[:,0].max()],[poly[:,1].min(),poly[:,1].max()],
                        color=INK if p['speed'] else '#897052',lw=1.6)


def fitted_palm(root,cfg):
    if cfg.plant['height_m']==0:return []
    polys=palm(root)
    # build_synthetic_palm.py declares 64 trunk triangles, 5*32 frond
    # triangles, then a fruit illustration. Keep the known trunk/frond asset;
    # omit the fruit illustration because no fruit trait is supplied here.
    if len(polys)!=544:raise ValueError('Illustrative palm topology changed; inspect its component mapping')
    p=cfg.plant;h=p['height_m'];result=[];frond=np.concatenate(polys[64:224]);radius=np.max(abs(frond[:,:2]),axis=0)
    for i,poly in enumerate(polys[:224]):
        q=poly.copy()
        if i<64:
            q[:,:2]*=p['trunk_radius_m']/.15;q[:,2]*=.45*h/1.8
        else:
            q[:,:2]*=np.array(p['crown_radii_m'])/radius
            q[:,2]=.45*h+(q[:,2]-1.8)*(.55*h)/(frond[:,2].max()-1.8)
        q[:,2]+=cfg.geometry['floor_z_m'];result.append(q)
    return result


def crop_indices(centres,low,high):
    ids=np.flatnonzero((centres>=low)&(centres<=high))
    if not len(ids):raise ValueError('Displayed section contains no saved cells')
    return ids


def draw_scene(ax,d,cfg,norm,tree):
    x,y,z=d['x'],d['y'],d['z'];c=d['co2'];u=d['velocity'];ex,ey,ez=extents(cfg)
    iy=int(np.argmin(abs(y)));iz=int(np.argmin(abs(z-(cfg.geometry['floor_z_m']+.25*cfg.height))))
    ix=crop_indices(x,-ex,ex);jy=crop_indices(y,-ey,ey);kz=crop_indices(z,0,ez)
    xf=d['x_faces'][ix[0]:ix[-1]+2];yf=d['y_faces'][jy[0]:jy[-1]+2];zf=d['z_faces'][kz[0]:kz[-1]+2]
    xx,zz=np.meshgrid(xf,zf,indexing='ij')
    ax.plot_surface(xx,np.full_like(xx,y[iy]),zz,facecolors=CO2(norm(c[np.ix_(ix,[iy],kz)][:,0,:])),
                    rstride=1,cstride=1,shade=False,antialiased=False,alpha=.52,linewidth=0,zorder=1)
    xx,yy=np.meshgrid(xf,yf,indexing='ij')
    ax.plot_surface(xx,yy,np.full_like(xx,z[iz]),facecolors=CO2(norm(c[np.ix_(ix,jy,[iz])][:,:,0])),
                    rstride=1,cstride=1,shade=False,antialiased=False,alpha=.48,linewidth=0,zorder=1)
    shells(ax,cfg,float(d['gap']))
    if tree:ax.add_collection3d(Poly3DCollection(tree,facecolor=(.05,.24,.12,.95),edgecolor='none',zorder=4))
    qx,qz=np.meshgrid(ix[::max(1,len(ix)//8)],kz[::max(1,len(kz)//5)],indexing='ij')
    ax.quiver(x[qx],np.full(qx.shape,y[iy]),z[qz],u[qx,iy,qz,0],u[qx,iy,qz,1],u[qx,iy,qz,2],
              length=.5,normalize=False,color=INK,linewidth=.65,arrow_length_ratio=.25,zorder=5)
    ax.set(xlim=(-ex,ex),ylim=(-ey,ey),zlim=(0,ez),xlabel='x (m)',ylabel='y (m)',zlabel='z (m)')
    ax.view_init(elev=22,azim=-64);ax.set_box_aspect((2*ex,2*ey,ez));ax.grid(False)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=False
    ax.set_title(f'Saved CO₂ slices · y={y[iy]:.2f} m, z={z[iz]:.2f} m',fontsize=10,pad=-7)


def draw_plan(ax,d,cfg,norm,title=None):
    x,y=d['x'],d['y'];ex,ey,_=extents(cfg)
    iz=int(np.argmin(abs(d['z']-(cfg.geometry['floor_z_m']+.25*cfg.height))))
    c=d['co2'][:,:,iz];u=d['velocity'][:,:,iz,:]
    ax.pcolormesh(d['x_faces'],d['y_faces'],c.T,cmap=CO2,norm=norm,shading='flat')
    ix=crop_indices(x,-ex,ex);jy=crop_indices(y,-ey,ey)
    qx,qy=np.meshgrid(ix[::max(1,len(ix)//11)],jy[::max(1,len(jy)//9)],indexing='ij')
    ax.quiver(x[qx],y[qy],u[qx,qy,0],u[qx,qy,1],color=INK,scale_units='xy',angles='xy',scale=2,width=.0025)
    plan_panels(ax,cfg,float(d['gap']),float(d['z'][iz]))
    theta=np.deg2rad(cfg.wind_angle_deg)
    if cfg.wind:
        a=np.array([.68*ex,.72*ey]);b=a+min(ex,ey)*.23*np.array([np.cos(theta),np.sin(theta)])
        ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle='->',color='#af5128',lw=2))
    ax.set(xlim=(-ex,ex),ylim=(-ey,ey),xlabel='x (m)',ylabel='y (m)');ax.set_aspect('equal',adjustable='box')
    ax.set_title(title or f'Top view at z={d["z"][iz]:.2f} m · {cfg.wind:g} m/s toward {cfg.wind_angle_deg:g}°',fontsize=10)


def draw_history(ax,case,t):
    ts=case['ts'];cfg=case['cfg'];tt=ts['time_s'];mean=ts['roi_mean_ppm'];sd=ts['roi_std_ppm']
    ax.fill_between(tt,mean-sd,mean+sd,color=ACCENT,alpha=.12,label='Spatial ±1σ')
    ax.plot(tt,mean,color=ACCENT,lw=1.8,label=f'Fixed {ts["roi_volume_m3"][0]:.1f} m³ mean')
    for key,color,ls,label in [('left_half_mean_ppm','#7659a5','--','Left sheltered mean'),('right_half_mean_ppm','#b57834',':','Right sheltered mean')]:
        ax.plot(tt,ts[key],color=color,ls=ls,lw=1.1,label=label)
    ax.axhline(cfg.ambient_ppm,color='#9f9b8c',ls='--',lw=1,label='Ambient')
    ax.axvline(t,color=INK,lw=.8);i=int(np.argmin(abs(tt-t)));ax.scatter([t],[mean[i]],s=24,color=INK,zorder=4)
    ax.set(xlim=(0,cfg.duration_s),xlabel='Physical time (s)',ylabel='CO₂ (ppm)')
    ax.legend(loc='lower right',fontsize=7,frameon=False,ncol=2);ax.grid(alpha=.16)


def finish_movie(tmp,path,expected):
    subprocess.run(['ffmpeg','-v','error','-i',str(tmp),'-f','null','-'],check=True,capture_output=True)
    data=json.loads(subprocess.check_output(['ffprobe','-v','error','-show_entries','stream=width,height,nb_frames:format=duration','-of','json',str(tmp)],text=True))
    assert int(data['streams'][0]['nb_frames'])==expected,(path,expected,data)
    tmp.replace(path)
    return dict(path=path.name,encoded_frames=expected,fps=FPS,duration_s=float(data['format']['duration']),
                width=data['streams'][0]['width'],height=data['streams'][0]['height'],full_decode='PASS')


def write_movie(case,out,root,name,limits,plan_only=False,preview=False):
    style();cfg=case['cfg'];s=case['scenario'];norm=Normalize(*limits);tree=fitted_palm(root,cfg)
    for folder in ['videos','figures']:(out/folder).mkdir(exist_ok=True,parents=True)
    title='Opening and air refresh' if cfg.tree_source_umol_s==0 else ('Closed → open → closed' if s['family']=='cycle' else 'Plant exchange in the open chamber' if cfg.start_open else 'Plant exchange during opening')
    fig=plt.figure(figsize=(14,8),facecolor=BG)
    gs=fig.add_gridspec(2,2,left=.025,right=.900,bottom=.16,top=.84,width_ratios=[1.18,1],hspace=.42,wspace=.16)
    ax3=fig.add_subplot(gs[:,0]) if plan_only else fig.add_subplot(gs[:,0],projection='3d',computed_zorder=False)
    axp=fig.add_subplot(gs[0,1]);axt=fig.add_subplot(gs[1,1])
    fig.text(.04,.955,'WHOLE-TREE CHAMBER  /  RECREATED STUDY',fontsize=10,weight='bold',color=ACCENT)
    fig.text(.04,.908,title,fontsize=23,weight='bold')
    fig.text(.04,.865,s['geometry']['label']+' · '+s['plant']['label'],fontsize=10)
    fig.text(.95,.955,'NUMERICAL DIAGNOSTIC • UNQUALIFIED',ha='right',fontsize=9,color='#996320')
    stamp=fig.text(.95,.906,'',ha='right',fontsize=11)
    cax=fig.add_axes([.925,.39,.01,.37]);fig.colorbar(ScalarMappable(norm=norm,cmap=CO2),cax=cax,label='CO₂ (ppm); fixed scale',ticks=np.linspace(*limits,5))
    cax.tick_params(labelsize=8);cax.yaxis.label.set_size(9)
    footer=fig.text(.04,.080,'',fontsize=9)
    fig.text(.04,.025,f'{cfg.width:g} × {cfg.depth:g} × {cfg.height:g} m moving envelope · '+('empty control' if not tree else 'illustrative palm; porous-canopy flow model')+' · saved physical time · travel playback slowed',fontsize=9,color='#60736d')
    frames=case['summary']['frames'];selected=[len(frames)//2] if preview else range(len(frames))
    path=out/'videos'/f'{name}.mp4';tmp=out/'videos'/f'{name}.encoding.mp4'
    writer=FFMpegWriter(fps=FPS,codec='libx264',extra_args=['-pix_fmt','yuv420p','-crf','20','-movflags','+faststart','-threads','1'])
    count=0;provenance=[]
    from contextlib import nullcontext
    with (nullcontext() if preview else writer.saving(fig,str(tmp),dpi=100)):
        for j in selected:
            frame=frames[j]
            with np.load(case['path']/frame['path']) as d:
                for ax in [ax3,axp,axt]:ax.clear();ax.set_facecolor(BG)
                if plan_only:draw_plan(ax3,d,cfg,norm)
                else:draw_scene(ax3,d,cfg,norm,tree)
                if plan_only:
                    iy=int(np.argmin(abs(d['y'])));axp.pcolormesh(d['x_faces'],d['z_faces'],d['co2'][:,iy,:].T,cmap=CO2,norm=norm,shading='flat')
                    axp.set(xlabel='x (m)',ylabel='z (m)',title=f'Vertical saved section y={d["y"][iy]:.2f} m');axp.set_aspect('equal')
                else:draw_plan(axp,d,cfg,norm)
                t=float(d['time_s']);gap=float(d['gap']);draw_history(axt,case,t)
            state='SEALED' if gap<1e-8 else ('FULLY OPEN' if gap>cfg.gap-1e-8 else 'PANELS MOVING')
            stamp.set_text(f'{t:.1f} s · {state}\ngap {gap:.2f} m')
            i=int(np.argmin(abs(case['ts']['time_s']-t)))
            rate=f'Prescribed net exchange {cfg.tree_source_umol_s:.2f} µmol/s' if tree else 'Source-free tracer refresh'
            source_note=case.get('source_note','') if tree else ''
            footer.set_text(f'{rate}'+('  |  '+source_note if source_note else '')+
                f'\nFixed mean {case["ts"]["roi_mean_ppm"][i]:.2f} ppm  |  gas-budget error {100*case["ts"]["budget_fraction"][i]:.2g}% of reference scale')
            if preview or j in {0,len(frames)//4,len(frames)//2,len(frames)-1}:
                fig.savefig(out/'figures'/f'{name}_{j:04d}.png',dpi=100,facecolor=BG)
            n=repeats(gap,cfg.gap,j);provenance.append(dict(saved_frame=frame['path'],physical_time_s=t,first_encoded_frame=count,repeats=n));count+=n
            if not preview:
                for _ in range(n):writer.grab_frame(facecolor=BG)
    plt.close(fig)
    if preview:return dict(preview=str(out/'figures'/f'{name}_{selected[0]:04d}.png'))
    result=finish_movie(tmp,path,count)
    result.update(source_case=s['id'],color_limits_ppm=list(limits),frame_mapping=provenance,
                  source_note=case.get('source_note',''),color_normalization='linear; fixed over the whole movie',
                  color_map=CO2.name,
                  color_palette=[matplotlib.colors.to_hex(CO2(v)) for v in np.linspace(0,1,9)],
                  data_treatment='Saved cell values on saved faces; no concentration/time interpolation',kind='2d' if plan_only else '3d')
    return result


def comparison_movie(cases,out,name,title,design_comparison=False):
    style();allframes=synchronized_frames(cases);limits=color_limits(cases);norm=Normalize(*limits)
    fig,axes=plt.subplots(2,3 if design_comparison else 4,figsize=(16,9),facecolor=BG)
    fig.subplots_adjust(left=.055,right=.92,bottom=.16,top=.82,wspace=.30,hspace=.46)
    fig.text(.045,.957,'WHOLE-TREE CHAMBER  /  '+('DESIGN COMPARISON' if design_comparison else 'EIGHT WIND DIRECTIONS'),fontsize=10,weight='bold',color=ACCENT)
    fig.text(.045,.911,title,fontsize=23,weight='bold')
    stamp=fig.text(.93,.951,'',ha='right',fontsize=12)
    fig.text(.045,.865,('Each design has its own spatial axes; all share the CO₂ scale and physical timestamps.' if design_comparison else cases[0]['scenario']['geometry']['label']+' · '+cases[0]['scenario']['plant']['label']),fontsize=10)
    fig.text(.045,.07,'Saved CO₂ and computed velocity · 0° = flow toward +x; 90° = +y · physical time stamped · playback slowed during travel',fontsize=9)
    if cases[0].get('source_note'):
        fig.text(.045,.048,f'Net exchange {cases[0]["cfg"].tree_source_umol_s:.2f} µmol/s · '+cases[0]['source_note'],fontsize=9)
    fig.text(.045,.027,'Diagnostic only · geometry, prescribed source and numerical sensitivity limits remain unchanged',fontsize=9,color='#60736d')
    cax=fig.add_axes([.945,.25,.012,.46]);fig.colorbar(ScalarMappable(norm=norm,cmap=CO2),cax=cax,label='CO₂ (ppm); fixed scale')
    path=out/'videos'/f'{name}.mp4';tmp=path.with_name(name+'.encoding.mp4')
    writer=FFMpegWriter(fps=FPS,codec='libx264',extra_args=['-pix_fmt','yuv420p','-crf','20','-movflags','+faststart','-threads','1'])
    count=0;mapping=[]
    with writer.saving(fig,str(tmp),dpi=100):
        for j,frame in enumerate(allframes[0]):
            t=frame['time_s'];moving=False
            for ax,case,frames in zip(axes.flat,cases,allframes):
                cfg=case['cfg'];ts=case['ts'];ax.clear();ax.set_facecolor(BG)
                with np.load(case['path']/frames[j]['path']) as d:
                    gap=float(d['gap']);moving|=1e-8<gap<cfg.gap-1e-8
                    mean=float(np.interp(t,ts['time_s'],ts['shell_mean_ppm']))
                    label=case['scenario']['geometry']['label'] if design_comparison else f'{cfg.wind_angle_deg:g}°'
                    draw_plan(ax,d,cfg,norm,title=f'{label}\nSheltered mean {mean:.2f} ppm')
            if design_comparison:
                ax=axes.flat[-1];ax.clear();ax.set_facecolor(BG)
                for case,color in zip(cases,COLORS):ax.plot(case['ts']['time_s'],case['ts']['roi_mean_ppm'],label=case['scenario']['geometry']['id'],color=color,lw=1.4)
                ax.axvline(t,color=INK,lw=.8);ax.set(title='Original-region means',xlabel='Physical time (s)',ylabel='CO₂ (ppm)');ax.legend(fontsize=7,frameon=False);ax.grid(alpha=.15)
            stamp.set_text(f't = {t:.1f} s')
            if j in {0,len(allframes[0])//2,len(allframes[0])-1}:fig.savefig(out/'figures'/f'{name}_{round(t):03d}s.png',dpi=100,facecolor=BG)
            n=8 if j==0 else 6 if moving else 2;mapping.append(dict(physical_time_s=t,first_encoded_frame=count,repeats=n));count+=n
            for _ in range(n):writer.grab_frame(facecolor=BG)
    plt.close(fig);result=finish_movie(tmp,path,count)
    result.update(source_cases=[c['scenario']['id'] for c in cases],color_limits_ppm=list(limits),frame_mapping=mapping,
                  color_map=CO2.name,color_normalization='linear; fixed over the whole movie',
                  color_palette=[matplotlib.colors.to_hex(CO2(v)) for v in np.linspace(0,1,9)],
                  kind='design_comparison' if design_comparison else 'wind_comparison')
    return result


def gif(video):
    path=video.with_suffix('.gif');tmp=path.with_name(path.stem+'.encoding.gif')
    subprocess.run(['ffmpeg','-v','error','-i',str(video),'-filter_complex','[0:v]fps=6,scale=960:-1:flags=lanczos,split[a][b];[a]palettegen[p];[b][p]paletteuse','-loop','0',str(tmp)],check=True,capture_output=True)
    subprocess.run(['ffmpeg','-v','error','-i',str(tmp),'-f','null','-'],check=True,capture_output=True)
    tmp.replace(path);return path
