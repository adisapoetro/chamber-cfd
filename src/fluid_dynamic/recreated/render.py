"""Render only saved solution fields, with fixed scales and explicit assumptions."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

BG='#f5f3eb';INK='#203c38';ACCENT='#1a7862'
# Owner-selected standard viridis: purple = low CO₂, yellow = high CO₂.
# The caller retains the fixed physical concentration limits and normalization.
CO2=matplotlib.colormaps['viridis']


def series(path):
    with path.open() as stream:rows=list(csv.DictReader(stream))
    return {k:np.array([float(r[k]) for r in rows]) for k in rows[0]}


def palm(root,cfg=None):
    vertices=[];faces=[]
    for line in (root/'meshes/palm_libz.obj').read_text().splitlines():
        parts=line.split()
        if not parts:continue
        if parts[0]=='v':vertices.append([float(x) for x in parts[1:4]])
        if parts[0]=='f':faces.append([int(x.split('/')[0])-1 for x in parts[1:]])
    v=np.array(vertices)
    if cfg and cfg.get('plant'):
        p=cfg['plant'];v[:,2]-=v[:,2].min()
        v[:,2]*=p['height_m']/v[:,2].max()
        v[:,2]+=cfg.get('geometry',{}).get('floor_z_m',0.)
        for axis in (0,1):v[:,axis]*=p['crown_radii_m'][axis]/max(abs(v[:,axis]))
    return [v[f] for f in faces]


def geometry_values(cfg=None):
    cfg=cfg or dict(width=4.,depth=4.,height=6.)
    return cfg['width'],cfg['depth'],cfg['height'],cfg.get('geometry',{}).get('floor_z_m',0.)


def shell_polygons(gap,cfg=None):
    width,depth,height,floor=geometry_values(cfg);d=depth/2;top=floor+height
    panels=[]
    for side in [-1,1]:
        a=side*gap/2;b=side*(width/2+gap/2)
        panels.extend([[(b,-d,floor),(b,d,floor),(b,d,top),(b,-d,top)],
                       [(a,-d,floor),(b,-d,floor),(b,-d,top),(a,-d,top)],
                       [(a,d,floor),(b,d,floor),(b,d,top),(a,d,top)],
                       [(a,-d,top),(b,-d,top),(b,d,top),(a,d,top)]])
        if floor>0:panels.append([(a,-d,floor),(b,-d,floor),(b,d,floor),(a,d,floor)])
    return panels


def shells(ax,gap,cfg=None):
    panels=shell_polygons(gap,cfg)
    ax.add_collection3d(Poly3DCollection(panels,facecolors=(.45,.62,.62,.035),edgecolors=(.18,.35,.35,.8),linewidths=1.,zorder=3))


def plan_panels(ax,gap,cfg=None):
    width,depth,_,_=geometry_values(cfg)
    for side in [-1,1]:
        a=side*gap/2;b=side*(width/2+gap/2)
        ax.plot([a,b,b,a],[-depth/2,-depth/2,depth/2,depth/2],color=INK,lw=2)


def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelcolor':INK,
                         'text.color':INK,'xtick.color':INK,'ytick.color':INK,
                         'axes.spines.top':False,'axes.spines.right':False})


def render_case(case,out,root,name):
    style();summary=json.loads((case/'summary.json').read_text());cfg=json.loads((case/'config.json').read_text())
    ts=series(case/'timeseries.csv');tree=palm(root,cfg)
    width,depth,height,floor=geometry_values(cfg);volume=width*depth*height
    xmax=(width+cfg['gap'])/2+2;ymax=depth/2+1.2;zmax=floor+height+1
    plane_z=floor+.5*height
    source=cfg['tree_source_umol_s']!=0;cycle=cfg['period_s']>0
    title='Opening and air refresh' if not source else ('Two measurement cycles' if cycle else 'Continuous plant exchange in the open chamber')
    low,high=(340.,400.) if not source else ((360.,400.) if cycle else (396.,400.))
    norm=Normalize(low,high);cmap=CO2
    fig=plt.figure(figsize=(14,8),facecolor=BG)
    gs=fig.add_gridspec(2,2,left=.025,right=.925,bottom=.16,top=.86,width_ratios=[1.18,1],hspace=.34,wspace=.14)
    ax3=fig.add_subplot(gs[:,0],projection='3d',computed_zorder=False);axp=fig.add_subplot(gs[0,1]);axt=fig.add_subplot(gs[1,1])
    fig.text(.04,.952,'WHOLE-TREE CHAMBER  /  RECREATED STUDY',fontsize=10,weight='bold',color=ACCENT)
    fig.text(.04,.905,title,fontsize=23,weight='bold')
    status=fig.text(.95,.951,'NUMERICAL DIAGNOSTIC • UNQUALIFIED',ha='right',fontsize=9,color='#996320')
    stamp=fig.text(.95,.9,'',ha='right',fontsize=13)
    cax=fig.add_axes([.945,.39,.01,.39]);fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=cax,label='CO₂ (ppm); fixed scale',ticks=np.linspace(low,high,5),extend='min' if source else 'neither')
    cax.tick_params(labelsize=9);cax.yaxis.label.set_size(9)
    footer=fig.text(.04,.067,'',fontsize=10)
    angle=cfg.get('wind_angle_deg',0.)
    wind_label=f'wind {cfg["wind"]:g} m/s toward {angle:g}°' if cfg['wind'] else 'calm ambient'
    fig.text(.04,.027,f'{width:g} × {depth:g} × {height:g} m airspace · floor {floor:g} m · synthetic palm · {wind_label} · diagnostic only',fontsize=9,color='#60736d')
    writer=FFMpegWriter(fps=6,codec='libx264',bitrate=2800,extra_args=['-pix_fmt','yuv420p','-movflags','+faststart','-threads','1'])
    (out/'videos').mkdir(exist_ok=True);(out/'figures').mkdir(exist_ok=True)
    video=out/'videos'/f'{name}.mp4'
    with writer.saving(fig,str(video),dpi=100):
        for j,frame in enumerate(summary['frames']):
            with np.load(case/frame['path']) as d:
                x=d['x'];y=d['y'];z=d['z'];c=d['co2'];u=d['velocity'];gap=float(d['gap']);t=float(d['time_s'])
                # New bundles retain true face coordinates on the unequal ALE cells.
                xf=d['x_faces'] if 'x_faces' in d else None
                yf=d['y_faces'] if 'y_faces' in d else None
            iy=int(np.argmin(abs(y)));iz=int(np.argmin(abs(z-plane_z)))
            for ax in [ax3,axp,axt]:ax.clear();ax.set_facecolor(BG)
            ix3=np.flatnonzero(abs(x)<xmax);iy3=np.flatnonzero(abs(y)<ymax);iz3=np.flatnonzero(z<zmax)
            xx,zz=np.meshgrid(x[ix3],z[iz3],indexing='ij')
            ax3.plot_surface(xx,np.zeros_like(xx),zz,facecolors=cmap(norm(c[np.ix_(ix3,[iy],iz3)][:,0,:])),rstride=1,cstride=1,
                             shade=False,antialiased=False,alpha=.48,linewidth=0,edgecolor='none',zorder=1)
            xx,yy=np.meshgrid(x[ix3],y[iy3],indexing='ij')
            ax3.plot_surface(xx,yy,np.full_like(xx,z[iz]),facecolors=cmap(norm(c[np.ix_(ix3,iy3,[iz])][:,:,0])),rstride=1,cstride=1,
                             shade=False,antialiased=False,alpha=.45,linewidth=0,edgecolor='none',zorder=1)
            shells(ax3,gap,cfg)
            ax3.add_collection3d(Poly3DCollection(tree,facecolor=(.05,.24,.12,.95),edgecolor='none',zorder=4))
            step=max(1,round(cfg['resolution']))
            selx=np.flatnonzero(abs(x)<xmax-.3)[::step];selz=np.flatnonzero((z>floor+.2)&(z<floor+height))[::2*step]
            qx,qz=np.meshgrid(selx,selz,indexing='ij')
            ax3.quiver(x[qx],np.full(qx.shape,-.03),z[qz],u[qx,iy,qz,0],u[qx,iy,qz,1],u[qx,iy,qz,2],
                       length=.5,normalize=False,color=INK,linewidth=.7,arrow_length_ratio=.25,zorder=5)
            ax3.set(xlim=(-xmax,xmax),ylim=(-ymax,ymax),zlim=(0,zmax),xlabel='x (m)',ylabel='y (m)',zlabel='z (m)')
            ax3.view_init(elev=22,azim=-64);ax3.set_box_aspect((2*xmax,2*ymax,zmax));ax3.grid(False)
            for axis in [ax3.xaxis,ax3.yaxis,ax3.zaxis]:axis.pane.fill=False
            ax3.set_title('Saved CO₂ slices and computed velocity',fontsize=11,pad=-7)
            if xf is None:
                axp.pcolormesh(x,y,c[:,:,iz].T,cmap=cmap,norm=norm,shading='nearest',rasterized=True)
            else:
                axp.pcolormesh(xf,yf,c[:,:,iz].T,cmap=cmap,norm=norm,shading='flat',rasterized=True)
            selx=np.flatnonzero(abs(x)<xmax)[::step];sely=np.flatnonzero(abs(y)<ymax)[::step]
            qx,qy=np.meshgrid(selx,sely,indexing='ij')
            axp.quiver(x[qx],y[qy],u[qx,qy,iz,0],u[qx,qy,iz,1],color=INK,scale=20,width=.0025)
            plan_panels(axp,gap,cfg);axp.set(xlim=(-xmax,xmax),ylim=(-ymax,ymax),xlabel='x (m)',ylabel='y (m)')
            axp.set_aspect('equal',adjustable='box');axp.set_title(f'Top view at z = {z[iz]:.2f} m  |  {wind_label}',fontsize=11)
            theta=np.deg2rad(angle)
            if cfg['wind']:
                axp.annotate('',xy=(4+0.65*np.cos(theta),2+0.65*np.sin(theta)),xytext=(4,2),
                             arrowprops=dict(arrowstyle='->',color='#af5128',lw=2))
            tt=ts['time_s'];mean=ts['roi_mean_ppm'];std=ts['roi_std_ppm']
            axt.fill_between(tt,mean-std,mean+std,color=ACCENT,alpha=.12,label='Spatial ±1σ')
            axt.plot(tt,mean,color=ACCENT,lw=1.8,label=f'Fixed {volume:g} m³ mean')
            if 'left_half_mean_ppm' in ts:
                axt.plot(tt,ts['left_half_mean_ppm'],color='#7659a5',ls='--',lw=1,label='Left shell mean')
                axt.plot(tt,ts['right_half_mean_ppm'],color='#b57834',ls=':',lw=1.2,label='Right shell mean')
            axt.axhline(400,color='#9f9b8c',ls='--',lw=1,label='Ambient')
            axt.axvline(t,color=INK,lw=.8)
            current=int(np.argmin(abs(tt-t)));axt.scatter([t],[mean[current]],s=24,color=INK,zorder=4)
            axt.set(xlim=(0,cfg['duration_s']),xlabel='Physical time (s)',ylabel='CO₂ (ppm)')
            axt.set_ylim(min(low,float((mean-std).min())-1),401.)
            axt.legend(loc='lower right',fontsize=7,frameon=False,ncol=2)
            axt.grid(alpha=.16)
            state='SEALED' if gap<1e-8 else ('FULLY OPEN' if gap>cfg['gap']-1e-8 else 'PANELS MOVING')
            stamp.set_text(f'{t:6.1f} s  ·  {state}  ·  gap {gap:.2f} m')
            source_text=f'Continuous synthetic source: {cfg["tree_source_umol_s"]:g} µmol/s' if source else 'Source-free tracer refresh'
            footer.set_text(f'{source_text}   |   mean {mean[current]:.2f} ppm   |   cumulative gas-budget error {100*ts["budget_fraction"][current]:.2g}% of reference scale')
            if j in {0,len(summary['frames'])//4,len(summary['frames'])//2,len(summary['frames'])-1}:
                fig.savefig(out/'figures'/f'{name}_{j:04d}.png',dpi=100,facecolor=BG)
            repeats=6 if j==0 else (3 if 1e-8<gap<cfg['gap']-1e-8 else 1)
            for _ in range(repeats):writer.grab_frame(facecolor=BG)
            if j%30==0:print(f'Render {name}: {j+1}/{len(summary["frames"])}',flush=True)
    plt.close(fig)
    return video


def legacy_comparison(out,root,case):
    """Nearest legacy output at stated times, identical spatial/colour axes."""
    style();legacy=root/'archive/legacy_inputs/20260608-112333'
    times=[0.,20.,40.,120.]
    snapshots=[]
    for p in sorted(legacy.glob('snap_*.npz')):
        with np.load(p) as d:snapshots.append((float(d['t']),p))
    summary=json.loads((case/'summary.json').read_text())
    fig,axes=plt.subplots(2,4,figsize=(15,7.5),facecolor=BG,sharex=True,sharey=True)
    norm=Normalize(340,400)
    for k,t in enumerate(times):
        oldt,p=min(snapshots,key=lambda pair:abs(pair[0]-t))
        with np.load(p) as d:
            c=d['c'];x=-7+(np.arange(c.shape[0])+.5)*.2;z=(np.arange(c.shape[2])+.5)*.2
            axes[0,k].pcolormesh(x,z,c[:,c.shape[1]//2,:].T,cmap=CO2,norm=norm,shading='nearest')
        f=min(summary['frames'],key=lambda f:abs(f['time_s']-t))
        with np.load(case/f['path']) as d:
            c=d['co2'];axes[1,k].pcolormesh(d['x'],d['z'],c[:,len(d['y'])//2,:].T,cmap=CO2,norm=norm,shading='nearest')
        for row in [0,1]:
            ax=axes[row,k];ax.set(xlim=(-5,5),ylim=(0,7));ax.set_aspect('equal');ax.set_facecolor(BG)
            gap=2*np.clip((t-10)/20,0,1)
            for side in [-1,1]:
                a=side*gap/2;b=side*(2+gap/2)
                ax.plot([a,b,b],[6,6,0],color=INK,lw=1)
            ax.set_title(f't = {oldt if row==0 else f["time_s"]:.2f} s',fontsize=10)
        axes[1,k].set_xlabel('x (m)')
    axes[0,0].set_ylabel('Original solver\nz (m)');axes[1,0].set_ylabel('Recreated solver\nz (m)')
    fig.suptitle('Same chamber hypothesis, wind and CO₂ scale — different numerical methods',fontsize=16,x=.47)
    fig.subplots_adjust(left=.065,right=.9,top=.87,bottom=.12,wspace=.15,hspace=.22)
    cax=fig.add_axes([.93,.2,.012,.57]);fig.colorbar(ScalarMappable(norm=norm,cmap=CO2),cax=cax,label='CO₂ (ppm)')
    fig.text(.065,.045,'Legacy snapshots use nearest saved time. Original solid-cell concentrations are shown as saved. New eddy closure and porous palm are additional assumptions.',fontsize=8)
    fig.savefig(out/'figures/legacy_comparison.png',dpi=120);plt.close(fig)


def render_campaign(out,root):
    preferred=next(name for name in ['refresh_r4','refresh_r3','refresh_r2','refresh_r1'] if (out/'data'/name/'summary.json').exists())
    for case,name in [(preferred,'cfd_refresh_3d'),('plume_r2','cfd_plume'),('cycles_r2','scenarios_timeseries')]:
        path=out/'data'/case
        if (path/'summary.json').exists():render_case(path,out,root,name)
    legacy_comparison(out,root,out/'data'/preferred)
    opening_detail(out)


def opening_detail(out):
    import subprocess
    subprocess.run(['ffmpeg','-v','error','-y','-t',str(31/6),'-i',str(out/'videos/cfd_refresh_3d.mp4'),
                    '-vf','setpts=2.5*PTS','-an','-c:v','libx264','-pix_fmt','yuv420p','-movflags','+faststart',
                    str(out/'videos/opening_detail.mp4')],check=True)
