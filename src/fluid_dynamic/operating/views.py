"""Dimension- and fan-count-aware panels; saved fields only."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from ..recreated.fan_forcing import FanArray
from ..recreated.render import series, style
from ..scenarios.geometry import ScenarioStudy
from ..scenarios import presentation_render as pr
WHITE='#FFFFFF'
INK='#202B30'
ON=GREEN='#007C66'
CO2=matplotlib.colormaps['viridis']
def white_style():
    style()
    plt.rcParams.update({'figure.facecolor':WHITE,'axes.facecolor':WHITE,
        'savefig.facecolor':WHITE,'savefig.edgecolor':WHITE,'savefig.transparent':False,
        'text.color':INK,'axes.labelcolor':INK,'font.size':10})


def load_cases(out):
    study = json.loads((out/'configs/study.json').read_text()); cases = {}
    for name, item in study['cases'].items():
        path = out/'data/cfd_runs'/name
        if not (path/'summary.json').exists():
            continue
        cases[name] = dict(name=name, path=path, cfg=ScenarioStudy(**item['config']),
            fans=FanArray(**item['fans']), ts=series(path/'timeseries.csv'),
            summary=json.loads((path/'summary.json').read_text()),
            label='Leaflets along fronds',
            color=GREEN)
    return study, cases


def marks(ax,case,gap,t,kind):
    fans=case['fans'];cfg=case['cfg']
    if not fans.enabled or not fans.heights_above_floor_m:return
    points=fans.positions(cfg,gap);on=fans.active(cfg,t);color=ON if on else '#666666'
    directed_marks(ax,case,points,color,kind)


def directed_marks(ax,case,points,color,kind):
    """Show actual per-fan directions; plan arrows project all configured heights."""
    fans=case['fans'];floor=case['cfg'].geometry['floor_z_m'];directions=fans.directions()
    x,y=points[0,:2];pole_top=points[:,2].max()+.12
    if kind=='3d':
        if fans.mounting_shell==0:
            ax.plot([x,x],[y,y],[floor,pole_top],color='#696F73',lw=2,zorder=7)
        ax.scatter(*points.T,s=30,color=color,edgecolor=WHITE,depthshade=False,zorder=8)
        for i,(p,a) in enumerate(zip(points,directions)):
            ax.quiver(*p,*(.8*a),color=color,arrow_length_ratio=.3,linewidth=2,zorder=8)
            end=p+.95*a
            ax.text(*end,f'F{i+1}',color=INK,fontsize=8,zorder=9)
    elif kind=='plan':
        ax.scatter(0,0,s=30,color='#27583B',zorder=8)
        ax.scatter(x,y,marker='s',s=35,color=color,edgecolor=WHITE,zorder=9)
        for i,(p,a) in enumerate(zip(points,directions)):
            end=p[:2]+1.05*a[:2]
            ax.annotate('',xy=end,xytext=p[:2],arrowprops=dict(arrowstyle='->',color=color,lw=2),zorder=9)
            label=p[:2]+1.3*a[:2]
            ax.text(*label,f'F{i+1}',ha='center',va='center',fontsize=8,color=INK,zorder=10,
                bbox=dict(facecolor=WHITE,edgecolor='none',alpha=.85,pad=.5))
        ax.text(.02,.02,f'Fan arrows: {len(points)} heights projected',transform=ax.transAxes,fontsize=7,
            bbox=dict(facecolor=WHITE,edgecolor='none',alpha=.85,pad=2))
    else:
        if fans.mounting_shell==0:
            ax.plot([y,y],[0,pole_top-floor],color='#696F73',lw=2,zorder=7)
        for i,(p,a) in enumerate(zip(points,directions)):
            z=p[2]-floor
            ax.scatter(p[1],z,marker='s',s=35,color=color,edgecolor=WHITE,zorder=8)
            if np.linalg.norm(a[1:])>1e-12:
                ax.annotate('',xy=(p[1]+.9*a[1],z+.9*a[2]),xytext=(p[1],z),
                    arrowprops=dict(arrowstyle='->',color=color,lw=2),zorder=9)
                label=f'F{i+1}'
            else:
                label=f'F{i+1} · +x (out of section)' if a[0]>0 else f'F{i+1} · −x (into section)'
            ax.text(p[1]+.12,z+.13,label,fontsize=7,color=INK,zorder=10,
                bbox=dict(facecolor=WHITE,edgecolor='none',alpha=.8,pad=.5))


def plan(ax,d,case,norm,t):
    cfg=case['cfg'];floor=cfg.geometry['floor_z_m']
    iz=int(np.argmin(abs(d['z']-(floor+cfg.height/2))))
    ax.pcolormesh(d['x_faces'],d['y_faces'],d['co2'][:,:,iz].T,cmap=CO2,norm=norm,shading='flat')
    x,y=d['x'],d['y'];qx,qy=np.meshgrid(np.arange(0,len(x),2),np.arange(0,len(y),2),indexing='ij')
    u=d['velocity'][:,:,iz,:]
    ax.quiver(x[qx],y[qy],u[qx,qy,0],u[qx,qy,1],color=INK,scale_units='xy',angles='xy',scale=2,width=.0025)
    pr.plan_panels(ax,cfg,float(d['gap']),float(d['z'][iz]))
    marks(ax,case,float(d['gap']),t,'plan')
    ex,ey,_=pr.extents(cfg)
    ax.set(xlim=(-ex,ex),ylim=(-ey,ey),xlabel='x along slide (m)',ylabel='y (m)',
        title=f'{case["label"]} · top view at {d["z"][iz]-floor:.2f} m above floor')
    ax.set_aspect('equal',adjustable='box')


def section(ax,d,case,norm,t):
    cfg=case['cfg'];floor=cfg.geometry['floor_z_m'];gap=float(d['gap'])
    xpos=case['fans'].closed_x_m;ix=int(np.argmin(abs(d['x']-xpos)))
    ax.pcolormesh(d['y_faces'],d['z_faces']-floor,d['co2'][ix,:,:].T,cmap=CO2,norm=norm,shading='flat')
    y,z=d['y'],d['z']-floor
    qy,qz=np.meshgrid(np.arange(len(y)),np.arange(len(z)),indexing='ij');u=d['velocity'][ix,:,:,:]
    ax.quiver(y[qy],z[qz],u[qy,qz,1],u[qy,qz,2],color=INK,scale_units='xy',angles='xy',scale=2,width=.002)
    ax.plot([-cfg.depth/2,cfg.depth/2,cfg.depth/2,-cfg.depth/2,-cfg.depth/2],
        [0,0,cfg.height,cfg.height,0],color=INK,lw=1.4)
    marks(ax,case,gap,t,'section')
    ax.set(xlim=(-cfg.depth/2-1,cfg.depth/2+1),ylim=(-floor,cfg.height+.65),
        xlabel='y (m); fan arrows projected' if case['fans'].axes is not None else 'y in blowing direction (m)',ylabel='Height above floor (m)',
        title=f'Fan-height section · saved x={d["x"][ix]:.2f} m')
    ax.set_aspect('equal',adjustable='box')


def history(ax, cases, t, key):
    for c in cases:
        ax.plot(c['ts']['time_s'], c['ts'][key], color=c['color'], lw=1.6, label=c['label'])
    cfg = cases[0]['cfg']; ax.axvline(t, color=INK, lw=.8)
    if key == 'roi_mean_ppm':
        ax.axhline(cfg.ambient_ppm, color='#999999', ls=':', lw=1)
        ax.set(ylabel='Mean CO₂ (ppm)', title=f'Mean across the fixed {cfg.width*cfg.depth*cfg.height:g} m³ region')
    else:
        ax.set(ylabel='Spatial standard deviation (ppm)', title='Variation between cells; lower means more uniform')
    ax.set(xlim=(0, cfg.duration_s), xlabel='Physical time (s)'); ax.grid(alpha=.15)
    ax.legend(frameon=False, fontsize=8, loc='best')


def operating(axes, c, t):
    ts = c['ts']; cfg = c['cfg']
    for ax, key in zip(axes[:2], ['roi_mean_ppm', 'roi_std_ppm']):
        history(ax, [c], t, key)
    axes[2].plot(ts['time_s'], ts['source_umol_s'], color=GREEN, label='Expanded leaflet total')
    axes[2].set(ylabel='Net exchange (µmol/s)', title='Uptake stays active during opening and while open')
    axes[3].plot(ts['time_s'], ts['gap_m']/cfg.gap, color='#385FC3', label='Opening fraction')
    axes[3].step(ts['time_s'], ts['fans_on'], where='post', color=GREEN, label='Fans on = 1')
    axes[3].set(ylabel='Operating state', title=f'{len(c["fans"].heights_above_floor_m)} fans · closed only')
    for ax in axes[2:]:
        ax.axvline(t, color=INK, lw=.8); ax.grid(alpha=.15)
        ax.set(xlim=(0, cfg.duration_s), xlabel='Physical time (s)'); ax.legend(frameon=False, fontsize=8)


def source_map(ax, case, frame_index, g, norm):
    frame = case['summary']['frames'][frame_index]
    with np.load(case['path']/frame['path']) as d:
        cfg = case['cfg']; gap = float(d['gap'])
        with np.load(case['path']/'leaf_area'/Path(frame['path']).name) as source:
            q = -source['source_umol_s']
        # Depth-integrated uptake density, not concentration or a leaf rate.
        area = np.diff(d['x_faces'])[:, None]*np.diff(d['z_faces'])[None, :]
        density = q.sum(axis=1)/area
        zf = d['z_faces']-cfg.geometry['floor_z_m']
        ax.pcolormesh(d['x_faces'], zf, density.T, cmap=CO2, norm=norm, shading='flat')
        projection = g['axes'][:, :, [0, 2]].copy(); projection[..., 1] -= cfg.geometry['floor_z_m']
        ax.add_collection(LineCollection(projection, colors='white', linewidths=.25, alpha=.5))
        ax.set(xlim=(-cfg.width/2, cfg.width/2), ylim=(0, cfg.height), xlabel='x (m)', ylabel='Height above floor (m)', title=case['label'])
        ax.set_aspect('equal')


