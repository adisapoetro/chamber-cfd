"""Preferred fan-video visual style; source surfaces are explained separately.

Movie layout derives from the preserved renderer. No shared renderer or earlier
bundle is changed. The main scene uses only the exact original display mesh.
"""
from contextlib import nullcontext
import json
from pathlib import Path
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from .views import white_style, marks, plan, section, WHITE, INK, CO2
from .views import history, operating, source_map
from .workflow import digest, write_json
from .views import load_cases as existing_load_cases
from ..scenarios import presentation_render as pr
GREEN = '#087E69'
ORANGE = '#AF661E'


def load_cases(out):
    study, cases = existing_load_cases(out)
    for name, c in cases.items():
        c['label'] = 'Leaves along fronds'
    return study, cases


def renderer_inputs(out, root, cases, kind):
    paths = [Path(__file__), root/'src/fluid_dynamic/operating/views.py',
        root/'src/fluid_dynamic/recreated/render.py', root/'src/fluid_dynamic/scenarios/presentation_render.py']
    return dict(kind=kind, background=WHITE, geometry_sha256=digest(out/'geometry/palm.npz'),
        study_sha256=digest(out/'configs/study.json'),
        renderer_hashes={str(p.relative_to(root)): digest(p) for p in paths},
        cases={c['name']: {name: digest(c['path']/name) for name in ('config.json', 'fan_config.json', 'summary.json', 'timeseries.csv')}
            for c in cases})


def movie(out, root, study, cases, g, name, kind, *, preview=False):
    if [c['name'] for c in cases] != ['main']:
        raise ValueError('Presentation must contain only the C2 leaflet-uptake case')
    white_style(); cfg = cases[0]['cfg']; frames = pr.synchronized_frames(cases)
    preferred = study['display']['preferred_limits_ppm']
    low = min(c['summary']['minimum_ppm'] for c in cases)
    high = max(c['summary']['maximum_ppm'] for c in cases)
    limits = [min(preferred[0], float(np.floor(low))), max(preferred[1], float(np.ceil(high)))]
    norm = Normalize(*limits)
    record = out/'verification/media'/(name+'.json')
    inputs = renderer_inputs(out, root, cases, kind)
    if record.exists() and not preview:
        old = json.loads(record.read_text())
        if old['inputs'] != inputs or any(digest(out/p) != sha for p, sha in old['output_hashes'].items()):
            raise ValueError('Changed rendering input or output requires a separate result folder')
        return old
    source_norm = Normalize(0, 9.)
    # Fixed source scale includes every saved depth-integrated source value.
    if kind == 'source':
        maximum = 0.; minimum = 0.
        for c in cases:
            for frame in c['summary']['frames']:
                with np.load(c['path']/frame['path']) as d:
                    with np.load(c['path']/'leaf_area'/Path(frame['path']).name) as s:
                        q = -s['source_umol_s']
                    density = q.sum(axis=1)/(np.diff(d['x_faces'])[:, None]*np.diff(d['z_faces'])[None, :])
                    maximum = max(maximum, float(density.max())); minimum = min(minimum, float(density.min()))
        source_norm = Normalize(np.floor(minimum), max(1., np.ceil(maximum)))
    fig = plt.figure(figsize=(16, 9), facecolor=WHITE)
    gs = fig.add_gridspec(2, 2, left=.065, right=.9, bottom=.16, top=.81,
        height_ratios=[1.7, 1], hspace=.4, wspace=.25)
    axes = [fig.add_subplot(gs[0, 0], projection='3d', computed_zorder=False) if kind in ('3d', 'source') else fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    titles = {'3d': study['title'], 'section': 'Fan circulation and leaf uptake',
        'operating': 'Leaf uptake and chamber operation', 'source': 'Where the CO₂ exchange is applied'}
    cy=study['inputs']['cycle']
    fig.text(.045, .95, 'WHOLE-TREE CHAMBER  /  OPERATING DESIGN', fontsize=10, weight='bold', color=GREEN)
    fig.text(.045, .905, titles[kind], fontsize=22, weight='bold')
    fig.text(.045, .862, f'{cfg.width:g} × {cfg.depth:g} × {cfg.height:g} m airspace · '
        f'{cy["closed_s"]/60:g} min closed / {cy["open_phase_s"]/60:g} min open phase · {cfg.duration_s:g} s', fontsize=11)
    fig.text(.955, .95, 'DIAGNOSTIC · FAN PERFORMANCE UNVERIFIED', ha='right', fontsize=9, color=ORANGE)
    stamp = fig.text(.955, .858, '', ha='right', fontsize=10)
    if kind != 'operating':
        cax = fig.add_axes([.932, .42, .011, .34]); cax.set_facecolor(WHITE)
        color_norm = source_norm if kind == 'source' else norm
        label = 'CO₂ removal integrated through y (µmol/m²/s)' if kind == 'source' else 'CO₂ (ppm), fixed linear viridis'
        fig.colorbar(ScalarMappable(norm=color_norm, cmap=CO2), cax=cax, label=label,
            ticks=np.linspace(color_norm.vmin, color_norm.vmax, 5))
    rate_label=f'{cfg.tree_source_umol_s:.2f}'.replace('-', '−')
    fig.text(.045, .081, f'{study["inputs"]["exchange"]["label"]} · net exchange {rate_label} µmol/s throughout\n'
        'Representative leaflets along fronds · trunk and bare bases excluded · schematic palm', fontsize=10)
    f=cases[0]['fans']
    fig.text(.045,.029,f'{len(f.heights_above_floor_m)} fans · fixed pole x={f.closed_x_m:g}, y={f.closed_y_m:g} m · '
        'prescribed force and source geometry · saved fields',fontsize=9,color='#48545A')
    path = out/'videos'/(name+'.mp4'); tmp = path.with_name(name+'.encoding.mp4')
    writer = FFMpegWriter(fps=study["display"]["fps"], codec='libx264', extra_args=['-pix_fmt', 'yuv420p', '-crf', '20', '-movflags', '+faststart', '-threads', '1'])
    count = 0; mapping = []; figures = []
    desired = [0, cfg.open_at, cfg.open_at+cfg.opening_s/2, cfg.period_s, cfg.duration_s]
    still_indices = {int(np.argmin([abs(f['time_s']-t) for f in frames[0]])) for t in desired}
    jpreview=min(still_indices, key=lambda j: abs(frames[0][j]['time_s']-(cfg.open_at+cfg.opening_s/2)))
    chosen = [(jpreview, frames[0][jpreview])] if preview else enumerate(frames[0])
    with (nullcontext() if preview else writer.saving(fig, str(tmp), dpi=100)):
        for j, frame in chosen:
            t = frame['time_s']
            with np.load(cases[0]['path']/frame['path']) as d:
                gap = float(d['gap'])
                for ax in axes:
                    ax.clear(); ax.set_facecolor(WHITE)
                if kind == 'operating':
                    operating(axes, cases[0], t)
                elif kind == 'source':
                    leaflet_surfaces(axes[0], g, cfg)
                    source_map(axes[1], cases[0], j, g, source_norm)
                    history(axes[2], cases, t, 'roi_mean_ppm'); history(axes[3], cases, t, 'roi_std_ppm')
                else:
                    if kind == '3d':
                        pr.draw_scene(axes[0], d, cfg, norm, list(g['display_triangles']))
                        marks(axes[0], cases[0], gap, t, '3d')
                        for axis in (axes[0].xaxis, axes[0].yaxis, axes[0].zaxis):
                            axis.set_pane_color((1, 1, 1, 1)); axis.pane.fill = False
                        axes[0].set_title('Same simple palm · saved CO₂ slices and velocity', fontsize=10, pad=-7)
                    else:
                        section(axes[0], d, cases[0], norm, t)
                    plan(axes[1], d, cases[0], norm, t)
                    history(axes[2], cases, t, 'roi_mean_ppm'); history(axes[3], cases, t, 'roi_std_ppm')
                state = 'SEALED' if gap < 1e-8 else 'FULLY OPEN' if gap > cfg.gap-1e-8 else 'MOVING'
                if abs(t % cfg.period_s-cfg.open_at) < 1e-8:
                    state = 'OPENING START'
                stamp.set_text(f'{t:.0f} s · {state} · fans {"ON" if cases[0]["fans"].active(cfg, t) else "OFF"}')
                n = pr.repeats(gap, cfg.gap, j)
                if preview or j in still_indices:
                    p = out/'figures'/f'{name}_{"preview_" if preview else ""}{round(t):04d}s.png'
                    fig.savefig(p, dpi=100, facecolor=WHITE, transparent=False); figures.append(str(p.relative_to(out)))
                mapping.append(dict(physical_time_s=t, saved_frames=[fs[j]['path'] for fs in frames],
                    first_encoded_frame=count, repeats=n)); count += n
                if not preview:
                    for _ in range(n):
                        writer.grab_frame(facecolor=WHITE, transparent=False)
    plt.close(fig)
    if preview:
        return dict(preview=figures)
    media = pr.finish_movie(tmp, path, count); gif = pr.gif(path)
    media.update(background=WHITE, color_map='viridis', color_limits_ppm=limits if kind in ('3d', 'section') else None,
        source_scale_umol_m2_s=[source_norm.vmin, source_norm.vmax] if kind == 'source' else None,
        source_cases=[c['name'] for c in cases], frame_mapping=mapping,
        data_treatment='Saved fields and physical times; motion frames repeated; no interpolation')
    files = [str(path.relative_to(out)), str(gif.relative_to(out)), *figures]
    result = dict(inputs=inputs, media=media, output_hashes={p: digest(out/p) for p in files})
    write_json(record, result); print('Rendered', name, flush=True)
    return result


def leaflet_surfaces(ax, g, cfg):
    """Explain the single case's uptake geometry without a control panel."""
    floor = cfg.geometry['floor_z_m']
    ax.add_collection3d(Poly3DCollection(g['display_triangles'], facecolor='#27583B', edgecolor='none', alpha=.22))
    ax.add_collection3d(Poly3DCollection(g['triangles'][g['source_mask']], facecolor='#1EA885', edgecolor='#166347', linewidths=.15, zorder=4))
    ax.set(xlim=(-cfg.width/2, cfg.width/2), ylim=(-cfg.depth/2, cfg.depth/2), zlim=(floor, floor+cfg.height), xlabel='x (m)', ylabel='y (m)', zlabel='z (m)')
    ax.view_init(elev=22, azim=-64); ax.set_box_aspect((cfg.width, cfg.depth, cfg.height)); ax.grid(False)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((1, 1, 1, 1)); axis.pane.fill=False
    ax.set_title('Representative leaf surfaces used for uptake', fontsize=10, pad=-7)


def static_figures(out, cases, g):
    white_style(); selected = [cases['main']]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), facecolor=WHITE)
    operating(list(axes.flat), selected[0], 0)
    fig.suptitle('Continuous leaf uptake and chamber operation')
    fig.tight_layout(); fig.savefig(out/'figures/operating_schedule.png', dpi=130, facecolor=WHITE); plt.close(fig)
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), facecolor=WHITE)
    for ax, key in zip(axes, ['roi_mean_ppm', 'roi_std_ppm']):
        history(ax, selected, 0, key)
    fig.tight_layout(); fig.savefig(out/'figures/co2_timeseries.png', dpi=130, facecolor=WHITE); plt.close(fig)
    fig = plt.figure(figsize=(12, 6), facecolor=WHITE)
    cfg = selected[0]['cfg']; floor = cfg.geometry['floor_z_m']
    for i in (1, 2):
        ax = fig.add_subplot(1, 2, i, projection='3d', computed_zorder=False)
        ax.add_collection3d(Poly3DCollection(g['display_triangles'], facecolor='#27583B', edgecolor='none', alpha=1 if i == 1 else .22))
        if i == 2:
            ax.add_collection3d(Poly3DCollection(g['triangles'][g['source_mask']], facecolor='#1EA885', edgecolor='#166347', linewidths=.15, zorder=4))
        ax.set(xlim=(-cfg.width/2, cfg.width/2), ylim=(-cfg.depth/2, cfg.depth/2), zlim=(floor, floor+cfg.height), xlabel='x (m)', ylabel='y (m)', zlabel='z (m)')
        ax.view_init(elev=22, azim=-64); ax.set_box_aspect((cfg.width, cfg.depth, cfg.height)); ax.grid(False)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.set_pane_color((1, 1, 1, 1)); axis.pane.fill=False
        ax.set_title('Unchanged palm in the main video' if i == 1 else 'Representative leaf surfaces used for uptake', fontsize=11)
    fig.suptitle('Same simple palm · uptake belongs to leaves along the fronds', fontsize=15)
    fig.text(.5, .035, 'Right: schematic source patches, shown only to explain the calculation; no trunk or bare-frond-base sink.', ha='center', fontsize=10)
    fig.savefig(out/'figures/where_uptake_happens.png', dpi=140, facecolor=WHITE); plt.close(fig)
    if selected[0]['fans'].mounting_shell == 0:
        fan_configuration(out, selected[0])


def fan_configuration(out, case):
    """Mount and direction schematic, separately labelled from saved airflow."""
    cfg=case['cfg'];fans=case['fans'];points=fans.positions(cfg,0)
    fig,axes=plt.subplots(1,2,figsize=(12,6),facecolor=WHITE)
    ax=axes[0]
    ax.plot(np.array([-1,1,1,-1,-1])*cfg.width/2,np.array([-1,-1,1,1,-1])*cfg.depth/2,color=INK,lw=1.2)
    ax.add_patch(plt.Circle((0,0),cfg.plant['trunk_radius_m'],facecolor='#27583B',label='Palm trunk'))
    marks(ax,case,0,0,'plan')
    ax.set(xlim=(-cfg.width/2-.3,cfg.width/2+.3),ylim=(-cfg.depth/2-.3,cfg.depth/2+.3),xlabel='x along slide (m)',ylabel='y (m)',
        title='Top view · prescribed discharge directions')
    ax.set_aspect('equal');ax.legend(loc='upper right',fontsize=8)
    ax=axes[1]
    ax.plot([0,0],[0,max(fans.heights_above_floor_m,default=0)+.12],color='#696F73',lw=3)
    for i,(height,axis) in enumerate(zip(fans.heights_above_floor_m,fans.directions())):
        angle=np.degrees(np.arctan2(axis[1],axis[0]))%360
        ax.scatter(0,height,marker='s',s=60,color=GREEN,zorder=3)
        ax.text(.16,height,f'F{i+1} · {height:.2f} m · {angle:.0f}°',va='center',fontsize=11)
    ax.axhline(0,color=INK,lw=1)
    ax.set(xlim=(-.25,1.5),ylim=(-.1,cfg.height),xticks=[],ylabel='Height above chamber floor (m)',
        title=f'{len(points)} fans · prescribed height above floor')
    ax.text(.04,.93,f'Pole: x={fans.closed_x_m:.2f} m, y={fans.closed_y_m:.2f} m',transform=ax.transAxes,fontsize=10)
    fig.suptitle('Fixed support beside the trunk · closed-only operation',fontsize=16)
    fig.text(.5,.045,'Position and absolute azimuth assumed · pole blockage unresolved · nominal fan performance unverified',
        ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.085,1,.94))
    fig.savefig(out/'figures/fan_configuration.png',dpi=140,facecolor=WHITE);plt.close(fig)


def render_bundle(out, root, *, preview=False):
    study, cases = load_cases(out)
    pr.FPS = study['display']['fps']
    if set(cases) != set(study['cases']):
        raise ValueError('Complete every simulation before rendering')
    with np.load(out/'geometry/palm.npz') as d:
        g = {k: d[k] for k in d.files}
    static_figures(out, cases, g)
    selected = [cases['main']]
    for name, kind in [('cycle', '3d'), ('fan_sections', 'section'),
                       ('operating_timeseries', 'operating'), ('leaf_source_location', 'source')]:
        movie(out, root, study, selected, g, name, kind, preview=preview)
    if not preview:
        for p in renderer_inputs(out, root, selected, 'all')['renderer_hashes']:
            dest = out/'configs/source_rendered'/p; dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(root/p, dest)
