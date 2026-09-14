#!/usr/bin/env python3
"""Render saved WaterLily velocities. This view contains no simulated CO₂."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import animation, colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np


def main(folder):
    folder = folder.resolve()
    summary = json.loads((folder / 'summary.json').read_text())
    config = json.loads((folder / 'inputs.json').read_text())
    assert summary['engine'] == 'WaterLily' and summary['co2_transport_enabled'] is False
    shape = tuple(summary['shape'])
    lower = np.array(summary['grid_lower_m'])
    dx = summary['grid_cell_m']
    xyz = [lower[i] + (np.arange(shape[i]) + .5) * dx for i in range(3)]
    ch = config['chamber']
    plant = json.loads((Path(__file__).resolve().parents[1] / 'assets/reference_palm.json').read_text())
    triangles = np.array(plant['display_triangles'])
    scale = np.array([config['palm']['crown_radius_x_m'] / 1.52,
                      config['palm']['crown_radius_y_m'] / 1.52,
                      config['palm']['height_m'] / 2.39])
    triangles = (triangles - [0, 0, .6]) * scale + [0, 0, ch['floor_m']]
    z_index = int(np.argmin(abs(xyz[2] - (ch['floor_m'] + .5 * ch['height_m']))))
    records = summary['snapshots']
    time = np.array([r['time_s'] for r in records])
    vmax = max(r['max_speed_m_s'] for r in records)
    norm = colors.Normalize(0, np.ceil(vmax * 10) / 10)
    for sub in ['videos', 'figures']:
        (folder / sub).mkdir(exist_ok=True)
    mp4 = folder / 'videos/airflow_cycle.mp4'
    gif = folder / 'videos/airflow_cycle.gif'
    if mp4.exists() or gif.exists():
        raise FileExistsError('Airflow media already exist; preserve them before rendering again.')
    plt.rcParams.update({'figure.facecolor': 'white', 'axes.facecolor': 'white',
                         'savefig.facecolor': 'white', 'font.family': 'DejaVu Sans',
                         'text.color': '#233334', 'axes.labelcolor': '#233334',
                         'xtick.color': '#233334', 'ytick.color': '#233334'})
    fig = plt.figure(figsize=(14, 8), facecolor='white')
    grid = fig.add_gridspec(2, 2, left=.065, right=.88, top=.82, bottom=.13,
                           height_ratios=[2.2, 1], hspace=.4, wspace=.3)
    ax3 = fig.add_subplot(grid[0, 0], projection='3d')
    plan = fig.add_subplot(grid[0, 1])
    speed_ax = fig.add_subplot(grid[1, 0])
    state_ax = fig.add_subplot(grid[1, 1])
    cax = fig.add_axes([.915, .47, .015, .32])
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='viridis'), cax=cax,
                 label='Air speed (m/s) · fixed scale')
    fig.text(.065, .965, 'WHOLE-TREE CHAMBER  /  WATERLILY 1.8.0', color='#007F71', weight='bold', size=10)
    fig.text(.065, .915, 'Central fans and moving chamber walls', weight='bold', size=23)
    fig.text(.065, .876, f"{ch['width_m']:g} × {ch['depth_m']:g} × {ch['height_m']:g} m nominal airspace · "
             f"{config['cycle']['closed_s']/60:g} min closed / {config['cycle']['open_phase_s']/60:g} min open phase", size=11)
    stamp = fig.text(.94, .876, '', ha='right', size=10)
    fig.text(.94, .963, 'AIRFLOW DIAGNOSTIC · CO₂ TRANSPORT NOT ENABLED', ha='right', color='#AD650C', size=8)
    fig.text(.065, .068, 'Native WaterLily velocity fields · schematic palm · fans operate only while sealed', size=10)
    fig.text(.065, .039, f"Numerical walls: {summary['numerical_wall_thickness_m']:g} m thick, extended outward · "
             'geometry and grid convergence unestablished', size=9, color='#5B6565')
    fps = config['display']['fps']
    writer = animation.FFMpegWriter(fps=fps, codec='libx264',
                                    extra_args=['-pix_fmt', 'yuv420p', '-crf', '20'])
    still_times = [0, config['cycle']['closed_s'],
                   config['cycle']['closed_s'] + config['cycle']['opening_travel_s']/2,
                   config['cycle']['closed_s'] + config['cycle']['open_phase_s'], summary['duration_s']]
    still_indices = {int(np.argmin(abs(time - t))) for t in still_times}
    with writer.saving(fig, str(mp4), dpi=100):
        for index, record in enumerate(records):
            velocity = np.fromfile(folder / record['velocity_file'], dtype='<f8').reshape(shape + (3,), order='F')
            distance = np.fromfile(folder / record['distance_file'], dtype='<f8').reshape(shape, order='F')
            speed = np.linalg.norm(velocity, axis=-1)
            gap = record['gap_m']
            for axis in [ax3, plan, speed_ax, state_ax]:
                axis.clear()
            for side in [-1, 1]:
                x0, x1 = sorted([side*gap/2, side*(gap/2+ch['width_m']/2)])
                y0, y1 = -ch['depth_m']/2, ch['depth_m']/2
                z0, z1 = ch['floor_m'], ch['floor_m']+ch['height_m']
                for z in [z0, z1]:
                    ax3.plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],[z]*5,color='#526B6C',lw=.8)
                for x in [x0,x1]:
                    for y in [y0,y1]:ax3.plot([x,x],[y,y],[z0,z1],color='#526B6C',lw=.8)
                inner=side*gap/2;outer=side*(gap/2+ch['width_m']/2)
                plan.plot([inner,outer,outer,inner],[y1,y1,y0,y0],color='#263C3D',lw=1.4)
            ax3.add_collection3d(Poly3DCollection(triangles,facecolors='#235238',alpha=.85,edgecolors='none'))
            X,Y=np.meshgrid(xyz[0],xyz[1],indexing='ij')
            plane=np.ma.masked_where(distance[:,:,z_index]<0,speed[:,:,z_index])
            cmap=plt.colormaps['viridis'].copy();cmap.set_bad('#F0F0F0')
            plan.pcolormesh(X,Y,plane,cmap=cmap,norm=norm,shading='nearest',rasterized=True)
            select=(slice(None,None,2),slice(None,None,2))
            u=np.ma.masked_where(distance[:,:,z_index]<0,velocity[:,:,z_index,0])
            v=np.ma.masked_where(distance[:,:,z_index]<0,velocity[:,:,z_index,1])
            plan.quiver(X[select],Y[select],u[select],v[select],color='#183334',scale=25,width=.002)
            for fan in summary['fan_configuration']:
                p=np.array(fan['position']);d=np.array(fan['direction'])
                ax3.quiver(*p,*d,length=.5,color='#008B78' if record['fans_on'] else '#777777')
                ax3.scatter(*p,color='#008B78',s=16)
                plan.arrow(*p[:2],*(.65*d[:2]),head_width=.12,color='#DBECE5',length_includes_head=True)
            ax3.set(xlim=(-4,4),ylim=(-4,4),zlim=(0,ch['floor_m']+ch['height_m']+.3),
                    xlabel='x (m)',ylabel='y (m)',zlabel='z (m)',title='Schematic palm and nominal inside surfaces')
            ax3.set_box_aspect((8,8,5));ax3.view_init(elev=22,azim=-64);ax3.grid(False)
            for axis in [ax3.xaxis,ax3.yaxis,ax3.zaxis]:axis.set_pane_color((1,1,1,1))
            plan.set(xlim=(-5,5),ylim=(-5,5),aspect='equal',xlabel='x along slide (m)',ylabel='y (m)',
                     title=f'Air speed at {xyz[2][z_index]-ch["floor_m"]:.2f} m above floor')
            speed_ax.plot(time,[r['roi_mean_speed_m_s'] for r in records],color='#008B78')
            speed_ax.set(ylabel='Mean speed (m/s)',title='Fixed 96 m³ reference region' if np.isclose(summary['nominal_reference_volume_m3'],96)
                         else f"Fixed {summary['nominal_reference_volume_m3']:g} m³ reference region")
            state_ax.plot(time,[r['gap_m']/ch['maximum_gap_m'] for r in records],label='Opening fraction',color='#3D60C8')
            state_ax.step(time,[r['fans_on'] for r in records],where='post',label='Fans on = 1',color='#008B78')
            state_ax.set(ylabel='Operating state',ylim=(-.05,1.1),title='Opening and fan schedule')
            state_ax.legend(frameon=False,loc='center right',fontsize=8)
            for axis in [speed_ax,state_ax]:
                axis.axvline(record['time_s'],color='#385555',lw=.8)
                axis.set(xlim=(0,max(time[-1],1)),xlabel='Physical time (s)');axis.grid(alpha=.15)
                axis.spines[['top','right']].set_visible(False)
            stamp.set_text(f"{record['time_s']:g} s · fans {'ON' if record['fans_on'] else 'OFF'}")
            if index in still_indices:
                fig.savefig(folder/f"figures/airflow_{record['time_s']:06.1f}s.png",dpi=140,facecolor='white')
            writer.grab_frame(facecolor='white')
    plt.close(fig)
    subprocess.run(['ffmpeg','-v','error','-y','-i',str(mp4),'-filter_complex',
                    f'fps={fps},scale=960:-1:flags=lanczos,split[a][b];[a]palettegen[p];[b][p]paletteuse',
                    '-loop','0',str(gif)],check=True)
    media = sorted([*folder.glob('videos/*'),*folder.glob('figures/*.png')])
    for path in media:
        subprocess.run(['ffmpeg','-v','error','-xerror','-i',str(path),'-f','null','-'],check=True,stdout=subprocess.DEVNULL)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (folder/'rendering.json').write_text(json.dumps({'renderer_sha256':sha(Path(__file__)),
        'media_decode':'PASS','media':[str(p.relative_to(folder)) for p in media],
        'co2_transport_enabled':False},indent=2)+'\n')
    files={str(p.relative_to(folder)):sha(p) for p in sorted(folder.rglob('*')) if p.is_file() and p.name!='manifest.json'}
    (folder/'manifest.json').write_text(json.dumps({'sha256':files},indent=2)+'\n')
    print(json.dumps({'media_decode':'PASS','media_files':len(media),'output':str(folder)},indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path)
    main(parser.parse_args().folder)
