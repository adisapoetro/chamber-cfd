#!/usr/bin/env python3
"""White-background 3D and sectional views of saved WaterLily + Julia CO₂ fields."""
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

VIEWS = ['c2_202502_cycle', 'chamber_3d', 'fan_sections',
         'operating_timeseries', 'leaf_source_location']


class Result:
    def __init__(self, folder):
        self.folder = folder.resolve()
        self.summary = json.loads((folder / 'summary.json').read_text())
        self.config = json.loads((folder / 'inputs.json').read_text())
        assert self.summary['engine'] == 'WaterLily' and self.summary['co2_transport_enabled']
        self.records = self.summary['snapshots']
        self.time = np.array([r['time_s'] for r in self.records])
        self.shape = tuple(self.summary['shape'])
        self.dx = self.summary['grid_cell_m']
        self.lower = np.array(self.summary['grid_lower_m'])
        self.xyz = [self.lower[k] + (np.arange(self.shape[k]) + .5)*self.dx for k in range(3)]
        self.edges = [self.lower[k] + np.arange(self.shape[k]+1)*self.dx for k in range(3)]
        self.ch = self.config['chamber']
        self.X, self.Y = np.meshgrid(*self.xyz[:2], indexing='ij')
        self.area = self.array('fields/leaflet_area.f64')
        frozen = folder / 'configs/source/assets/reference_palm.json'
        plant = json.loads(frozen.read_text())
        env = plant['reference_envelope']
        p = self.config['palm']
        scale = [p['crown_radius_x_m']/env['crown_radii_m'][0],
                 p['crown_radius_y_m']/env['crown_radii_m'][1], p['height_m']/env['height_m']]
        transform = lambda a: (np.array(a)-[0,0,env['floor_m']])*scale+[0,0,self.ch['floor_m']]
        self.palm = transform(plant['display_triangles'])
        self.leaves = transform(plant['triangles'])[np.array(plant['source_mask'])]
        self.z = int(np.argmin(abs(self.xyz[2]-(self.ch['floor_m']+.8*p['height_m']))))
        lo = min(self.config['display']['co2_min_ppm'], np.floor(min(r['co2_min_ppm'] for r in self.records)))
        hi = max(self.config['display']['co2_max_ppm'], np.ceil(max(r['co2_max_ppm'] for r in self.records)))
        self.norm = colors.Normalize(lo, hi)
        self.cmap = plt.colormaps['viridis'].copy()
        self.cmap.set_bad('white')

    def array(self, name, velocity=False):
        shape = self.shape + ((3,) if velocity else ())
        return np.fromfile(self.folder/name, dtype='<f8').reshape(shape, order='F')

    def load(self, r):
        self.C = self.array(r['co2_file'])
        self.volume = self.array(r['volume_file'])
        self.u = self.array(r['velocity_file'], True)

    def walls(self, ax, gap, three=False):
        ch = self.ch
        for side in [-1, 1]:
            x0,x1 = sorted([side*gap/2, side*(gap/2+ch['width_m']/2)])
            y0,y1 = -ch['depth_m']/2, ch['depth_m']/2
            z0,z1 = ch['floor_m'], ch['floor_m']+ch['height_m']
            if three:
                for z in [z0,z1]: ax.plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],[z]*5,color='#526B6C',lw=.8)
                for x in [x0,x1]:
                    for y in [y0,y1]: ax.plot([x,x],[y,y],[z0,z1],color='#526B6C',lw=.8)
            else:
                inner=side*gap/2;outer=side*(gap/2+ch['width_m']/2)
                ax.plot([inner,outer,outer,inner],[y1,y1,y0,y0],color='#263C3D',lw=1.3)

    def plan(self, ax, r, z, source=False):
        if source:
            values = self.area[:,:,z] / self.area.sum()*abs(self.config['exchange']['net_co2_umol_s'])
            values = np.ma.masked_where(values==0, values)
            ax.pcolormesh(self.X,self.Y,values,cmap='Greens',shading='nearest',vmin=0,
                          vmax=max(self.area.max()/self.area.sum()*abs(self.config['exchange']['net_co2_umol_s']),1e-12))
        else:
            mask = self.volume[:,:,z]<=0
            ax.pcolormesh(self.X,self.Y,np.ma.masked_where(mask,self.C[:,:,z]),cmap=self.cmap,
                          norm=colors.Normalize(self.norm.vmin,self.norm.vmax),shading='nearest')
            pick=(slice(None,None,2),slice(None,None,2))
            ax.quiver(self.X[pick],self.Y[pick],np.ma.masked_where(mask,self.u[:,:,z,0])[pick],
                      np.ma.masked_where(mask,self.u[:,:,z,1])[pick],color='#173A37',scale=24,width=.002)
        self.walls(ax,r['gap_m'])
        for i,f in enumerate(self.summary['fan_configuration']):
            p=np.array(f['position']);d=np.array(f['direction'])
            if int(np.argmin(abs(self.xyz[2]-p[2]))) != z:
                continue
            colour='#D88426' if r['fans_on'] else '#92979D'
            ax.arrow(*p[:2],*(.5*d[:2]),head_width=.12,color=colour,length_includes_head=True)
            ax.text(p[0]+.12,p[1]-.25,f'F{i+1}',color=colour,fontsize=7)
        xlim=self.ch['width_m']/2+self.ch['maximum_gap_m']/2+1
        ylim=self.ch['depth_m']/2+1
        ax.set(xlim=(-xlim,xlim),ylim=(-ylim,ylim),aspect='equal',xlabel='x along slide (m)',ylabel='y (m)')
        ax.set_title(f"{'Leaflet source' if source else 'CO₂ + native airflow'} · {self.xyz[2][z]-self.ch['floor_m']:.2f} m above floor",fontsize=10)

    def concentration_sections(self, ax):
        """Plot cell colours on their actual faces; no invented volume interpolation."""
        indices = [int(np.argmin(abs(self.xyz[k]))) for k in (0, 1)] + [self.z]
        ch = self.ch
        limits = [(-(ch['width_m']+ch['maximum_gap_m']+1)/2, (ch['width_m']+ch['maximum_gap_m']+1)/2),
                  (-(ch['depth_m']+1)/2, (ch['depth_m']+1)/2), (0, ch['floor_m']+ch['height_m'])]
        selected = [np.flatnonzero((x >= lo) & (x <= hi)) for x, (lo, hi) in zip(self.xyz, limits)]
        for fixed, index in enumerate(indices):
            axes = [k for k in range(3) if k != fixed]
            edge = [self.edges[k][selected[k][0]:selected[k][-1]+2] for k in axes]
            a, b = np.meshgrid(*edge, indexing='ij')
            coords = [None]*3
            coords[axes[0]], coords[axes[1]] = a, b
            coords[fixed] = np.full_like(a, self.xyz[fixed][index])
            cut = np.ix_(*[selected[k] for k in axes])
            C = np.take(self.C, index, axis=fixed)[cut]
            V = np.take(self.volume, index, axis=fixed)[cut]
            rgba = self.cmap(self.norm(C))
            rgba[..., 3] = np.where(V > 0, .38, 0.)
            ax.plot_surface(*coords, facecolors=rgba, shade=False,
                            rstride=1, cstride=1, linewidth=0, antialiased=False)

    def velocity_glyphs(self, ax):
        # Sample saved cell-centred native velocity, with fixed arrow length scaling.
        stride = max(1, int(round(1./self.dx)))
        pick = (slice(None, None, stride),)*3
        x, y, z = np.meshgrid(*self.xyz, indexing='ij')
        u = self.u[pick]; speed = np.linalg.norm(u, axis=-1)
        ch = self.ch
        visible = ((self.volume[pick] > .5*self.dx**3) & (speed > .02)
                   & (abs(x[pick]) <= ch['width_m']/2+ch['maximum_gap_m']/2+.5)
                   & (abs(y[pick]) <= ch['depth_m']/2+.5)
                   & (z[pick] >= ch['floor_m']) & (z[pick] <= ch['floor_m']+ch['height_m']))
        ax.quiver(x[pick][visible], y[pick][visible], z[pick][visible],
                  *[u[..., k][visible] for k in range(3)], length=.35,
                  normalize=False, color='#223F49', alpha=.7, linewidth=.5,
                  arrow_length_ratio=.25)

    def palm_view(self, ax, r, source=False, large=False):
        self.walls(ax,r['gap_m'],True)
        ax.add_collection3d(Poly3DCollection(self.palm,facecolors='#235238',alpha=.85,edgecolors='none'))
        if source:
            ax.add_collection3d(Poly3DCollection(self.leaves,facecolors='#76B947',edgecolors='#377939',linewidths=.2))
        else:
            self.concentration_sections(ax)
            self.velocity_glyphs(ax)
        for i,f in enumerate(self.summary['fan_configuration']):
            p=np.array(f['position']);d=np.array(f['direction'])
            ax.scatter(*p,color='#008B78',s=15)
            ax.quiver(*p,*d,length=.5,color='#008B78' if r['fans_on'] else '#777777')
            if source: ax.text(*p,f' F{i+1}',fontsize=7)
        width=self.ch['width_m']+self.ch['maximum_gap_m']+1
        depth=self.ch['depth_m']+1;height=self.ch['floor_m']+self.ch['height_m']+.3
        ax.set(xlim=(-width/2,width/2),ylim=(-depth/2,depth/2),zlim=(0,height),xlabel='x (m)',ylabel='y (m)',zlabel='z (m)')
        ax.set_box_aspect((width,depth,height),zoom=1.06 if large else 1.)
        ax.view_init(elev=25,azim=-58);ax.grid(False)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]: axis.set_pane_color((1,1,1,1))
        ax.set_title('Uptake on representative leaflets' if source else 'Three CO₂ sections · sampled 3D airflow',fontsize=10)

    def line(self, ax, r, key, label, title, colour='#008579'):
        ax.plot(self.time,[v[key] for v in self.records],color=colour,lw=1.5)
        ax.set(ylabel=label,title=title)
        ax.axvline(r['time_s'],color='#395D59',lw=.8)
        ax.set(xlim=(0,max(self.time[-1],1)),xlabel='Physical time (s)')
        ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
        ax.tick_params(labelsize=8);ax.title.set_size(10)
        ax.ticklabel_format(axis='y',style='plain',useOffset=False)

    def draw(self, view, r):
        fig=plt.figure(figsize=(14,8),facecolor='white')
        fig.text(.055,.963,'WHOLE-TREE CHAMBER  /  WATERLILY 1.8.0 + JULIA CO₂',color='#007F71',weight='bold',size=10)
        names={'c2_202502_cycle':'C2 · central fans and leaflet uptake','chamber_3d':'The chamber in three dimensions',
               'fan_sections':'CO₂ at the three fan levels',
               'operating_timeseries':'Chamber operation and CO₂ exchange','leaf_source_location':'Where the palm removes CO₂'}
        count=len(self.summary['fan_configuration'])
        if count!=3:
            names['fan_sections']=f'CO₂ at {min(count,3)} selected fan levels' if count else 'CO₂ at chamber mid-height · no fans'
        fig.text(.055,.914,names[view],weight='bold',size=23)
        status=('DIAGNOSTIC · AIRFLOW–SCALAR COUPLING SCREEN FAILED'
                if self.summary['checks']['coupling_correction_10pct_gate']=='FAIL'
                else 'DIAGNOSTIC · FAN PERFORMANCE AND GRID CONVERGENCE UNVERIFIED')
        fig.text(.945,.963,status,ha='right',size=7,color='#A46917')
        c=self.config
        fig.text(.055,.875,f"{self.ch['width_m']:g} × {self.ch['depth_m']:g} × {self.ch['height_m']:g} m · "
                 f"{c['cycle']['closed_s']/60:g} min sealed / {c['cycle']['open_phase_s']/60:g} min unsealed",size=11)
        fig.text(.945,.875,f"{r['time_s']:g} s · fans {'ON' if r['fans_on'] else 'OFF'}",ha='right',size=10)
        grid=fig.add_gridspec(2,2,left=.075,right=.88,top=.80,bottom=.16,hspace=.55,wspace=.30,
                             height_ratios=[2,1] if view!='operating_timeseries' else [1,1])
        if view=='c2_202502_cycle':
            self.palm_view(fig.add_subplot(grid[0,0],projection='3d'),r)
            self.plan(fig.add_subplot(grid[0,1]),r,self.z)
            self.line(fig.add_subplot(grid[1,0]),r,'co2_mean_ppm','Mean CO₂ (ppm)',f"Mean in fixed {self.summary['nominal_reference_volume_m3']:g} m³ region")
            self.line(fig.add_subplot(grid[1,1]),r,'co2_spatial_sd_ppm','Spatial SD (ppm)','Variation between cells; lower is more uniform')
        elif view=='chamber_3d':
            ax=fig.add_axes([.035,.17,.59,.65],projection='3d')
            self.palm_view(ax,r,large=True)
            small=fig.add_gridspec(2,1,left=.66,right=.87,top=.76,bottom=.25,hspace=.65)
            self.line(fig.add_subplot(small[0]),r,'co2_mean_ppm','Mean CO₂ (ppm)',
                      f"Fixed {self.summary['nominal_reference_volume_m3']:g} m³ chamber region")
            self.line(fig.add_subplot(small[1]),r,'co2_spatial_sd_ppm','Spatial SD (ppm)','Variation across the region')
            fig.text(.08,.16,'Coloured planes: saved CO₂ sections in x, y and z.\n'
                     'Dark arrows: native 3D velocity; 1 m/s is drawn as 0.35 m.',size=9,linespacing=1.5)
        elif view=='fan_sections':
            top=fig.add_gridspec(1,3,left=.055,right=.88,top=.8,bottom=.44,wspace=.32)
            fans=self.summary['fan_configuration']
            heights=[f['position'][2] for f in fans] or [self.ch['floor_m']+self.ch['height_m']/2]
            # Show up to three representative levels; all fan definitions stay in the metadata.
            selected=np.linspace(0,len(heights)-1,min(3,len(heights)),dtype=int)
            for i,j in enumerate(selected):
                z=int(np.argmin(abs(self.xyz[2]-heights[j])))
                self.plan(fig.add_subplot(top[0,i]),r,z)
            self.line(fig.add_subplot(grid[1,:]),r,'co2_spatial_sd_ppm','Spatial SD (ppm)','Concentration variation across the fixed chamber region')
        elif view=='operating_timeseries':
            self.line(fig.add_subplot(grid[0,0]),r,'co2_mean_ppm','Mean CO₂ (ppm)','Mean concentration')
            self.line(fig.add_subplot(grid[0,1]),r,'co2_spatial_sd_ppm','Spatial SD (ppm)','Spatial variation')
            ax=fig.add_subplot(grid[1,0]);self.line(ax,r,'gap_m','Gap (m)','Opening and fan schedule','#405DA8')
            ax.step(self.time,[v['fans_on']*self.ch['maximum_gap_m'] for v in self.records],where='post',color='#008579',label='Fans on (scaled)')
            ax.legend(frameon=False,fontsize=8)
            self.line(fig.add_subplot(grid[1,1]),r,'net_exchange_umol_s','Net CO₂ (µmol/s)','Prescribed net exchange in every phase')
        else:
            self.palm_view(fig.add_subplot(grid[0,0],projection='3d'),r,True)
            self.plan(fig.add_subplot(grid[0,1]),r,self.z,True)
            self.line(fig.add_subplot(grid[1,0]),r,'net_exchange_umol_s','Net CO₂ (µmol/s)','Whole-tree rate distributed by leaflet area')
            ax=fig.add_subplot(grid[1,1]);ax.axis('off')
            ax.text(0,.95,'240 source triangles · 120 representative leaflet patches\n\n'
                    'Trunk and bare frond axes: no CO₂ source\n\n'
                    'Triangles are clipped to cells; total uptake is conserved.\n'
                    'Leaf area is schematic, not a measured canopy area.',va='top',fontsize=10,linespacing=1.4)
        if view in ('c2_202502_cycle','chamber_3d','fan_sections'):
            cax=fig.add_axes([.925,.42,.012,.36])
            fig.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(self.norm.vmin,self.norm.vmax),cmap=self.cmap),
                         cax=cax,label='CO₂ (ppm) · fixed viridis')
        elif view=='leaf_source_location':
            maximum=self.area.max()/self.area.sum()*abs(c['exchange']['net_co2_umol_s'])
            cax=fig.add_axes([.925,.47,.012,.30])
            fig.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(0,max(maximum,1e-12)),cmap='Greens'),
                         cax=cax,label='Uptake strength (µmol/s per cell)')
        fig.text(.055,.077,f"{c['exchange']['label']} · net exchange {c['exchange']['net_co2_umol_s']:.2f} µmol/s throughout",size=10)
        fig.text(.055,.045,f"Numerical walls {self.summary['numerical_wall_thickness_m']:g} m outward · conservative first-order transport · schematic palm",size=9,color='#586764')
        return fig


def main(folder):
    # Retain the actual postprocessing code separately from the solver snapshot.
    renderer_source=Path(__file__).read_bytes()
    source_path=folder/'configs/postprocessing/render_co2.py'
    source_path.parent.mkdir(parents=True,exist_ok=True)
    if source_path.exists() and source_path.read_bytes()!=renderer_source:
        raise FileExistsError(f'Preserve the previous renderer source: {source_path}')
    source_path.write_bytes(renderer_source)
    data=Result(folder)
    plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white',
                         'font.family':'DejaVu Sans','font.size':9,'text.color':'#233334',
                         'axes.labelcolor':'#233334','xtick.color':'#233334','ytick.color':'#233334'})
    for sub in ['videos','figures']: (folder/sub).mkdir(exist_ok=True)
    fps=data.config['display']['fps'];cy=data.config['cycle'];closed=cy['closed_s']
    still_times=[0,closed,closed+cy['opening_travel_s']/2,closed+cy['opening_travel_s'],closed+cy['open_phase_s'],data.time[-1]]
    still_indices={int(np.argmin(abs(data.time-t))) for t in still_times}
    for view in VIEWS:
        mp4=folder/f'videos/{view}.mp4';gif=folder/f'videos/{view}.gif'
        if mp4.exists() or gif.exists(): raise FileExistsError(f'Preserve existing media before rerendering: {mp4}')
        writer=animation.FFMpegWriter(fps=fps,codec='libx264',extra_args=['-pix_fmt','yuv420p','-crf','20'])
        data.load(data.records[0]);placeholder=data.draw(view,data.records[0])
        with writer.saving(placeholder,str(mp4),dpi=100):
            for i,r in enumerate(data.records):
                data.load(r);fig=data.draw(view,r);writer.fig=fig
                if i in still_indices: fig.savefig(folder/f"figures/{view}_{r['time_s']:06.1f}s.png",dpi=120,facecolor='white')
                writer.grab_frame(facecolor='white');plt.close(fig)
        plt.close(placeholder)
        subprocess.run(['ffmpeg','-v','error','-y','-i',str(mp4),'-filter_complex',
                        f'fps={fps},scale=960:-1:flags=lanczos,split[a][b];[a]palettegen[p];[b][p]paletteuse',
                        '-loop','0',str(gif)],check=True)
        print(f'Rendered {view}',flush=True)
    media=sorted([*folder.glob('videos/*'),*folder.glob('figures/*.png')])
    for p in media: subprocess.run(['ffmpeg','-v','error','-xerror','-i',str(p),'-f','null','-'],check=True,stdout=subprocess.DEVNULL)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (folder/'rendering.json').write_text(json.dumps({'renderer_sha256':hashlib.sha256(renderer_source).hexdigest(),
        'renderer_source':str(source_path.relative_to(folder)),
        'media_decode':'PASS','media':[str(p.relative_to(folder)) for p in media],
        'co2_transport_enabled':True,'colour_limits_ppm':[data.norm.vmin,data.norm.vmax],'background':'white',
        'views':VIEWS,'scalar_view':'Three orthogonal saved-field sections; no volume interpolation',
        'velocity_view':'Saved native 3D cell-centred velocities; fixed glyph scale 0.35 s'},indent=2)+'\n')
    files={str(p.relative_to(folder)):sha(p) for p in sorted(folder.rglob('*')) if p.is_file() and p.name!='manifest.json'}
    (folder/'manifest.json').write_text(json.dumps({'sha256':files},indent=2)+'\n')
    print(json.dumps({'media_decode':'PASS','media_files':len(media),'output':str(folder)},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path)
    main(parser.parse_args().folder)
