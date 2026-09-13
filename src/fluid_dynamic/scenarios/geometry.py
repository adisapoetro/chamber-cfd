"""Fitted ALE geometry for sliding shells and nested open-ended sleeves.

The numerical conservation operators are the existing recreated solver. Sorted
moving/fixed planes define continuous piecewise-linear mesh trajectories. Cells
between crossing planes collapse and are born with zero inventory. All plane
crossings are time-step events. No voxel overwrite or ambient mixing source.
"""
from dataclasses import dataclass, field
import numpy as np
from ..recreated.mesh import Mesh, Study


@dataclass(frozen=True)
class ScenarioStudy(Study):
    geometry: dict = field(default_factory=dict)
    plant: dict = field(default_factory=dict)
    target_cell_m: float = 1.

    def validate(self):
        # Portable Scenario validation owns parameter types and units. Validate
        # solver scalars here as well for callers constructing this class directly.
        from dataclasses import fields
        vals={f.name:getattr(self,f.name) for f in fields(Study)}
        if not np.isfinite(list(vals.values())).all():raise ValueError('Nonfinite solver input')
        if min(self.target_cell_m,self.width,self.depth,self.height,self.gap,self.dt_s,
               self.duration_s,self.opening_s,self.closing_s,self.temperature_k,
               self.pressure_pa,self.viscosity,self.turbulent_schmidt,self.save_every_s)<=0:
            raise ValueError('Nonpositive physical scale')
        if self.wind<0 or not 0<=self.wind_angle_deg<360:raise ValueError('Invalid wind')
        if self.close_at<self.open_at+self.opening_s:raise ValueError('Overlapping motion')
        if self.period_s and self.period_s<self.close_at+self.closing_s:raise ValueError('Period cuts off closure')
        if min(self.diffusivity,self.smagorinsky,self.canopy_drag_per_m,self.trunk_drag_per_m,self.open_at)<0:
            raise ValueError('Negative coefficient')
        g=self.geometry
        if g['kind'] not in ('sliding_shells','nested_sleeves'):raise ValueError('Unsupported geometry')
        extent=self.width/2+self.gap/2 if g['kind']=='sliding_shells' else g['fixed_end_m']
        if self.outer_x/2<=extent or self.outer_y<=g.get('fixed_depth_m',self.depth) or self.outer_z<=g.get('fixed_roof_z_m',g['floor_z_m']+self.height):
            raise ValueError('Exterior domain cuts chamber travel')
        if g['kind']=='nested_sleeves':
            if not 0<g['fixed_start_m']<self.width/2 or self.width/2+self.gap/2>=g['fixed_end_m']:
                raise ValueError('Sleeves must overlap fixed enclosures and stop before end walls')


def planes(cfg,gap):
    a=cfg.width/2; mouth=gap/2
    x=[-cfg.outer_x/2,-a-mouth,-mouth,0.,mouth,a+mouth,cfg.outer_x/2]
    if cfg.geometry['kind']=='nested_sleeves':
        b=cfg.geometry['fixed_start_m'];e=cfg.geometry['fixed_end_m']
        x += [-e,-b,b,e]
    return np.sort(x)


def crossing_gaps(cfg):
    if cfg.geometry['kind']=='nested_sleeves':
        return [v for v in [2*cfg.geometry['fixed_start_m']] if 0<v<cfg.gap]
    return []


def nodes(cfg,gap):
    # Fixed number of master cells in each sorted interval, including intervals
    # that vanish at endpoints/crossings. Counts depend on maximum interval width.
    widths=np.max([np.diff(planes(cfg,g)) for g in [0.,cfg.gap,*crossing_gaps(cfg)]],axis=0)
    counts=np.maximum(1,np.ceil(widths/cfg.target_cell_m-1e-9).astype(int))
    p=planes(cfg,gap)
    return np.concatenate([np.linspace(a,b,n+1)[:-1] for a,b,n in zip(p,p[1:],counts)]+[p[-1:]]),counts


def fitted_axis(points,dx):
    p=np.unique(np.round(points,10))
    return np.concatenate([np.linspace(a,b,max(1,int(np.ceil((b-a)/dx-1e-9)))+1)[:-1] for a,b in zip(p,p[1:])]+[p[-1:]])


def panel(axis,position,lo,hi,speed=0.,hole=None):
    return dict(axis=axis,position=position,lo=lo,hi=hi,speed=speed,hole=hole)


def panels(cfg,gap):
    g=cfg.geometry;a=cfg.width/2;d=cfg.depth/2;z0=g['floor_z_m'];z1=z0+cfg.height
    out=[]
    for s in [-1,1]:
        ends=sorted([s*gap/2,s*(gap/2+a)]);lo=[ends[0],-d,z0];hi=[ends[1],d,z1]
        if g['kind']=='sliding_shells':out.append(panel(0,s*(gap/2+a),lo,hi,s*.5))
        for sign in [-1,1]:out.append(panel(1,sign*d,lo,hi,s*.5))
        for z in [z0,z1]:out.append(panel(2,z,lo,hi,s*.5))
    if g['kind']=='nested_sleeves':
        b=g['fixed_start_m'];e=g['fixed_end_m'];D=g['fixed_depth_m']/2
        Z0=g['fixed_floor_z_m'];Z1=g['fixed_roof_z_m']
        for s in [-1,1]:
            ends=sorted([s*b,s*e]);lo=[ends[0],-D,Z0];hi=[ends[1],D,Z1]
            out.append(panel(0,s*e,lo,hi))
            for sign in [-1,1]:out.append(panel(1,sign*D,lo,hi))
            for z in [Z0,Z1]:out.append(panel(2,z,lo,hi))
            # Explicit idealized stationary transition diaphragm / sliding seal.
            # The central rectangular aperture is air, not a solid box.
            out.append(panel(0,s*b,lo,hi,hole=([-1e9,-d,z0],[1e9,d,z1])))
    return out


def on_panel(pos,p):
    mask=np.isclose(pos[:,p['axis']],p['position'],atol=1e-8,rtol=0)
    for a in range(3):
        if a!=p['axis']:mask &= (pos[:,a]>p['lo'][a]-1e-9)&(pos[:,a]<p['hi'][a]+1e-9)
    if p['hole']:
        inside=np.ones(len(pos),bool)
        for a in range(3):
            if a!=p['axis']:inside&=(pos[:,a]>p['hole'][0][a])&(pos[:,a]<p['hole'][1][a])
        mask &= ~inside
    return mask


class ScenarioMesh(Mesh):
    def __init__(self,cfg,gap):
        self.cfg=cfg;self.gap=gap;g=cfg.geometry
        full,self.counts=nodes(cfg,gap);self.master_nx=sum(self.counts)
        self.columns=np.flatnonzero(np.diff(full)>1e-12)
        self.xf=np.r_[full[self.columns],full[-1]]
        yp=[-cfg.outer_y/2,-cfg.depth/2,0,cfg.depth/2,cfg.outer_y/2]
        zp=[0,g['floor_z_m'],g['floor_z_m']+cfg.height,cfg.outer_z]
        if g['kind']=='nested_sleeves':
            yp.extend([-g['fixed_depth_m']/2,g['fixed_depth_m']/2]);zp.extend([g['fixed_floor_z_m'],g['fixed_roof_z_m']])
        self.yf=fitted_axis(yp,cfg.target_cell_m);self.zf=fitted_axis(zp,cfg.target_cell_m)
        self.nodes=[self.xf,self.yf,self.zf];self.widths=[np.diff(a) for a in self.nodes]
        self.centers=[(a[1:]+a[:-1])/2 for a in self.nodes];self.shape=tuple(map(len,self.centers))
        self.master_shape=(self.master_nx,*self.shape[1:])
        self.xyz=np.stack(np.meshgrid(*self.centers,indexing='ij'),axis=-1)
        self.volume=(self.widths[0][:,None,None]*self.widths[1][None,:,None]*self.widths[2][None,None,:]).ravel()
        self.size=len(self.volume);ids=np.arange(self.size).reshape(self.shape);self.faces=[]
        self.panels=panels(cfg,gap)
        for axis in range(3):
            sl=[slice(None)]*3;sr=sl.copy();sl[axis]=slice(None,-1);sr[axis]=slice(1,None)
            left=ids[tuple(sl)].ravel();right=ids[tuple(sr)].ravel()
            pos=(self.xyz[tuple(sl)]+self.xyz[tuple(sr)])/2;shp=[1]*3;shp[axis]=-1
            pos[...,axis]=np.broadcast_to(self.nodes[axis][1:-1].reshape(shp),pos.shape[:-1]);pos=pos.reshape(-1,3)
            distance=(self.xyz[tuple(sr)][...,axis]-self.xyz[tuple(sl)][...,axis]).ravel()
            area=(self.volume.reshape(self.shape)/self.widths[axis].reshape(shp))[tuple(sl)].ravel()
            wall=np.zeros(len(left),bool)
            for p in self.panels:
                if p['axis']==axis:wall|=on_panel(pos,p)
            self.faces.append(dict(left=left,right=right,axis=axis,area=area,distance=distance,wall=wall,position=pos,normal=np.ones(len(left))))
        self.boundaries=[]
        for axis in range(3):
            for side in [-1,1]:
                sl=[slice(None)]*3;sl[axis]=0 if side<0 else -1;cells=ids[tuple(sl)].ravel()
                pos=self.xyz[tuple(sl)].reshape(-1,3).copy();pos[:,axis]=self.nodes[axis][0 if side<0 else -1]
                w=self.widths[axis][0 if side<0 else -1]
                kind='wall' if axis==2 and side==-1 else ('inlet' if cfg.wind_vector[axis]*side<0 else 'open')
                self.boundaries.append(dict(cell=cells,axis=axis,side=side,area=self.volume[cells]/w,distance=np.full(len(cells),w/2),kind=kind,position=pos))

    def volumes_at(self,gap):
        full,_=nodes(self.cfg,gap)
        return (np.diff(full)[self.columns,None,None]*self.widths[1][None,:,None]*self.widths[2][None,None,:]).ravel()

    def sweep(self,g0,g1,dt):
        old,_=nodes(self.cfg,g0);new,_=nodes(self.cfg,g1);speed=(new-old)/dt
        local=np.r_[speed[self.columns],speed[-1]]
        qi=[np.broadcast_to(local[1:-1,None,None],(self.shape[0]-1,*self.shape[1:])).ravel()*f['area'] if f['axis']==0 else np.zeros(len(f['left'])) for f in self.faces]
        return qi,[np.zeros(len(f['cell'])) for f in self.boundaries]

    def wall_velocity(self,positions,gap_rate):
        result=np.zeros_like(positions)
        for p in self.panels:
            if p['speed']:result[on_panel(positions,p),0]=p['speed']*gap_rate
        return result

    def boundary_wall_velocity(self,f,rate):
        if f['kind']=='inlet':return np.broadcast_to(self.cfg.wind_vector,f['position'].shape).copy()
        return self.wall_velocity(f['position'],rate) if f['kind']=='wall' else np.zeros_like(f['position'])

    def box_weights(self,lo,hi):
        w=[np.maximum(0,np.minimum(a[1:],hi[i])-np.maximum(a[:-1],lo[i])) for i,a in enumerate(self.nodes)]
        return (w[0][:,None,None]*w[1][None,:,None]*w[2][None,None,:]).ravel()

    def reporting_weights(self,kind):
        cfg=self.cfg;g=cfg.geometry;z=g['floor_z_m'];d=cfg.depth/2;a=cfg.width/2
        if g['kind']=='sliding_shells':
            if kind=='fixed':return self.box_weights([-a,-d,z],[a,d,z+cfg.height])
            return sum(self.box_weights([lo,-d,z],[hi,d,z+cfg.height]) for lo,hi in [(-a-self.gap/2,-self.gap/2),(self.gap/2,a+self.gap/2)])
        b=g['fixed_start_m'];e=g['fixed_end_m'];D=g['fixed_depth_m']/2
        fixed=sum(self.box_weights([lo,-D,g['fixed_floor_z_m']],[hi,D,g['fixed_roof_z_m']]) for lo,hi in [(-e,-b),(b,e)])
        if kind=='fixed':return fixed+self.box_weights([-b,-d,z],[b,d,z+cfg.height])
        # Union of stationary side enclosures and sleeve-covered central region.
        if self.gap/2<b:
            fixed+=sum(self.box_weights([lo,-d,z],[hi,d,z+cfg.height]) for lo,hi in [(-b,-self.gap/2),(self.gap/2,b)])
        return fixed

    def probes(self):
        z=self.cfg.geometry['floor_z_m'];h=self.cfg.height
        return [('centre_low',(0,0,z+.2*h)),('centre_mid',(0,0,z+.5*h)),('centre_high',(0,0,z+.8*h)),
                ('left_mid',(-.35*self.cfg.width,0,z+.5*h)),('right_mid',(.35*self.cfg.width,0,z+.5*h))]

    def mouth_selector(self,positions):
        x,y,z=positions.T;floor=self.cfg.geometry['floor_z_m']
        return np.isclose(abs(x),self.gap/2,rtol=0,atol=1e-8)&(abs(y)<self.cfg.depth/2)&(z>floor)&(z<floor+self.cfg.height)

    def canopy(self):
        return self.canopy_at(self.gap)

    def canopy_at(self,gap):
        p=self.cfg.plant
        if not p or p['height_m']==0:return np.zeros(self.size),np.zeros(self.size)
        full,_=nodes(self.cfg,gap)
        xc=((full[:-1]+full[1:])/2)[self.columns]
        x,y,z=(a.ravel() for a in np.meshgrid(xc,self.centers[1],self.centers[2],indexing='ij'))
        z=z-self.cfg.geometry['floor_z_m'];h=p['height_m'];rx,ry=p['crown_radii_m']
        crown=np.exp(-((x/(.8*rx))**4+(y/(.8*ry))**4+((z-.62*h)/(.32*h))**4))
        crown*=(abs(x)<rx)&(abs(y)<ry)&(z>.15*h)&(z<h)
        trunk=np.exp(-((x/p['trunk_radius_m'])**2+(y/p['trunk_radius_m'])**2))*((z>0)&(z<.45*h))
        return crown,trunk

    def source_weights_at(self,gap):
        """Endpoint inventory weights on the current master-column mapping.

        A vanishing/newborn cell has exactly zero source inventory at that
        endpoint. Each half source step still integrates the prescribed total.
        """
        return self.canopy_at(gap)[0]*self.volumes_at(gap)

    def motion_events(self):
        c=self.cfg;out=[]
        for shift in (np.arange(0,c.duration_s+c.period_s,c.period_s) if c.period_s else [0.]):
            for g in crossing_gaps(c):
                out.extend([shift+c.open_at+c.opening_s*g/c.gap,shift+c.close_at+c.closing_s*(1-g/c.gap)])
        return out
