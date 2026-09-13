"""Rectilinear ALE mesh fitted to the two zero-thickness translating shells.

The gap block collapses exactly at closure. Its cells have zero inventory before
birth/after collapse; the implicit transport equation supplies/removes inventory
through space-time face fluxes, never through an ambient fill operation.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Study:
    resolution: int = 2  # cells per metre in the chamber and full gap
    outer_x: float = 14.
    outer_y: float = 10.
    outer_z: float = 9.
    width: float = 4.
    depth: float = 4.
    height: float = 6.
    gap: float = 2.
    wind: float = 1.5
    wind_angle_deg: float = 0.  # flow TO, counterclockwise from +x toward +y
    inside_ppm: float = 340.
    ambient_ppm: float = 400.
    temperature_k: float = 298.15
    pressure_pa: float = 101325.
    viscosity: float = 1.5e-5
    diffusivity: float = 1.6e-5
    smagorinsky: float = .12
    turbulent_schmidt: float = .7
    canopy_drag_per_m: float = .3
    trunk_drag_per_m: float = 20.
    tree_source_umol_s: float = 0.
    open_at: float = 10.
    opening_s: float = 20.
    close_at: float = 1.e9
    closing_s: float = 20.
    period_s: float = 0.
    duration_s: float = 240.
    dt_s: float = .2
    save_every_s: float = 2.
    flux_correction: bool = True
    local_scalar_bounds: bool = True
    recovery_dwell_s: float = 30.
    start_open: bool = False

    def validate(self):
        from dataclasses import asdict
        if not np.isfinite(list(asdict(self).values())).all():
            raise ValueError("All case inputs must be finite")
        if not isinstance(self.resolution,int) or self.resolution<1:
            raise ValueError("Resolution must be a positive integer")
        if min(self.width,self.depth,self.height,self.gap,self.temperature_k,
               self.pressure_pa,self.viscosity,self.turbulent_schmidt,
               self.opening_s,self.closing_s,self.duration_s,self.dt_s,self.save_every_s)<=0:
            raise ValueError("Invalid physical or time scale")
        if min(self.diffusivity,self.smagorinsky,self.canopy_drag_per_m,self.trunk_drag_per_m,
               self.recovery_dwell_s,self.open_at,self.period_s)<0:
            raise ValueError("Negative coefficient or clock")
        if self.wind<0 or not 0<=self.wind_angle_deg<360:
            raise ValueError("Wind speed must be nonnegative and flow-to angle in [0, 360)")
        if self.close_at<self.open_at+self.opening_s:
            raise ValueError("Travel windows must not overlap")
        if self.period_s and self.period_s<self.close_at+self.closing_s:
            raise ValueError("Cycle period cuts off closure")
        if self.outer_x<=self.width+self.gap or self.outer_y<=self.depth or self.outer_z<=self.height:
            raise ValueError("Ambient domain does not contain chamber travel")
        for v in (self.width/2,self.depth,self.height,self.gap,self.outer_z,
                  (self.outer_x-self.width)/2,(self.outer_y-self.depth)/2):
            if not np.isclose(v*self.resolution,round(v*self.resolution)):
                raise ValueError("Block lengths must align with the selected grid family")

    @property
    def wind_vector(self):
        angle=np.deg2rad(self.wind_angle_deg)
        direction=np.array([np.cos(angle),np.sin(angle),0.])
        direction[np.abs(direction)<1e-14]=0.  # exact cardinal boundary classification
        return self.wind*direction


def opening(cfg,t):
    if cfg.period_s:t=t%cfg.period_s
    if cfg.start_open and t<=cfg.close_at:return cfg.gap
    if t<=cfg.open_at:return 0.
    if t<cfg.open_at+cfg.opening_s:return cfg.gap*(t-cfg.open_at)/cfg.opening_s
    if t<=cfg.close_at:return cfg.gap
    if t<cfg.close_at+cfg.closing_s:return cfg.gap*(1-(t-cfg.close_at)/cfg.closing_s)
    return 0.


def master_nodes(cfg,gap):
    r=cfg.resolution; half=cfg.width/2; outer=cfg.outer_x/2
    limits=[-outer,-half-gap/2,-gap/2,gap/2,half+gap/2,outer]
    counts=[round((outer-half)*r),round(half*r),round(cfg.gap*r),round(half*r),round((outer-half)*r)]
    nodes=np.concatenate([np.linspace(a,b,n+1)[:-1] for a,b,n in zip(limits,limits[1:],counts)]+[np.array([outer])])
    return nodes,counts


class Mesh:
    def __init__(self,cfg,gap):
        self.cfg=cfg;self.gap=gap
        full,counts=master_nodes(cfg,gap);self.counts=counts
        self.master_nx=sum(counts)
        self.columns=np.flatnonzero(np.diff(full)>1e-13)
        self.xf=np.r_[full[self.columns],full[-1]]
        self.yf=np.linspace(-cfg.outer_y/2,cfg.outer_y/2,round(cfg.outer_y*cfg.resolution)+1)
        self.zf=np.linspace(0,cfg.outer_z,round(cfg.outer_z*cfg.resolution)+1)
        self.nodes=[self.xf,self.yf,self.zf]
        self.widths=[np.diff(a) for a in self.nodes]
        self.centers=[.5*(a[1:]+a[:-1]) for a in self.nodes]
        self.shape=tuple(len(a) for a in self.centers)
        self.master_shape=(self.master_nx,*self.shape[1:])
        self.xyz=np.stack(np.meshgrid(*self.centers,indexing='ij'),axis=-1)
        self.volume=(self.widths[0][:,None,None]*self.widths[1][None,:,None]*self.widths[2][None,None,:]).ravel()
        self.size=len(self.volume)
        ids=np.arange(self.size).reshape(self.shape)
        self.faces=[]
        half=cfg.width/2
        for axis in range(3):
            lower=[slice(None)]*3;upper=lower.copy();lower[axis]=slice(None,-1);upper[axis]=slice(1,None)
            left=ids[tuple(lower)].ravel();right=ids[tuple(upper)].ravel()
            pos=.5*(self.xyz[tuple(lower)]+self.xyz[tuple(upper)])
            # Actual face position, not midpoint of unequal cell centres.
            shp=[1]*3;shp[axis]=-1
            pos[...,axis]=np.broadcast_to(self.nodes[axis][1:-1].reshape(shp),pos.shape[:-1])
            pos=pos.reshape(-1,3)
            dist=(self.xyz[tuple(upper)][...,axis]-self.xyz[tuple(lower)][...,axis]).ravel()
            area=(self.volume.reshape(self.shape)/self.widths[axis].reshape(shp))[tuple(lower)].ravel()
            x,y,z=pos.T
            in_half=(np.abs(x)>gap/2-1e-10)&(np.abs(x)<half+gap/2+1e-10)
            if axis==0:
                wall=np.isclose(abs(x),half+gap/2,atol=1e-10)&(abs(y)<cfg.depth/2)&(z<cfg.height)
            elif axis==1:
                wall=np.isclose(abs(y),cfg.depth/2,atol=1e-10)&in_half&(z<cfg.height)
            else:
                wall=np.isclose(z,cfg.height,atol=1e-10)&in_half&(abs(y)<cfg.depth/2)
            self.faces.append(dict(left=left,right=right,axis=axis,area=area,distance=dist,
                                   wall=wall,position=pos,normal=np.ones(len(left))))
        self.boundaries=[]
        for axis in range(3):
            for side in (-1,1):
                sl=[slice(None)]*3;sl[axis]=0 if side<0 else -1
                cells=ids[tuple(sl)].ravel()
                position=self.xyz[tuple(sl)].reshape(-1,3).copy()
                position[:,axis]=self.nodes[axis][0 if side<0 else -1]
                cell_width=(self.widths[axis][0 if side<0 else -1])
                area=self.volume[cells]/cell_width
                kind='wall' if axis==2 and side==-1 else ('inlet' if cfg.wind_vector[axis]*side<0 else 'open')
                self.boundaries.append(dict(cell=cells,axis=axis,side=side,area=area,
                    distance=np.full(len(cells),cell_width/2),kind=kind,position=position))

    def volumes_at(self,gap):
        xf,_=master_nodes(self.cfg,gap)
        return (np.diff(xf)[self.columns,None,None]*self.widths[1][None,:,None]*self.widths[2][None,None,:]).ravel()

    def sweep(self,gap0,gap1,dt):
        old,_=master_nodes(self.cfg,gap0);new,_=master_nodes(self.cfg,gap1)
        speed=(new-old)/dt
        local=np.r_[speed[self.columns],speed[-1]]
        internal=[]
        for f in self.faces:
            if f['axis']==0:
                w=np.broadcast_to(local[1:-1,None,None],(self.shape[0]-1,*self.shape[1:])).ravel()*f['area']
            else:w=np.zeros(len(f['left']))
            internal.append(w)
        boundary=[]
        for f in self.boundaries:
            w=local[0 if f['side']<0 else -1]*f['side']*f['area'] if f['axis']==0 else np.zeros(len(f['cell']))
            boundary.append(w)
        return internal,boundary

    def wall_velocity(self,positions,gap_rate):
        velocity=np.zeros_like(positions)
        velocity[:,0]=np.sign(positions[:,0])*gap_rate/2
        return velocity

    def boundary_wall_velocity(self,f,gap_rate):
        result=np.zeros_like(f['position'])
        if f['kind']=='inlet':result[:]=self.cfg.wind_vector
        if f['kind']=='wall':
            x,y,_=f['position'].T
            inside=(abs(x)>self.gap/2)&(abs(x)<self.cfg.width/2+self.gap/2)&(abs(y)<self.cfg.depth/2)
            result[inside,0]=np.sign(x[inside])*gap_rate/2
        return result

    def collect(self,master):
        return master[self.columns].reshape(self.size,*master.shape[3:])

    def store(self,values):
        a=np.zeros((*self.master_shape,*values.shape[1:]))
        a[self.columns]=values.reshape(*self.shape,*values.shape[1:])
        return a

    def net_outflow(self,internal,boundary):
        result=np.zeros(self.size)
        for f,q in zip(self.faces,internal):
            np.add.at(result,f['left'],q);np.add.at(result,f['right'],-q)
        for f,q in zip(self.boundaries,boundary):np.add.at(result,f['cell'],q)
        return result

    def canopy(self):
        x,y,z=self.xyz.reshape(-1,3).T
        # Explicit synthetic porous crown + wood drag; no solid gas-storage cells.
        crown=np.exp(-((x/1.3)**4+(y/1.3)**4+((z-2.2)/1.1)**4))
        crown*=((abs(x)<1.6)&(abs(y)<1.6)&(z>.5)&(z<3.5))
        trunk=np.exp(-((x/.15)**2+(y/.15)**2))*((z>0)&(z<1.8))
        return crown,trunk
