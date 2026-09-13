"""Conservative, regularized internal fan momentum forcing.

The nominal free-jet momentum scale rho Q_free**2/A is an assumption, not an
installed pressure-flow curve. Force is integrated over finite-volume cells;
fans create neither air volume nor CO2. No fitted scalar mixing time is used.
"""
from dataclasses import asdict, dataclass
import numpy as np
from scipy.special import erf
from .mesh import opening


@dataclass(frozen=True)
class FanArray:
    model: str = 'Delta AFC1212DE'
    heights_above_floor_m: tuple = (3.7/6, 3.7/2, 5*3.7/6)
    closed_x_m: float = -1.
    closed_y_m: float = -2.75
    axis: tuple = (0., 1., 0.)
    mounting_shell: int = -1
    diameter_m: float = .12
    nominal_free_air_cfm: float = 148.3
    momentum_factor: float = 1.
    kernel_sigma_m: float = .3
    enabled: bool = True
    schedule: str = 'closed_only'
    # Optional per-height discharge vectors. None preserves the original array.
    axes: tuple | None = None

    def validate(self, cfg):
        values=[*self.heights_above_floor_m,self.closed_x_m,self.closed_y_m,
                *self.axis,self.diameter_m,self.nominal_free_air_cfm,
                self.momentum_factor,self.kernel_sigma_m]
        if not np.isfinite(values).all():
            raise ValueError('Nonfinite fan parameter')
        if (not self.heights_above_floor_m and self.enabled) or len(self.axis)!=3:
            raise ValueError('Fan heights and three-component axis are required')
        if not np.isclose(np.linalg.norm(self.axis),1.,rtol=0,atol=1e-12):
            raise ValueError('Fan axis must be a unit vector')
        if self.axes is not None:
            axes=np.asarray(self.axes,dtype=float)
            if axes.shape!=(len(self.heights_above_floor_m),3) or not np.isfinite(axes).all():
                raise ValueError('One finite three-component axis is required per fan height')
            if not np.allclose(np.linalg.norm(axes,axis=1),1.,rtol=0,atol=1e-12):
                raise ValueError('Every fan axis must be a unit vector')
        if self.schedule!='closed_only' or cfg.geometry['kind']!='sliding_shells':
            raise ValueError('This extension supports closed-only sliding-shell fans')
        # Zero denotes a stationary support, independent of either moving shell.
        if self.mounting_shell not in (-1,0,1) or (self.mounting_shell!=0 and self.mounting_shell*self.closed_x_m<=0):
            raise ValueError('Fan position and mounting shell disagree')
        if min(self.diameter_m,self.nominal_free_air_cfm,self.kernel_sigma_m)<=0 or self.momentum_factor<0:
            raise ValueError('Invalid fan force or length scale')
        if self.kernel_sigma_m<self.diameter_m/2:
            raise ValueError('Kernel narrower than fan radius is unsupported')
        if not isinstance(self.enabled,bool):
            raise ValueError('Fan enabled must be boolean')
        if not self.heights_above_floor_m:
            return  # Explicit disabled, zero-fan configuration.
        r=self.diameter_m/2
        if (abs(self.closed_x_m)+r>=cfg.width/2 or abs(self.closed_y_m)+r>=cfg.depth/2
            or min(self.heights_above_floor_m)<=r or max(self.heights_above_floor_m)+r>=cfg.height):
            raise ValueError('A fan lies outside the chamber airspace')

    def active(self, cfg, time_s):
        # Test the right-hand phase at an event: off at opening start, on when
        # closure is complete. Solver steps are already split at every event.
        t=time_s % cfg.period_s if cfg.period_s else time_s
        closed=(t<cfg.open_at or t>=cfg.close_at+cfg.closing_s)
        return bool(self.enabled and self.heights_above_floor_m and self.momentum_factor>0 and closed
                    and opening(cfg,time_s)==0)

    def positions(self, cfg, gap):
        z=np.asarray(self.heights_above_floor_m)+cfg.geometry['floor_z_m']
        return np.column_stack((np.full(len(z),self.closed_x_m+self.mounting_shell*gap/2),
                                np.full(len(z),self.closed_y_m),z))

    def directions(self):
        return (np.broadcast_to(np.asarray(self.axis),(len(self.heights_above_floor_m),3))
                if self.axes is None else np.asarray(self.axes,dtype=float))

    def properties(self, cfg):
        q=self.nominal_free_air_cfm*.0004719474432
        area=np.pi*self.diameter_m**2/4
        rho=cfg.pressure_pa/(287.05*cfg.temperature_k)
        inputs=asdict(self)
        if self.axes is None:
            inputs.pop('axes')  # Keep historical serialized properties unchanged.
        result=dict(inputs,nominal_free_air_m3_s_per_fan=q,disc_area_m2=area,
                    nominal_disc_velocity_m_s=q/area,air_density_kg_m3=rho,
                    force_n_per_fan=rho*q*q/area*self.momentum_factor,
                    force_interpretation='Prescribed free-jet momentum scale; installed thrust/flow unverified',
                    placement_status='Owner three height bands; horizontal placement and direction assumed',
                    world_positions_closed_m=self.positions(cfg,0).tolist())
        if self.mounting_shell==0:
            result['placement_status']='Fixed pole; prescribed heights and discharge vectors; mounting and installed flow require field verification'
            result['support_geometry']='Schematic support only; pole diameter and blockage not resolved'
        if self.axes is not None:
            result['world_discharge_unit_vectors']=self.directions().tolist()
            result['azimuth_deg_ccw_from_positive_x']=(np.degrees(np.arctan2(self.directions()[:,1],self.directions()[:,0]))%360).tolist()
        return result


class FanForcing:
    def __init__(self, cfg, fans):
        fans.validate(cfg)
        self.cfg=cfg;self.fans=fans;self._key=None;self._source=None

    def source(self, mesh, time_s):
        """Return cell-integrated acceleration [m4/s2] for momentum inventory."""
        f=self.fans;c=self.cfg
        if not f.active(c,time_s):
            return np.zeros((mesh.size,3))
        if mesh.gap!=0:
            raise ValueError('Active fans require a sealed mesh')
        key=(mesh.shape,mesh.gap)
        if self._key!=key:
            source=np.zeros((mesh.size,3))
            floor=c.geometry['floor_z_m']
            bounds=((-c.width/2,c.width/2),(-c.depth/2,c.depth/2),(floor,floor+c.height))
            props=f.properties(c)
            force_over_rho=props['force_n_per_fan']/props['air_density_kg_m3']
            for centre,axis in zip(f.positions(c,mesh.gap),f.directions()):
                integrals=[]
                for edges,mu,(lo,hi) in zip(mesh.nodes,centre,bounds):
                    lower=np.maximum(edges[:-1],lo);upper=np.minimum(edges[1:],hi)
                    w=.5*(erf((upper-mu)/(np.sqrt(2)*f.kernel_sigma_m))
                          -erf((lower-mu)/(np.sqrt(2)*f.kernel_sigma_m)))
                    integrals.append(np.where(upper>lower,w,0.))
                weights=(integrals[0][:,None,None]*integrals[1][None,:,None]*integrals[2][None,None,:]).ravel()
                if weights.sum()<=0:
                    raise ValueError('No fluid support for fan force')
                source+=weights[:,None]/weights.sum()*force_over_rho*axis
            self._key=key;self._source=source
        return self._source
