"""Conservative leaflet-area allocation on the existing moving ALE grid.

Sutherland–Hodgman clips each triangle to fixed y/z cell strips. A triangular
area CDF integrates the remaining x cuts exactly (to floating-point precision).
No centroid binning, trunk sink, fitted diffusion, or extra CO2 source is used.
"""
from collections import OrderedDict
import numpy as np
from .leaf_geometry import triangle_area
from ..scenarios.geometry import ScenarioMesh, nodes


def clip_polygon(polygon, axis, bound, keep_above):
    p = np.asarray(polygon, dtype=float)
    if len(p) == 0:
        return p.reshape(0, 3)
    out = []
    previous = p[-1]
    previous_in = previous[axis] >= bound if keep_above else previous[axis] <= bound
    for current in p:
        current_in = current[axis] >= bound if keep_above else current[axis] <= bound
        if current_in != previous_in:
            fraction = (bound - previous[axis]) / (current[axis] - previous[axis])
            point = previous + fraction * (current - previous)
            point[axis] = bound
            out.append(point)
        if current_in:
            out.append(current)
        previous, previous_in = current, current_in
    return np.asarray(out, dtype=float).reshape(-1, 3)


def intersected_bins(lo, hi, faces):
    # A polygon exactly in a grid plane belongs once, to the upper cell.
    if lo == hi:
        i = min(len(faces)-2, int(np.searchsorted(faces, lo, side='right')-1))
        return range(max(0, i), max(0, i)+1)
    start = max(0, int(np.searchsorted(faces, lo, side='right')-1))
    stop = min(len(faces)-1, int(np.searchsorted(faces, hi, side='left')))
    return range(start, stop)


class LeafSurfaceGrid:
    def __init__(self, triangles, yf, zf):
        triangles = np.asarray(triangles, dtype=float)
        self.yf, self.zf = np.array(yf), np.array(zf)
        if triangles.ndim != 3 or triangles.shape[1:] != (3, 3) or not np.isfinite(triangles).all():
            raise ValueError('Expected finite triangles (N, 3, 3)')
        if any(not np.all(np.diff(f) > 0) for f in (self.yf, self.zf)):
            raise ValueError('Fixed grid faces must be increasing')
        if any(triangles[..., axis].min() < f[0] or triangles[..., axis].max() > f[-1]
               for axis, f in ((1, self.yf), (2, self.zf))):
            raise ValueError('Leaf geometry leaves the CFD domain')
        fragments, bins = [], []
        for tri in triangles:
            for j in intersected_bins(tri[:, 1].min(), tri[:, 1].max(), self.yf):
                p = clip_polygon(clip_polygon(tri, 1, self.yf[j], True), 1, self.yf[j+1], False)
                if len(p) < 3:
                    continue
                for k in intersected_bins(p[:, 2].min(), p[:, 2].max(), self.zf):
                    q = clip_polygon(clip_polygon(p, 2, self.zf[k], True), 2, self.zf[k+1], False)
                    for v in range(1, len(q)-1):
                        fragments.append([q[0], q[v], q[v+1]])
                        bins.append(j*(len(self.zf)-1)+k)
        f = np.asarray(fragments, dtype=float)
        if not len(f):
            raise ValueError('No resolved leaflet area')
        areas = triangle_area(f)
        positive = areas > 0
        self.area = areas[positive]
        self.x = np.sort(f[positive, :, 0], axis=1)
        self.bins = np.array(bins)[positive]
        self.total_area = float(triangle_area(triangles).sum())
        if not np.isclose(self.area.sum(), self.total_area, rtol=2e-12, atol=1e-13):
            raise ValueError('Clipping did not preserve surface area')
        self.cache = OrderedDict()

    def areas(self, xf):
        xf = np.asarray(xf, dtype=float)
        if not np.isfinite(xf).all() or np.any(np.diff(xf) < 0):
            raise ValueError('Moving faces must be finite and nondecreasing')
        key = xf.tobytes()
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key]
        a, b, c = self.x.T
        if xf[0] > a.min() or xf[-1] < c.max():
            raise ValueError('Leaf geometry leaves the x domain')
        cumulative = []
        nbin = (len(self.yf)-1)*(len(self.zf)-1)
        for cut in xf:
            fraction = np.zeros(len(a))
            fraction[cut >= c] = 1
            low = (cut > a) & (cut < b)
            high = (cut >= b) & (cut < c)
            fraction[low] = (cut-a[low])**2 / ((b[low]-a[low])*(c[low]-a[low]))
            fraction[high] = 1-(c[high]-cut)**2 / ((c[high]-a[high])*(c[high]-b[high]))
            # For x-constant triangles use a right-cell, half-open convention.
            fraction[(a == c) & (cut == a)] = 0
            cumulative.append(np.bincount(self.bins, weights=self.area*fraction, minlength=nbin))
        result = np.diff(cumulative, axis=0)
        if result.min() < -1e-12:
            raise ValueError('Negative intersected area')
        result = np.maximum(result, 0).reshape(len(xf)-1, len(self.yf)-1, len(self.zf)-1)
        if not np.isclose(result.sum(), self.total_area, rtol=3e-12, atol=1e-12):
            raise ValueError('Moving-grid source lost leaflet area')
        result.setflags(write=False)
        self.cache[key] = result
        if len(self.cache) > 256:
            self.cache.popitem(last=False)
        return result


class LeafMesh(ScenarioMesh):
    def __init__(self, cfg, gap, *, leaf_grid):
        super().__init__(cfg, gap)
        self.leaf_grid = leaf_grid
        if not np.array_equal(self.yf, leaf_grid.yf) or not np.array_equal(self.zf, leaf_grid.zf):
            raise ValueError('Leaf map and CFD grid differ')

    def source_weights_at(self, gap):
        full, _ = nodes(self.cfg, gap)
        return self.leaf_grid.areas(full)[self.columns].ravel()


def mesh_factory(cfg, geometry):
    """One immutable surface clip per run, reused at every moving endpoint."""
    with np.load(geometry) as g:
        triangles = g['triangles'][g['source_mask']]
    mesh = ScenarioMesh(cfg, 0.)
    grid = LeafSurfaceGrid(triangles, mesh.yf, mesh.zf)
    def factory(config, gap):
        return LeafMesh(config, gap, leaf_grid=grid)
    return factory
