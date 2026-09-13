"""Trait-conditioned XPalm.VPalm export; no plant-model code is modified."""
import json
from pathlib import Path
import numpy as np


def triangle_area(triangles):
    t = np.asarray(triangles, dtype=float)
    return np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2


def condition_export(raw, destination, *, height_m, radius_m, trunk_radius_m, floor_m):
    """Apply and report an envelope transformation, not an allometry fit.

    Measured frond lengths prescribe the raw VPalm input. This last display/source
    transformation changes those lengths; both the factors and final axis lengths
    are exported so the fitted envelope is not mistaken for exact organ geometry.
    """
    raw, destination = Path(raw), Path(destination)
    a = np.genfromtxt(raw / 'triangles.csv', delimiter=',', names=True, dtype=None, encoding='utf8')
    b = np.genfromtxt(raw / 'axes.csv', delimiter=',', names=True, dtype=None, encoding='utf8')
    triangles = np.stack([a[k] for k in a.dtype.names[3:]], axis=-1).reshape(-1, 3, 3)
    axes = np.stack([b[k] for k in b.dtype.names[3:]], axis=-1).reshape(-1, 2, 3)
    original_axes = axes.copy()
    leaf = a['tissue'] == 'Leaflet'
    source = leaf & (a['frond_rank'] >= 1)
    base = min(0., float(triangles[..., 2].min()))
    raw_radius = float(np.linalg.norm(triangles[leaf, :, :2], axis=-1).max())
    raw_height = float(triangles[..., 2].max() - base)
    scale = np.array([radius_m / raw_radius, radius_m / raw_radius, height_m / raw_height])
    triangles[..., 2] -= base
    axes[..., 2] -= base
    triangles *= scale
    axes *= scale
    trunk = a['tissue'] == 'Internode'
    stem_radius = float(np.linalg.norm(triangles[trunk, :, :2], axis=-1).max())
    triangles[trunk, :, :2] *= trunk_radius_m / stem_radius
    triangles[..., 2] += floor_m
    axes[..., 2] += floor_m
    area = triangle_area(triangles[source])
    if not np.isfinite(triangles).all() or area.sum() <= 0:
        raise ValueError('Invalid reconstructed geometry')
    destination.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(destination / 'palm.npz', triangles=triangles,
        tissue=a['tissue'], frond_rank=a['frond_rank'], node_id=a['node_id'],
        source_mask=source, axes=axes, axis_tissue=b['tissue'], axis_rank=b['frond_rank'])
    # Portable, labelled triangle OBJ: groups preserve tissue and frond rank.
    with (destination / 'palm.obj').open('w') as f:
        f.write('# XPalm.VPalm conditional C2 reconstruction; coordinates m, chamber floor included\n')
        for point in triangles.reshape(-1, 3):
            f.write('v ' + ' '.join(f'{v:.10g}' for v in point) + '\n')
        previous = None
        for i, (tissue, rank) in enumerate(zip(a['tissue'], a['frond_rank'])):
            group = f'{tissue}_rank_{rank}'
            if group != previous:
                f.write(f'g {group}\n'); previous = group
            f.write(f'f {3*i+1} {3*i+2} {3*i+3}\n')
    organ_lengths = []
    for rank in range(1, 34):
        row = dict(frond_rank=rank)
        for tissue in ('PetioleSegment', 'RachisSegment'):
            mask = (b['frond_rank'] == rank) & (b['tissue'] == tissue)
            for label, coordinates in [('raw', original_axes), ('conditioned', axes)]:
                row[f'{label}_{tissue}_length_m'] = float(np.linalg.norm(
                    coordinates[mask, 1] - coordinates[mask, 0], axis=1).sum())
        organ_lengths.append(row)
    meta = dict(status='CONDITIONAL_STATIC_RECONSTRUCTION_NOT_VALIDATED',
        coordinate_unit='m', floor_z_m=floor_m, height_above_floor_m=height_m,
        maximum_crown_radius_proxy_m=radius_m, raw_height_m=raw_height,
        raw_maximum_radius_m=raw_radius, xyz_scale_factors=scale.tolist(),
        stem_radius_proxy_m=trunk_radius_m, stem_xy_extra_factor=trunk_radius_m/stem_radius,
        triangles=len(triangles), expanded_leaflet_triangles=int(source.sum()),
        expanded_leaflet_nodes=len(np.unique(a['node_id'][source])),
        expanded_frond_ranks=np.unique(a['frond_rank'][source]).tolist(),
        excluded_spear_ranks=np.unique(a['frond_rank'][leaf & ~source]).tolist(),
        reconstructed_one_sided_leaf_area_m2=float(area.sum()),
        per_frond_leaf_area_m2={str(r): float(area[a['frond_rank'][source] == r].sum())
            for r in np.unique(a['frond_rank'][source])},
        axis_lengths=organ_lengths,
        bounds_m=[triangles.min((0, 1)).tolist(), triangles.max((0, 1)).tolist()],
        limitation='Envelope conditioning changes organ lengths. Fine leaflets, standing-frond interpretation, azimuth, inclination and stem height are reconstruction assumptions. Net exchange is allocated uniformly per reconstructed expanded leaflet area; no measured per-leaf physiology.')
    (destination / 'geometry.json').write_text(json.dumps(meta, indent=2, allow_nan=False)+'\n')
    return meta
