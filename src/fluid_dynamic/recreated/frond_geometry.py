"""Representative leaflets on the preferred illustration's five fronds.

The display mesh is unchanged. These small planar source patches represent
leaf-bearing frond regions; their count and area are not plant measurements.
"""
import json
from pathlib import Path
import numpy as np
from ..scenarios.presentation_render import fitted_palm
from .leaf_geometry import triangle_area


def build_geometry(root, cfg, parameters):
    display = np.asarray(fitted_palm(root, cfg))
    if display.shape != (224, 3, 3):
        raise ValueError('Preferred illustration topology changed')
    if not 0 < parameters['bare_base_fraction'] < parameters['last_leaflet_fraction'] < 1:
        raise ValueError('Invalid leaf-bearing interval')
    n = parameters['pairs_per_frond']
    if not isinstance(n, int) or n < 2:
        raise ValueError('At least two leaflet pairs required')
    axes, patches, ranks, ids = [], [], [], []
    for rank in range(1, 6):
        vertices = np.unique(np.round(display[64+(rank-1)*32:64+rank*32].reshape(-1, 3), 12), axis=0)
        # Affine-fitted eight-sided cylinders retain two rings and cap centres.
        centre = vertices.mean(axis=0)
        _, _, vh = np.linalg.svd(vertices-centre, full_matrices=False)
        direction = vh[0]
        if direction[2] < 0:
            direction = -direction
        projection = (vertices-centre)@direction
        bottom = vertices[projection < 0].mean(axis=0)
        top = vertices[projection > 0].mean(axis=0)
        if len(vertices) != 18:
            raise ValueError('Illustrated frond cylinder changed')
        d = top-bottom; unit = d/np.linalg.norm(d)
        tangent = np.cross([0., 0., 1.], unit); tangent /= np.linalg.norm(tangent)
        axes.append([bottom, top])
        for i, t in enumerate(np.linspace(parameters['bare_base_fraction'], parameters['last_leaflet_fraction'], n)):
            a = bottom+t*d
            length = parameters['minimum_leaflet_length_m']+parameters['length_amplitude_m']*np.sin(np.pi*(t-.2)/.8)**.8
            for side in (-1, 1):
                middle = a+side*.55*length*tangent
                half_width = parameters['leaflet_half_width_m']*unit
                end = a+side*length*tangent+parameters['tip_forward_m']*unit
                patches.extend([[a, middle-half_width, end], [a, end, middle+half_width]])
                ranks.extend([rank, rank]); ids.extend([2*((rank-1)*n+i)+(side+1)//2]*2)
    patches = np.asarray(patches)
    # This diagnostic deliberately excludes even the upper central stem column.
    if np.linalg.norm(patches[..., :2], axis=-1).min() <= cfg.plant['trunk_radius_m']+.05:
        raise ValueError('Representative leaflets overlap the central trunk exclusion')
    tri = np.concatenate([display, patches])
    mask = np.arange(len(tri)) >= len(display)
    tissue = np.array(['Trunk']*64+['IllustratedFrond']*160+['Leaflet']*len(patches))
    rank = np.r_[np.full(64, -1), np.repeat(np.arange(1, 6), 32), ranks]
    return dict(triangles=tri, source_mask=mask, tissue=tissue, frond_rank=rank,
        source_leaflet_id=np.asarray(ids), display_triangles=display, axes=np.asarray(axes))


def export_geometry(root, cfg, parameters, destination):
    destination = Path(destination); destination.mkdir(parents=True, exist_ok=False)
    g = build_geometry(root, cfg, parameters)
    np.savez_compressed(destination/'palm.npz', **g)
    area = triangle_area(g['triangles'][g['source_mask']])
    meta = dict(status='ILLUSTRATIVE_SOURCE_GEOMETRY_NOT_MEASURED',
        display='Exact original fitted_palm triangle coordinates; representative leaflets shown only in source-explanation figure',
        coordinate_unit='m', parameters=parameters, expanded_leaflet_nodes=len(np.unique(g['source_leaflet_id'])),
        reconstructed_one_sided_leaf_area_m2=float(area.sum()),
        height_above_floor_m=cfg.plant['height_m'], crown_radii_proxy_m=cfg.plant['crown_radii_m'],
        floor_z_m=cfg.geometry['floor_z_m'], source_triangle_count=int(g['source_mask'].sum()),
        frond_axes_m=g['axes'].tolist(), uptake='Expanded representative leaflets only; no trunk or bare-frond-axis sink',
        minimum_source_radial_distance_m=float(np.linalg.norm(g['triangles'][g['source_mask'], :, :2], axis=-1).min()),
        limitation='Five fronds, leaflet count/shape/area and source allocation are schematic assumptions. Whole-tree net exchange is prescribed; leaf-level physiology is not measured or dynamically simulated.')
    (destination/'geometry.json').write_text(json.dumps(meta, indent=2)+'\n')
    return meta
