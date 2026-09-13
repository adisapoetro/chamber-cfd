"""Generate a synthetic juvenile oil-palm mesh and save as OBJ.

The result is a watertight triangle mesh approximating:
  - a trunk (cylinder)
  - 5 radial fronds (cylinders tilted ~30° from vertical)
  - a small fruit bunch (sphere) near the trunk top

Coordinate frame: x and y centered on origin; z=0 = floor; z+ = up.
Total height ~4 m, fits within a 4×4 m chamber footprint.

Usage:
    python scripts/build_synthetic_palm.py
"""

from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import trimesh


def build_palm() -> trimesh.Trimesh:
    parts: list[trimesh.Trimesh] = []

    # Trunk: cylinder of radius 0.15 m, height 1.8 m, base at z=0.
    trunk = trimesh.creation.cylinder(radius=0.15, height=1.8, sections=16)
    trunk.apply_translation([0.0, 0.0, 0.9])  # center at z=0.9 so base sits at 0
    parts.append(trunk)

    # 5 fronds: cylinders rooted at trunk top (z=1.8), tilted 30° from vertical, length 1.5 m.
    n_fronds = 5
    frond_length = 1.5
    tilt_deg = 30.0
    for i in range(n_fronds):
        azimuth = 2 * np.pi * i / n_fronds
        # Build a cylinder along +z, then rotate by tilt around y, then by azimuth around z.
        f = trimesh.creation.cylinder(radius=0.05, height=frond_length, sections=8)
        # Move so the bottom is at origin, then translate to trunk top.
        f.apply_translation([0.0, 0.0, frond_length / 2])
        # Tilt: rotate about y-axis by tilt_deg.
        tilt_rad = np.deg2rad(tilt_deg)
        Ry = trimesh.transformations.rotation_matrix(tilt_rad, [0, 1, 0])
        f.apply_transform(Ry)
        # Azimuth: rotate about z-axis.
        Rz = trimesh.transformations.rotation_matrix(azimuth, [0, 0, 1])
        f.apply_transform(Rz)
        # Translate to trunk top.
        f.apply_translation([0.0, 0.0, 1.8])
        parts.append(f)

    # Fruit bunch: a sphere near the trunk top, at z=1.8 just south of trunk axis.
    bunch = trimesh.creation.icosphere(radius=0.20, subdivisions=2)
    bunch.apply_translation([0.0, -0.30, 1.8])
    parts.append(bunch)

    # Boolean union would be cleanest but is slow and requires extra deps; concatenate is fine
    # for voxelization purposes (each part is independently watertight, so contains() will
    # report a point as inside if it's inside any part).
    palm = trimesh.util.concatenate(parts)
    return palm


def main() -> int:
    palm = build_palm()
    out_path = Path(__file__).resolve().parent.parent / "meshes" / "palm_libz.obj"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    palm.export(out_path)
    print(f"Wrote {out_path}")
    print(f"  vertices: {len(palm.vertices)}")
    print(f"  faces:    {len(palm.faces)}")
    print(f"  bounds:   {palm.bounds.tolist()}")
    print(f"  watertight: {palm.is_watertight}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
