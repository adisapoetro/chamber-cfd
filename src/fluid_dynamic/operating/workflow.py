"""Freeze, run, render and audit a self-contained input configuration."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import redirect_stdout, redirect_stderr
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import numpy as np
from .config import load, compile_study
from ..scenarios.geometry import ScenarioStudy
from ..recreated.fan_forcing import FanArray
from ..recreated.fan_solver import run_case
from ..recreated.frond_geometry import export_geometry
from ..recreated.leaf_source import mesh_factory


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write_json(p, data):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def source_paths(root):
    paths = list((root/'src/fluid_dynamic/operating').glob('*.py'))
    paths += [root/'src/fluid_dynamic/recreated'/n for n in (
        '__init__.py','mesh.py','operators.py','solver.py','fan_forcing.py','fan_solver.py',
        'leaf_geometry.py','leaf_source.py','frond_geometry.py','render.py')]
    paths += [root/'src/fluid_dynamic/scenarios'/n for n in ('__init__.py','geometry.py','presentation_render.py')]
    paths += [root/p for p in ('src/fluid_dynamic/__init__.py','meshes/palm_libz.obj',
        'scripts/operating_chamber.py','requirements-operating.lock','configs/operating/current.yaml',
        'docs/operating/METHOD.md','docs/operating/INPUTS.md','docs/operating/REPRODUCING.md')]
    return sorted(paths)


def prepare(root, out, inputs):
    study = compile_study(inputs)
    if out.exists():
        raise ValueError('Output already exists. Use a new directory; use run/render/audit to resume frozen inputs.')
    # Validate source geometry before creating any output; failed inputs leave no partial bundle.
    from ..recreated.frond_geometry import build_geometry
    cfg = ScenarioStudy(**study['cases']['main']['config'])
    g = build_geometry(root, cfg, study['geometry'])
    tri = g['triangles']
    if (abs(tri[...,0]).max() >= cfg.width/2 or abs(tri[...,1]).max() >= cfg.depth/2
        or tri[...,2].max() >= cfg.geometry['floor_z_m']+cfg.height):
        raise ValueError('Leaflet geometry does not fit within the closed chamber')
    for folder in ('configs/source','data/cfd_runs','figures','videos','verification','logs'):
        (out/folder).mkdir(parents=True, exist_ok=True)
    write_json(out/'configs/inputs.json', inputs.model_dump())
    write_json(out/'configs/study.json', study)
    export_geometry(root, cfg, study['geometry'], out/'geometry')
    for p in source_paths(root):
        dest = out/'configs/source'/p.relative_to(root); dest.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(p,dest)
    write_json(out/'configs/source_hashes.json', {str(p.relative_to(root)):digest(p) for p in source_paths(root)})
    frozen = ['configs/inputs.json','configs/study.json','geometry/palm.npz','geometry/geometry.json']
    write_json(out/'configs/input_hashes.json', {p:digest(out/p) for p in frozen})
    write_json(out/'configs/environment.json', dict(python=sys.version,platform=platform.platform(),
        packages={p:version(p) for p in ('numpy','scipy','matplotlib','pydantic','PyYAML','Pillow')},
        blas_threads=1, ffmpeg=subprocess.check_output(['ffmpeg','-version'],text=True).splitlines()[0]))
    for name in ('METHOD.md','INPUTS.md','REPRODUCING.md'):
        shutil.copy2(root/'docs/operating'/name,out/name)
    return study


def verify(out, root):
    for p,sha in json.loads((out/'configs/input_hashes.json').read_text()).items():
        if digest(out/p) != sha:
            raise ValueError('Frozen input changed: '+p)
    for p,sha in json.loads((out/'configs/source_hashes.json').read_text()).items():
        if digest(root/p) != sha or digest(out/'configs/source'/p) != sha:
            raise ValueError('Code differs from frozen run: '+p+'; use a fresh output or its frozen source')


def execute(job):
    out, name, item = job; out=Path(out); path=out/'data/cfd_runs'/name
    marker=path/'COMPLETE.json'
    if marker.exists():
        record=json.loads(marker.read_text())
        for p,sha in record['files'].items():
            if digest(path/p)!=sha: raise ValueError('Completed case changed: '+name+'/'+p)
        return name, 'verified existing case'
    if path.exists():
        raise ValueError('Incomplete case exists: '+str(path)+'. Preserve it and retry with a new output.')
    with (out/'logs'/(name+'.log')).open('w',buffering=1) as log, redirect_stdout(log), redirect_stderr(log):
        cfg=ScenarioStudy(**item['config']); factory=mesh_factory(cfg,out/item['source_geometry'])
        result=run_case(cfg,path,mesh_factory=factory,fans=FanArray(**item['fans']))
        source=dict(model='Prescribed net exchange on representative leaflets',total_umol_s=cfg.tree_source_umol_s,
            continuous_open_and_closed=True, excluded='Trunk and bare frond axes', maps=[])
        (path/'leaf_area').mkdir()
        for frame in result['frames']:
            with np.load(path/frame['path']) as d:
                m=factory(cfg,float(d['gap'])); area=m.source_weights_at(float(d['gap'])).reshape(m.shape)
                rel='leaf_area/'+Path(frame['path']).name
                np.savez_compressed(path/rel,time_s=d['time_s'],gap=d['gap'],leaf_area_m2=area,
                    source_umol_s=cfg.tree_source_umol_s*area/area.sum(),
                    x_faces=d['x_faces'],y_faces=d['y_faces'],z_faces=d['z_faces'])
                source['maps'].append(dict(time_s=frame['time_s'],path=rel,area_m2=float(area.sum())))
        write_json(path/'leaf_source.json',source)
        write_json(marker,dict(files={str(p.relative_to(path)):digest(p) for p in path.rglob('*') if p.is_file()}))
    return name,dict(wall_seconds=result['wall_seconds'],max_budget_fraction=result['max_budget_fraction'])


def run(out, study, workers):
    jobs=[(str(out),name,item) for name,item in study['cases'].items()]
    if workers==1:
        for job in jobs: print(execute(job),flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for future in as_completed([pool.submit(execute,j) for j in jobs]): print(future.result(),flush=True)


def main(root=None):
    root=Path(root or Path(__file__).resolve().parents[3])
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['validate','prepare','run','render','audit','all'])
    p.add_argument('--input',type=Path,default=root/'configs/operating/current.yaml')
    p.add_argument('--set',action='append',default=[],metavar='GROUP.KEY=VALUE')
    p.add_argument('--output',type=Path,help='Fresh output required for prepare/all; frozen bundle for run/render/audit')
    p.add_argument('--workers',type=int,default=2)
    p.add_argument('--preview',action='store_true',help='Render still previews only; not a complete media bundle')
    args=p.parse_args()
    if args.workers<1: p.error('workers must be positive')
    if args.action=='validate':
        print(json.dumps(compile_study(load(args.input,args.set)),indent=2)); return
    if args.output is None: p.error('--output is required')
    out=args.output.resolve()
    if args.action in ('prepare','all'):
        prepare(root,out,load(args.input,args.set))
    elif args.set or args.input != root/'configs/operating/current.yaml':
        p.error('run/render/audit read frozen configs; set inputs only with prepare/all')
    verify(out,root); study=json.loads((out/'configs/study.json').read_text())
    if args.action in ('run','all'): run(out,study,args.workers)
    if args.action in ('render','all'):
        from .render import render_bundle
        render_bundle(out,root,preview=args.preview)
    if args.action in ('audit','all'):
        from .audit import audit
        print(json.dumps(audit(out,root,media=not args.preview),indent=2))


if __name__=='__main__':
    main()
