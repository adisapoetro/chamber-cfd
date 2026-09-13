"""Independent inventory, geometry, saved-field and media checks for any inputs."""
import csv
import json
from pathlib import Path
import subprocess
import numpy as np
from PIL import Image
from .workflow import digest, write_json, verify
from .views import load_cases
from ..recreated.leaf_geometry import triangle_area


def table(p, rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def require(ok, message):
    if not ok: raise ValueError(message)


def audit(out, root, *, media=True):
    verify(out,root)
    study,cases=load_cases(out)
    require(set(cases)==set(study['cases']),'Incomplete cases')
    with np.load(out/'geometry/palm.npz') as g:
        mask=g['source_mask']; tri=g['triangles'][mask]
        np.testing.assert_array_equal(mask,g['tissue']=='Leaflet')
        require(not mask[:224].any(),'Trunk/frond axes included in source')
        area=float(triangle_area(tri).sum())
    results={}; total=0
    for name,c in cases.items():
        cfg,ts,s=c['cfg'],c['ts'],c['summary']; volume=cfg.width*cfg.depth*cfg.height
        require((c['path']/'COMPLETE.json').exists(),'Unfinished case '+name)
        for rel,sha in json.loads((c['path']/'COMPLETE.json').read_text())['files'].items():
            require(digest(c['path']/rel)==sha,'Changed case artifact: '+rel)
        require(np.isfinite(np.concatenate(list(ts.values()))).all(),'Nonfinite history')
        require(ts['time_s'][-1]==cfg.duration_s,'Incomplete time coverage')
        np.testing.assert_array_equal(ts['source_umol_s'],cfg.tree_source_umol_s)
        np.testing.assert_array_equal(ts['fans_on'],[c['fans'].active(cfg,t) for t in ts['time_s']])
        np.testing.assert_allclose(ts['roi_volume_m3'],volume,atol=1e-10,rtol=0)
        require(s['max_budget_fraction']<.001,'CO2 inventory conservation failure')
        require(s['max_divergence_s_inv']<1e-7,'Pressure projection failure')
        require(s['max_gcl_m3']<1e-8,'Moving-volume conservation failure')
        floor=cfg.geometry['floor_z_m']; error=0.; maxsource=0.
        for frame in s['frames']:
            total+=1
            with np.load(c['path']/frame['path']) as d, np.load(c['path']/'leaf_area'/Path(frame['path']).name) as q:
                require(all(np.isfinite(d[k]).all() for k in d.files),'Nonfinite field')
                require(d['co2'].min()>=0,'Negative CO2 concentration')
                require((d['cell_volume_m3']>0).all(),'Nonpositive saved cell volume')
                ix=np.flatnonzero(abs(ts['time_s']-float(d['time_s']))<1e-9)
                require(len(ix)==1 and float(d['time_s'])==frame['time_s'],'Snapshot time mismatch')
                w=[]
                for key,(lo,hi) in zip(('x_faces','y_faces','z_faces'),
                    ((-cfg.width/2,cfg.width/2),(-cfg.depth/2,cfg.depth/2),(floor,floor+cfg.height))):
                    edges=d[key]; w.append(np.maximum(0,np.minimum(edges[1:],hi)-np.maximum(edges[:-1],lo)))
                    np.testing.assert_array_equal(q[key],d[key])
                w=w[0][:,None,None]*w[1][None,:,None]*w[2][None,None,:]
                np.testing.assert_allclose(w.sum(),volume,atol=1e-10,rtol=0)
                mean=float((w*d['co2']).sum()/volume)
                sd=float(np.sqrt((w*(d['co2']-mean)**2).sum()/volume))
                error=max(error,abs(mean-ts['roi_mean_ppm'][ix[0]]),abs(sd-ts['roi_std_ppm'][ix[0]]))
                require(np.isfinite(q['leaf_area_m2']).all() and (q['leaf_area_m2']>=0).all(),'Invalid leaf map')
                np.testing.assert_allclose(q['leaf_area_m2'].sum(),area,rtol=1e-10,atol=1e-12)
                np.testing.assert_allclose(q['source_umol_s'],cfg.tree_source_umol_s*q['leaf_area_m2']/area,rtol=1e-10,atol=1e-12)
                np.testing.assert_array_equal(q['source_umol_s'][q['leaf_area_m2']==0],0.)
                maxsource=max(maxsource,abs(float(q['source_umol_s'].sum())-cfg.tree_source_umol_s))
        require(error<1e-8,'History differs from saved fields')
        moles=cfg.pressure_pa*volume/(8.31446261815324*cfg.temperature_k)
        closed=ts['time_s']<=cfg.open_at
        expected=cfg.inside_ppm+cfg.tree_source_umol_s*ts['time_s'][closed]/moles
        analytic=float(np.max(abs(ts['roi_mean_ppm'][closed]-expected)))
        require(analytic<1e-5,'Sealed mean differs from analytical mass balance')
        results[name]=dict(status='PASS',volume_m3=volume,frames=len(s['frames']),
            minimum_ppm=s['minimum_ppm'],maximum_ppm=s['maximum_ppm'],
            sealed_error_ppm=analytic,field_history_error_ppm=error,source_error_umol_s=maxsource,
            max_budget_fraction=s['max_budget_fraction'],max_divergence_s_inv=s['max_divergence_s_inv'],max_gcl_m3=s['max_gcl_m3'])
    screens={}
    for name in ('grid','dt'):
        if name not in cases: continue
        ts=cases[name]['ts']; baseline=cases['main']['ts']
        # Use actual common saved times. Do not silently interpolate diagnostics.
        sample=np.array([f['time_s'] for f in cases[name]['summary']['frames'] if f['time_s']>0])
        ia=np.searchsorted(baseline['time_s'],sample); ib=np.searchsorted(ts['time_s'],sample)
        np.testing.assert_allclose(baseline['time_s'][ia],sample,rtol=0,atol=1e-9)
        np.testing.assert_allclose(ts['time_s'][ib],sample,rtol=0,atol=1e-9)
        metrics={k:float(np.linalg.norm(ts[k][ib]-baseline[k][ia])/max(np.linalg.norm(baseline[k][ia]),1e-15))
                 for k in ('roi_std_ppm','roi_mean_speed_m_s')}
        screens[name]=dict(status='PASS' if max(metrics.values())<=study['inputs']['screens']['relative_tolerance'] else 'FAIL',
            duration_s=float(sample[-1]),relative_rms_change=metrics)
    result=dict(implementation_status='PASS',scientific_status='DIAGNOSTIC_NOT_VALIDATED',cases=results,
        sensitivity=screens,screen_scope='Configured screen interval only; first closure by default. Opening convergence is not established.',
        saved_fields=total,leaf_source_maps=total,media_status='NOT_CHECKED')
    if media:
        records=list((out/'verification/media').glob('*.json')); require(len(records)==4,'Four media views required')
        for p in records:
            r=json.loads(p.read_text()); require(r['media']['background']=='#FFFFFF','Wrong background')
            for rel,sha in r['output_hashes'].items(): require(digest(out/rel)==sha,'Changed media: '+rel)
            limits=r['media'].get('color_limits_ppm')
            if limits:
                require(limits[0]<=results['main']['minimum_ppm'] and limits[1]>=results['main']['maximum_ppm'],'Clipped CO2 colour scale')
        files=list((out/'videos').glob('*.mp4'))+list((out/'videos').glob('*.gif'))
        require(len(files)==8,'Four MP4/GIF pairs required')
        for p in files:
            raw=subprocess.check_output(['ffmpeg','-v','error','-i',str(p),'-vf','crop=8:8:0:0,format=rgb24','-f','rawvideo','-'])
            require(len(raw)>0 and np.frombuffer(raw,dtype=np.uint8).min()>=250,'Media decode/white corner failure: '+p.name)
        for p in (out/'figures').glob('*.png'):
            with Image.open(p) as im:
                a=np.asarray(im.convert('RGBA'))
                for block in (a[:8,:8],a[-8:,:8],a[:8,-8:],a[-8:,-8:]):
                    require((block==255).all(),'Nonwhite figure: '+p.name)
        result['media_status']='PASS'; result['media_files']=len(files)
    write_json(out/'verification/checks.json',result)
    table(out/'data/cases_summary.csv',[dict(case=k,**v) for k,v in results.items()])
    if screens:
        table(out/'data/numerical_sensitivity.csv',[dict(case=k,status=v['status'],duration_s=v['duration_s'],**v['relative_rms_change']) for k,v in screens.items()])
    f=cases['main']['fans']; cfg=cases['main']['cfg']
    rows=[dict(fan=i+1,x_m=point[0],y_m=point[1],height_above_floor_m=height,
        direction_x=axis[0],direction_y=axis[1],direction_z=axis[2])
        for i,(point,height,axis) in enumerate(zip(f.positions(cfg,0),f.heights_above_floor_m,f.directions()))]
    if rows: table(out/'data/fan_configuration.csv',rows)
    failures=[name for name,v in screens.items() if v['status']=='FAIL']
    (out/'README.md').write_text(f'''# {study['title']}

Geometry: **{cfg.width:g} × {cfg.depth:g} × {cfg.height:g} m**, volume **{cfg.width*cfg.depth*cfg.height:g} m³**.
Prescribed net exchange: **{cfg.tree_source_umol_s:g} µmol/s**, continuous on representative leaflets.
Negative means removal from air. Fans: **{len(f.heights_above_floor_m)}**, fixed support, closed only.

| View | Video | GIF |
|---|---|---|
| Chamber cycle | [MP4](videos/cycle.mp4) | [GIF](videos/cycle.gif) |
| Fan section | [MP4](videos/fan_sections.mp4) | [GIF](videos/fan_sections.gif) |
| Operating timeseries | [MP4](videos/operating_timeseries.mp4) | [GIF](videos/operating_timeseries.gif) |
| Leaf exchange location | [MP4](videos/leaf_source_location.mp4) | [GIF](videos/leaf_source_location.gif) |

[Figures](figures/) · [Fields and histories](data/cfd_runs/) · [Resolved inputs](configs/inputs.json).
[Methods](METHOD.md) · [Inputs and units](INPUTS.md) · [Reproduction](REPRODUCING.md).

Engineering checks: PASS. Media: {result['media_status']}.
Sensitivity failures: {', '.join(failures) or ('none in configured screens' if screens else 'NOT RUN')}.
[Detailed checks](verification/checks.json). This is a diagnostic; installed fan performance,
pole blockage, leaf geometry and physiological forcing remain conditional. A first-closure
screen does not establish opening convergence. No dynamic plant physiology is solved.
''')
    manifest=out/'verification/bundle_manifest.json'
    write_json(manifest,dict(files={str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*'))
        if p.is_file() and p!=manifest}))
    return result
