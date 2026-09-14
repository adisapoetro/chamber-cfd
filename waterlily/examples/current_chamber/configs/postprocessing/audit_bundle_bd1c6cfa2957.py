#!/usr/bin/env python3
"""Independently replay saved CO₂ inventory and verify the generated media."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from PIL import Image


def audit(folder):
    summary=json.loads((folder/'summary.json').read_text())
    config=json.loads((folder/'inputs.json').read_text())
    assert summary['co2_transport_enabled'] and summary['engine']=='WaterLily'
    shape=tuple(summary['shape']);dx=summary['grid_cell_m'];lower=np.array(summary['grid_lower_m'])
    ch=config['chamber'];air=config['air'];records=summary['snapshots']
    roi=[]
    for k,(lo,hi) in enumerate([(-ch['width_m']/2,ch['width_m']/2),(-ch['depth_m']/2,ch['depth_m']/2),
                              (ch['floor_m'],ch['floor_m']+ch['height_m'])]):
        cell=lower[k]+np.arange(shape[k])*dx
        roi.append(np.maximum(0,np.minimum(cell+dx,hi)-np.maximum(cell,lo)))
    weights=roi[0][:,None,None]*roi[1][None,:,None]*roi[2][None,None,:]
    volume=ch['width_m']*ch['depth_m']*ch['height_m'];assert abs(weights.sum()-volume)<1e-9
    with (folder/'transport_audit.csv').open() as f:
        steps=list(csv.DictReader(f))
    t=np.array([float(r['time_s']) for r in steps]);dt=np.array([float(r['dt_s']) for r in steps])
    src=np.array([float(r['source_ppm_m3_s']) for r in steps])
    outward=np.array([float(r['boundary_outward_ppm_m3_s']) for r in steps])
    balance=np.cumsum(dt*(src-outward));initial=None;errors=[];sealed_errors=[]
    max_vertical_velocity=0.
    molar=air['pressure_pa']/(8.31446261815324*air['temperature_k'])
    assert np.max(abs(src-config['exchange']['net_co2_umol_s']/molar))<1e-10
    def array(name):
        data=np.fromfile(folder/name,dtype='<f8')
        assert data.size==np.prod(shape)
        return data.reshape(shape,order='F')
    for r in records:
        C=array(r['co2_file']);V=array(r['volume_file'])
        velocity=np.fromfile(folder/r['velocity_file'],dtype='<f8')
        assert velocity.size==3*np.prod(shape)
        velocity=velocity.reshape(shape+(3,),order='F')
        assert np.isfinite(velocity).all()
        max_vertical_velocity=max(max_vertical_velocity,float(np.max(abs(velocity[...,2][V>0]))))
        assert np.isfinite(C[V>0]).all() and C[V>0].min()>=-1e-9
        assert V.min()>=0 and V.max()<=dx**3+1e-12
        inventory=float(np.sum(C*V));mean=float(np.sum(C*weights)/volume)
        assert abs(mean-r['co2_mean_ppm'])<1e-9
        if initial is None: initial=inventory
        index=np.searchsorted(t,r['time_s']+1e-8)-1
        expected=initial+(balance[index] if index>=0 else 0)
        errors.append(abs(inventory-expected))
        if r['time_s']<=config['cycle']['closed_s']+1e-8:
            analytic=air['initial_co2_ppm']+config['exchange']['net_co2_umol_s']*r['time_s']/(molar*volume)
            sealed_errors.append(abs(mean-analytic))
    assert max(errors)<1e-5
    assert max(sealed_errors)<1e-8
    assert abs(records[-1]['time_s']-summary['duration_s'])<1e-8
    rendering=json.loads((folder/'rendering.json').read_text());media=rendering['media']
    pngs=[p for p in media if p.endswith('.png')];movies=[p for p in media if p.endswith('.mp4')]
    views=rendering.get('views',['c2_202502_cycle','fan_sections','operating_timeseries','leaf_source_location'])
    cycle=config['cycle'];closed=cycle['closed_s']
    still_times=[0,closed,closed+cycle['opening_travel_s']/2,closed+cycle['opening_travel_s'],
                 closed+cycle['open_phase_s'],records[-1]['time_s']]
    saved_time=np.array([r['time_s'] for r in records])
    still_indices={int(np.argmin(abs(saved_time-v))) for v in still_times}
    expected={f'videos/{view}.{ext}' for view in views for ext in ['mp4','gif']}
    expected.update(f"figures/{view}_{records[i]['time_s']:06.1f}s.png" for view in views for i in still_indices)
    assert set(media)==expected and len(media)==len(expected)
    for name in pngs:
        im=np.asarray(Image.open(folder/name).convert('RGB'))
        assert all(np.min(im[y,x])==255 for y,x in [(0,0),(0,-1),(-1,0),(-1,-1)])
    for name in movies:
        p=subprocess.run(['ffprobe','-v','error','-count_frames','-select_streams','v:0','-show_entries',
                          'stream=nb_read_frames','-of','json',str(folder/name)],check=True,capture_output=True,text=True)
        assert int(json.loads(p.stdout)['streams'][0]['nb_read_frames'])==len(records)
    for name in media:
        subprocess.run(['ffmpeg','-v','error','-xerror','-i',str(folder/name),'-f','null','-'],check=True,stdout=subprocess.DEVNULL)
    # Decode all four video/GIF corners at every frame. Allow only a small lossy
    # codec deviation from white; PNG corners above must be exactly white.
    corner_min=255
    corner_filter=('[0:v]format=rgb24,split=4[a][b][c][d];'
                   '[a]crop=2:2:0:0[a0];[b]crop=2:2:iw-2:0[b0];'
                   '[c]crop=2:2:0:ih-2[c0];[d]crop=2:2:iw-2:ih-2[d0];'
                   '[a0][b0][c0][d0]hstack=inputs=4')
    for name in [p for p in media if p.endswith(('.mp4','.gif'))]:
        corners=subprocess.run(['ffmpeg','-v','error','-i',str(folder/name),'-filter_complex',corner_filter,
                                '-f','rawvideo','-pix_fmt','rgb24','pipe:1'],check=True,capture_output=True)
        pixels=np.frombuffer(corners.stdout,dtype=np.uint8)
        assert pixels.size>0 and pixels.size%48==0
        corner_min=min(corner_min,int(pixels.min()))
    assert corner_min>=250,corner_min
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=json.loads((folder/'manifest.json').read_text())['sha256']
    for name,digest in manifest.items(): assert sha(folder/name)==digest,name
    auditor=Path(__file__).read_bytes();auditor_sha=hashlib.sha256(auditor).hexdigest()
    auditor_path=folder/f'configs/postprocessing/audit_bundle_{auditor_sha[:12]}.py'
    auditor_path.parent.mkdir(parents=True,exist_ok=True);auditor_path.write_bytes(auditor)
    result={'implementation_checks':'PASS','scientific_status':'DIAGNOSTIC_NOT_VALIDATED',
            'snapshots':len(records),'media_files':len(media),'pngs':len(pngs),'mp4s':len(movies),
            'independent_inventory_error_max_ppm_m3':max(errors),
            'first_sealed_mean_error_max_ppm':max(sealed_errors),
            'coupling_correction_10pct_gate':summary['checks']['coupling_correction_10pct_gate'],
            'source_and_output_checksums':'PASS','auditor_sha256':auditor_sha,
            'auditor_source':str(auditor_path.relative_to(folder)),
            'media_decode':'PASS','white_png_corners':'PASS',
            'white_video_gif_corners':'PASS','video_gif_corner_min_rgb':corner_min,
            'spatial_dimensions':3,'native_velocity_components':3,
            'saved_fluid_vertical_velocity_max_m_s':max_vertical_velocity}
    (folder/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    hashes={str(p.relative_to(folder)):sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}
    (folder/'manifest.json').write_text(json.dumps({'sha256':hashes},indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path)
    audit(parser.parse_args().folder)
