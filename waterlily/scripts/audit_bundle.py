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
    molar=air['pressure_pa']/(8.31446261815324*air['temperature_k'])
    assert np.max(abs(src-config['exchange']['net_co2_umol_s']/molar))<1e-10
    def array(name):
        data=np.fromfile(folder/name,dtype='<f8')
        assert data.size==np.prod(shape)
        return data.reshape(shape,order='F')
    for r in records:
        C=array(r['co2_file']);V=array(r['volume_file'])
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
    media=json.loads((folder/'rendering.json').read_text())['media']
    pngs=[p for p in media if p.endswith('.png')];movies=[p for p in media if p.endswith('.mp4')]
    if summary['duration_s']==1800: assert len(pngs)==24 and len(movies)==4 and len(media)==32
    for name in pngs:
        im=np.asarray(Image.open(folder/name).convert('RGB'))
        assert all(np.min(im[y,x])==255 for y,x in [(0,0),(0,-1),(-1,0),(-1,-1)])
    for name in movies:
        p=subprocess.run(['ffprobe','-v','error','-count_frames','-select_streams','v:0','-show_entries',
                          'stream=nb_read_frames','-of','json',str(folder/name)],check=True,capture_output=True,text=True)
        assert int(json.loads(p.stdout)['streams'][0]['nb_read_frames'])==len(records)
    for name in media:
        subprocess.run(['ffmpeg','-v','error','-xerror','-i',str(folder/name),'-f','null','-'],check=True,stdout=subprocess.DEVNULL)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=json.loads((folder/'manifest.json').read_text())['sha256']
    for name,digest in manifest.items(): assert sha(folder/name)==digest,name
    result={'implementation_checks':'PASS','scientific_status':'DIAGNOSTIC_NOT_VALIDATED',
            'snapshots':len(records),'media_files':len(media),'pngs':len(pngs),'mp4s':len(movies),
            'independent_inventory_error_max_ppm_m3':max(errors),
            'first_sealed_mean_error_max_ppm':max(sealed_errors),
            'coupling_correction_10pct_gate':summary['checks']['coupling_correction_10pct_gate'],
            'source_and_output_checksums':'PASS','media_decode':'PASS','white_png_corners':'PASS'}
    (folder/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    hashes={str(p.relative_to(folder)):sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}
    (folder/'manifest.json').write_text(json.dumps({'sha256':hashes},indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path)
    audit(parser.parse_args().folder)
