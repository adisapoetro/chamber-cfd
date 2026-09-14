#!/usr/bin/env python3
"""Compare complete bundles with identical inputs; retain failures and regressions."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def read(folder):
    summary=json.loads((folder/'summary.json').read_text())
    with (folder/'transport_audit.csv').open() as f:
        rows=list(csv.DictReader(f))
    values={key:np.array([float(row[key]) for row in rows]) for key in rows[0]}
    return summary,values


def main(before,after,output):
    if output.exists():
        raise FileExistsError(f'Choose a fresh comparison folder: {output}')
    assert json.loads((before/'inputs.json').read_text())==json.loads((after/'inputs.json').read_text())
    old,a=read(before);new,b=read(after)
    assert old['shape']==new['shape'] and old['duration_s']==new['duration_s']
    assert len(old['snapshots'])==len(new['snapshots'])
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    unchanged=True
    differences=[]
    for x,y in zip(old['snapshots'],new['snapshots']):
        assert x['time_s']==y['time_s']
        for field in ['velocity_file','distance_file']:
            unchanged &= sha(before/x[field])==sha(after/y[field])
        differences.append({
            'time_s':x['time_s'],
            'mean_difference_ppm':y['co2_mean_ppm']-x['co2_mean_ppm'],
            'spatial_sd_difference_ppm':y['co2_spatial_sd_ppm']-x['co2_spatial_sd_ppm']})
    report={
        'scientific_status':'DIAGNOSTIC_NOT_VALIDATED',
        'identical_inputs':True,'native_velocity_and_geometry_bitwise_identical':bool(unchanged),
        'before_summary_sha256':sha(before/'summary.json'),'after_summary_sha256':sha(after/'summary.json'),
        'comparator_sha256':sha(Path(__file__)),
        'before_checks':old['checks'],'after_checks':new['checks'],
        'snapshot_differences':differences,
        'max_abs_mean_difference_ppm':max(abs(r['mean_difference_ppm']) for r in differences),
        'max_abs_spatial_sd_difference_ppm':max(abs(r['spatial_sd_difference_ppm']) for r in differences),
        'worst_steps':{},'endpoints':[]}
    for name,v in [('before',a),('after',b)]:
        i=int(np.argmax(v['roi_correction_relative_l2']))
        report['worst_steps'][name]={key:float(v[key][i]) for key in
            ['time_s','dt_s','roi_correction_relative_l2','correction_max_m_s','raw_max_m_s']}
    for t in [300,320,900,1200,1220,1800]:
        if t>new['duration_s']:continue
        values={'time_s':t}
        for name,v in [('before',a),('after',b)]:
            i=int(np.argmin(abs(v['time_s']-t)))
            assert abs(v['time_s'][i]-t)<1e-8
            values[name+'_roi_correction_percent']=float(v['roi_correction_relative_l2'][i]*100)
        report['endpoints'].append(values)
    plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white',
                         'font.family':'DejaVu Sans','font.size':10})
    fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
    fig.suptitle('WaterLily chamber · iteration audit\nBoth versions remain diagnostic',fontsize=17)
    for name,v,s,col in [('Iteration 1',a,old,'#7A8191'),('Iteration 2',b,new,'#007F71')]:
        axes[0,0].plot(v['time_s'],100*v['roi_correction_relative_l2'],color=col,lw=1,label=name)
        axes[0,1].plot(v['time_s'],v['correction_max_m_s'],color=col,lw=1,label=name)
        axes[1,0].plot([r['time_s'] for r in s['snapshots']],
                       [r['co2_mean_ppm'] for r in s['snapshots']],color=col,lw=1.5,label=name)
    axes[0,0].axhline(10,color='#B06518',ls='--',lw=1,label='Unchanged 10% screen')
    axes[0,0].set(title='Chamber-region flux correction',ylabel='Relative correction (%)',yscale='log')
    axes[0,0].legend(frameon=False,fontsize=8)
    axes[0,1].set(title='Correction in the sharp aperture',ylabel='Velocity correction (m/s)',yscale='log')
    axes[0,1].text(.02,.94,'Coupling diagnostic; not predicted air speed',transform=axes[0,1].transAxes,
                   va='top',fontsize=8,color='#925B1B')
    axes[1,0].set(title='Mean CO₂ across the fixed 96 m³ region',ylabel='CO₂ (ppm)')
    axes[1,1].plot([r['time_s'] for r in differences],[r['mean_difference_ppm'] for r in differences],
                   color='#007F71',lw=1.5,label='Mean')
    axes[1,1].plot([r['time_s'] for r in differences],[r['spatial_sd_difference_ppm'] for r in differences],
                   color='#735AA2',lw=1,label='Spatial SD')
    axes[1,1].set(title='Iteration 2 minus iteration 1',ylabel='Difference (ppm)')
    axes[1,1].legend(frameon=False)
    for ax in axes.flat:
        ax.set_xlabel('Physical time (s)');ax.grid(alpha=.15)
        ax.spines[['top','right']].set_visible(False)
    output.mkdir(parents=True)
    fig.savefig(output/'comparison.png',dpi=140,facecolor='white');plt.close(fig)
    (output/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['snapshot_differences','before_checks','after_checks']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['before','after','output']:parser.add_argument(name,type=Path)
    args=parser.parse_args();main(args.before,args.after,args.output)
