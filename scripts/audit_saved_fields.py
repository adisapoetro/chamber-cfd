#!/usr/bin/env python3
"""Read-only replay of both full reference bundles; write a separate audit JSON.

Exit 0 means the replay completed, not that the chamber is validated. The
report separates file integrity, conservation, resolution and physical evidence.
No solver or plotting module is imported and no input bundle is modified.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def verify_manifest(folder, records):
    for name, digest in records.items():
        assert sha(folder / name) == digest, name
    return len(records)


def source_identity(folder, current_root):
    records = read_json(folder / 'configs/source_hashes.json')
    for name, digest in records.items():
        assert sha(folder / 'configs/source' / name) == digest, name
    core = {name: digest for name, digest in records.items() if name.startswith('src/')}
    current_matches = {name: sha(current_root / name) == digest for name, digest in core.items()}
    return {'frozen_files_verified': len(records), 'generating_core_sha256': core,
            'current_core_matches_generating_source': all(current_matches.values()),
            'different_current_core_files': [name for name, same in current_matches.items() if not same]}


def weights(edges, bounds):
    overlap = [np.maximum(0, np.minimum(e[1:], hi) - np.maximum(e[:-1], lo))
               for e, (lo, hi) in zip(edges, bounds)]
    return overlap[0][:, None, None] * overlap[1][None, :, None] * overlap[2][None, None, :]


def python_bundle(folder, repo):
    study = read_json(folder / 'configs/study.json')
    identity = source_identity(folder, repo)
    geometry = np.load(folder / 'geometry/palm.npz')
    mask = geometry['source_mask']
    np.testing.assert_array_equal(mask, geometry['tissue'] == 'Leaflet')
    triangles = geometry['triangles'][mask]
    leaflet_area = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0],
                                          triangles[:, 2] - triangles[:, 0]), axis=1).sum() / 2
    results = {}; histories = {}; summaries = {}
    for name, item in study['cases'].items():
        cfg = item['config']; case = folder / 'data/cfd_runs' / name
        verified = verify_manifest(case, read_json(case / 'COMPLETE.json')['files'])
        summary = read_json(case / 'summary.json'); summaries[name] = summary
        history = np.genfromtxt(case / 'timeseries.csv', delimiter=',', names=True)
        histories[name] = history
        volume = cfg['width'] * cfg['depth'] * cfg['height']
        floor = cfg['geometry']['floor_z_m']
        bounds = ((-cfg['width']/2, cfg['width']/2), (-cfg['depth']/2, cfg['depth']/2),
                  (floor, floor + cfg['height']))
        mean_error = sd_error = inventory_error = source_error = 0.
        for frame in summary['frames']:
            with np.load(case / frame['path']) as data, np.load(case / 'leaf_area' / Path(frame['path']).name) as leaf:
                C = data['co2']; V = data['cell_volume_m3']; U = data['velocity']
                assert C.ndim == 3 and U.shape == C.shape + (3,)
                assert np.isfinite(C).all() and np.isfinite(U).all() and np.isfinite(V).all()
                assert C.min() >= 0 and V.min() > 0
                w = weights([data[k] for k in ('x_faces', 'y_faces', 'z_faces')], bounds)
                assert abs(w.sum() - volume) < 1e-9
                matches = np.flatnonzero(abs(history['time_s'] - float(data['time_s'])) < 1e-9)
                assert len(matches) == 1
                row = history[matches[0]]
                mean = np.sum(w*C) / volume
                sd = np.sqrt(np.sum(w*(C-mean)**2) / volume)
                mean_error = max(mean_error, abs(mean-row['roi_mean_ppm']))
                sd_error = max(sd_error, abs(sd-row['roi_std_ppm']))
                inventory_error = max(inventory_error, abs(np.sum(V*C)-row['inventory_ppm_m3']))
                assert np.all(leaf['leaf_area_m2'] >= 0)
                assert abs(leaf['leaf_area_m2'].sum()-leaflet_area) < 1e-9
                expected_source = cfg['tree_source_umol_s'] * leaf['leaf_area_m2'] / leaflet_area
                source_error = max(source_error, float(np.max(abs(leaf['source_umol_s']-expected_source))))
        closed = history['time_s'] <= cfg['open_at']
        molar = cfg['pressure_pa'] / (8.31446261815324 * cfg['temperature_k'])
        expected = cfg['inside_ppm'] + cfg['tree_source_umol_s'] * history['time_s'][closed] / (molar*volume)
        closed_error = float(np.max(abs(history['roi_mean_ppm'][closed]-expected)))
        assert max(mean_error, sd_error, source_error) < 1e-8
        assert inventory_error < 1e-6 and closed_error < 1e-5
        results[name] = {'verified_case_files': verified, 'snapshots': len(summary['frames']),
            'reference_volume_m3': volume, 'saved_mean_error_ppm': float(mean_error),
            'saved_sd_error_ppm': float(sd_error), 'saved_inventory_history_error_ppm_m3': float(inventory_error),
            'source_distribution_error_umol_s': float(source_error), 'sealed_analytic_error_ppm': closed_error,
            'recorded_budget_fraction_max': float(np.max(abs(history['budget_fraction'])))}
    sensitivity = {}
    for name in ('grid', 'dt'):
        if name not in histories:
            continue
        sample = np.array([f['time_s'] for f in summaries[name]['frames'] if f['time_s'] > 0])
        a, b = histories['main'], histories[name]
        ia, ib = np.searchsorted(a['time_s'], sample), np.searchsorted(b['time_s'], sample)
        np.testing.assert_allclose(a['time_s'][ia], sample, rtol=0, atol=1e-9)
        np.testing.assert_allclose(b['time_s'][ib], sample, rtol=0, atol=1e-9)
        changes = {k: float(np.linalg.norm(b[k][ib]-a[k][ia])/np.linalg.norm(a[k][ia]))
                   for k in ('roi_std_ppm', 'roi_mean_speed_m_s')}
        threshold = study['inputs']['screens']['relative_tolerance']
        sensitivity[name] = {'relative_rms_changes': changes, 'threshold': threshold,
                             'end_time_s': float(sample[-1]),
                             'status': 'PASS' if max(changes.values()) <= threshold else 'FAIL'}
    return {'backend': 'Project-owned NumPy/SciPy ALE finite-volume solver',
        'source_identity': identity, 'cases': results, 'sensitivity': sensitivity,
        'leaflet_triangles': len(triangles), 'leaflet_proxy_area_m2': float(leaflet_area),
        'independent_open_phase_inventory_budget': 'NOT_REPLAYABLE: boundary scalar flux is not saved separately',
        'opening_grid_convergence': 'NOT_ESTABLISHED', 'experimental_validation': 'NOT_ESTABLISHED'}


def waterlily_bundle(folder, repo):
    summary = read_json(folder / 'summary.json'); config = read_json(folder / 'inputs.json')
    count = verify_manifest(folder, read_json(folder / 'manifest.json')['sha256'])
    identity = source_identity(folder, repo / 'waterlily')
    shape = tuple(summary['shape']); dx = summary['grid_cell_m']; lower = np.array(summary['grid_lower_m'])
    ch, air = config['chamber'], config['air']
    bounds = ((-ch['width_m']/2, ch['width_m']/2), (-ch['depth_m']/2, ch['depth_m']/2),
              (ch['floor_m'], ch['floor_m']+ch['height_m']))
    w = weights([lower[k]+np.arange(shape[k]+1)*dx for k in range(3)], bounds)
    volume = ch['width_m']*ch['depth_m']*ch['height_m']
    assert abs(w.sum()-volume) < 1e-9
    audit = np.genfromtxt(folder / 'transport_audit.csv', delimiter=',', names=True)
    molar = air['pressure_pa']/(8.31446261815324*air['temperature_k'])
    rate = config['exchange']['net_co2_umol_s']
    assert np.max(abs(audit['source_ppm_m3_s']-rate/molar)) < 1e-10
    integral = np.cumsum(audit['dt_s']*(audit['source_ppm_m3_s']-audit['boundary_outward_ppm_m3_s']))
    initial = None; budget_error = mean_error = sd_error = closed_error = 0.
    for row in summary['snapshots']:
        def array(name, vector=False):
            value = np.fromfile(folder / row[name], dtype='<f8').reshape(shape+((3,) if vector else ()), order='F')
            assert np.isfinite(value).all()
            return value
        C, V, U = array('co2_file'), array('volume_file'), array('velocity_file', True)
        assert C[V > 0].min() >= 0 and V.min() >= 0 and V.max() <= dx**3+1e-12
        inventory = float(np.sum(C*V))
        if initial is None:
            initial = inventory
        i = np.searchsorted(audit['time_s'], row['time_s']+1e-8)-1
        budget_error = max(budget_error, abs(inventory-initial-(integral[i] if i >= 0 else 0)))
        mean = float(np.sum(C*w)/volume)
        sd = float(np.sqrt(np.sum(w*(C-mean)**2)/volume))
        mean_error = max(mean_error, abs(mean-row['co2_mean_ppm']))
        sd_error = max(sd_error, abs(sd-row['co2_spatial_sd_ppm']))
        if row['time_s'] <= config['cycle']['closed_s']+1e-8:
            exact = air['initial_co2_ppm']+rate*row['time_s']/(molar*volume)
            closed_error = max(closed_error, abs(mean-exact))
    assert budget_error < 1e-5 and max(mean_error, sd_error, closed_error) < 1e-8
    worst = int(np.argmax(audit['roi_correction_relative_l2']))
    relative = float(audit['roi_correction_relative_l2'][worst])
    return {'backend': 'WaterLily '+summary['version']+' airflow plus project-owned Julia scalar extension',
        'source_identity': identity, 'manifest_files_verified': count, 'shape': list(shape),
        'snapshots': len(summary['snapshots']), 'transport_steps': len(audit),
        'independent_inventory_error_ppm_m3': float(budget_error), 'saved_mean_error_ppm': float(mean_error),
        'saved_sd_error_ppm': float(sd_error), 'sealed_analytic_error_ppm': float(closed_error),
        'reference_volume_m3': volume, 'numerical_wall_m': summary['numerical_wall_thickness_m'],
        'physical_wall_input_m': ch['physical_wall_thickness_m'],
        'max_roi_flux_correction_fraction': relative, 'max_roi_flux_correction_time_s': float(audit['time_s'][worst]),
        'max_aperture_velocity_correction_m_s': float(np.max(audit['correction_max_m_s'])),
        'max_gcl_residual_m3_s': float(np.max(audit['gcl_max_m3_s'])),
        'coupling_screen': 'PASS' if relative <= 0.1 else 'FAIL',
        'full_cycle_grid_convergence': 'NOT_ESTABLISHED', 'experimental_validation': 'NOT_ESTABLISHED',
        'budget_scope': 'Saved concentration/volume replay against logged source and boundary flux; native face-flux reconstruction is not independently replayed'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python-run', required=True, type=Path)
    parser.add_argument('--waterlily-run', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Choose a fresh audit output')
    repo = Path(__file__).resolve().parents[1]
    report = {'scope': 'Read-only numerical replay; no field validation or changes to existing results',
              'script_sha256': sha(Path(__file__)),
              'python': python_bundle(args.python_run, repo),
              'waterlily': waterlily_bundle(args.waterlily_run, repo),
              'implementation_replay': 'PASS', 'physical_qualification': 'NOT_ESTABLISHED'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
