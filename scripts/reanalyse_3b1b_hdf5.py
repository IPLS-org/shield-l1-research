#!/usr/bin/env python3
"""Independent raw-HDF5 reread of the 2026-11-19 PHARE dispersion benchmark.

No PHARE import, no upstream frequency analysis imported. This does not run plasma.
Input: extracted evidence ZIP directory with shield-evidence/ and PHARE/build/...
Output: independent-frequency-results.csv, independent-evaluation.json, baseline-errors.png
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import h5py

RUNS = {
    'low_k': dict(cells=128, dx=0.8, modes=(1, 2, 4, 8)),
    'high_k': dict(cells=64, dx=0.1, modes=(1, 2, 3, 4)),
}

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def read_grid(h5file, cells):
    times, buffers = [], []
    with h5py.File(h5file, 'r') as f:
        for tkey in sorted(f['t'], key=float):
            patches = f['t'][tkey]['pl0']
            arr = np.full(cells, complex(float('nan'), float('nan')), dtype=complex)
            coverage = np.zeros(cells, dtype=int)
            for patch in patches.values():
                start = int(patch.attrs['lower'][0])
                n = int(patch.attrs['nbrCells'][0])
                by = patch['EM_B_y'][...]
                bz = patch['EM_B_z'][...]
                ng = (len(by) - n) // 2
                assert len(by) == len(bz) and (len(by) - n) % 2 == 0
                arr[start:start+n] = by[ng:ng+n] + 1j * bz[ng:ng+n]
                coverage[start:start+n] += 1
            assert np.all(coverage == 1), 'grid incomplete or overlapping patches'
            times.append(float(tkey))
            buffers.append(arr)
    t = np.asarray(times)
    assert np.all(np.diff(t) > 0), 'non-increasing times'
    return t, np.stack(buffers)

def cold(k, branch):
    return k * (math.sqrt(1+k*k/4) + (1 if branch=='whistler' else -1)*k/2)

def scheme(k, branch, dx, order):
    kgrid = 2 / dx * math.sin(k*dx/2)
    deff = 1 - (order+1)*dx*dx/12  # Empirical PHARE source reference, not physics.
    return kgrid * (math.sqrt(1+(kgrid*deff)**2/4) + (1 if branch=='whistler' else -1)*kgrid*deff/2)

def evaluate(root, output):
    output.mkdir(parents=True, exist_ok=True)
    with (root/'shield-evidence/dispersion-frequency-results.csv').open(newline='') as f:
        reported = {(int(r['order']), r['run'], r['branch'], int(r['mode'])):r for r in csv.DictReader(f)}
    records=[]
    anomalies=[]
    seen=[]
    for order in (1,2,3):
        for name, p in RUNS.items():
            path = root/f'PHARE/build/tests/functional/dispersion/dispersion_{name}_o{order}/EM_B.h5'
            t, b = read_grid(path, p['cells'])
            fft = np.fft.fft(b, axis=1) / p['cells']
            for branch in ('aic','whistler'):
                for mode in p['modes']:
                    key=(order,name,branch,mode)
                    seen.append(key)
                    idx = mode if branch == 'aic' else -mode
                    coeff = fft[:,idx]
                    phase = np.unwrap(np.angle(coeff))
                    # Unlike PHARE's direct mode projection this uses an FFT, independently.
                    slope = np.polyfit(t, phase, 1)[0]
                    measured = abs(float(slope))
                    ratio=float(abs(coeff[-1])/abs(coeff[0]))
                    k=2*math.pi*mode/(p['cells']*p['dx'])
                    w_cold=cold(k,branch)
                    w_scheme=scheme(k,branch,p['dx'],order)
                    r=reported[key]
                    diff=abs(measured-float(r['omega_measured_printed']))
                    if diff>5.1e-7:
                        anomalies.append(f'frequency differs from printed value: {key} {diff}')
                    if not (0.3<ratio<3.0):
                        anomalies.append(f'amplitude outside prior window: {key}: {ratio}')
                    records.append({
                        'order':order,'run':name,'branch':branch,'mode':mode,
                        'cells':p['cells'],'dx_di':p['dx'],'k_di_exact':k,
                        'omega_fft_reread':measured,
                        'omega_reported_rounded':float(r['omega_measured_printed']),
                        'abs_frequency_difference':diff,
                        'omega_cold':w_cold,'omega_scheme_empirical':w_scheme,
                        'rel_error_vs_cold':(measured-w_cold)/w_cold,
                        'rel_error_vs_scheme_empirical':(measured-w_scheme)/w_scheme,
                        'amp_end_over_start':ratio,
                        'phase_fit_residual_rms_radians':float(np.sqrt(np.mean((phase-np.polyval([slope, np.polyfit(t,phase,1)[1]],t))**2))),
                        'sample_count':len(t),
                        'source_hdf5_sha256':digest(path)
                    })
    missing = set(reported)-set(seen)
    if missing: anomalies.append(f'CSV keys not measured: {sorted(missing)}')
    fields=list(records[0])
    with (output/'independent-frequency-results.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=fields)
        wr.writeheader();wr.writerows(records)
    summary={
        'input_sha256':digest(root/'shield-evidence/dispersion-frequency-results.csv'),
        'measurements':len(records),
        'max_abs_diff_frequency_reread_vs_printed':max(x['abs_frequency_difference'] for x in records),
        'max_abs_relative_error_vs_exact_cold':max(abs(x['rel_error_vs_cold']) for x in records),
        'max_abs_relative_error_vs_empirical_scheme':max(abs(x['rel_error_vs_scheme_empirical']) for x in records),
        'max_phase_fit_rms_radians':max(x['phase_fit_residual_rms_radians'] for x in records),
        'amplitude_ratio_min':min(x['amp_end_over_start'] for x in records),
        'amplitude_ratio_max':max(x['amp_end_over_start'] for x in records),
        'anomalies':anomalies,
        'assessment': 'PASS_BASELINE_RAW_DIAGNOSTIC_REREAD' if len(records)==48 and not anomalies else 'FAIL_OR_INCOMPLETE',
        'note': 'This independently extracts frequencies from archived PHARE magnetic-field HDF5. It does not run PHARE or demonstrate convergence, external-code agreement, or Project Shield feasibility.'
    }
    (output/'independent-evaluation.json').write_text(json.dumps(summary, indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.3,5.4))
    markers={'aic':'o','whistler':'s'}
    for branch in ('aic','whistler'):
        for order in (1,2,3):
            ss=[r for r in records if r['branch']==branch and r['order']==order]
            ax.scatter([r['k_di_exact'] for r in ss],[100*r['rel_error_vs_cold'] for r in ss],marker=markers[branch],s=24,label=f'{branch} / p={order}',alpha=0.75)
    ax.axhline(0, color='black',lw=0.8)
    ax.set(xscale='log',xlabel=r'$k d_i$', ylabel='Relative frequency error vs cold theory (%)',title='Phase 3B.1B: independently reconstructed raw-HDF5 baseline (not a convergence plot)')
    ax.grid(True,which='both',alpha=0.25)
    ax.legend(fontsize=8,ncol=2)
    fig.tight_layout();fig.savefig(output/'baseline-cold-theory-errors.png',dpi=150);plt.close(fig)
    print(json.dumps(summary,indent=2))
    return 0 if summary['assessment'].startswith('PASS') else 1

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--archive-extracted',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args()
    raise SystemExit(evaluate(a.archive_extracted, a.output))
