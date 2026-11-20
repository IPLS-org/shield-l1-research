#!/usr/bin/env python3
"""Prospective 1D PHARE dt-half sensitivity screen, not convergence or Shield validation.

Import the already qualified v0.4 HDF5 reader for identical raw-data conventions,
then add independent paired criteria frozen before generating a new PHARE run.
The post-run scientific audit must independently re-read raw HDF5 again.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import h5py
import numpy as np
from p3b1d_perturbed_analyse import read_rows, screen, digest

PHARE_SHA = '7e5c72c1e5eec5452c8e3cfd000ce26c01e084ef'
PREVIOUS_READER_SHA256 = 'edbbed49b5ff737caa7bc0426852263ac87a117733d7baac30d526e5bac72898'
BASELINE_ARTIFACT_SHA256='4eeb3a5d1c554c4b407dee120c94fb9c94c098ed207d73c8dd723097df8c8f07'
FIELDS0={
    'EM_B.h5':('EM_B_x','EM_B_y','EM_B_z'),
    'ions_pop_protons_domain.h5':('weight','iCell','delta','v'),
}


def validate_contract(protocol, reference, manifest, protocol_path, reference_path):
    if protocol['protocol_id']!='P3B1D-DT-HALF-20261120-v0.5':
        raise ValueError('Wrong prospective study ID')
    if protocol['phare_commit']!=PHARE_SHA or reference['baseline_phare_commit']!=PHARE_SHA:
        raise ValueError('Pinned PHARE source has changed')
    if reference['baseline_raw_artifact_sha256']!=BASELINE_ARTIFACT_SHA256 or protocol['baseline_artifact_sha256']!=BASELINE_ARTIFACT_SHA256:
        raise ValueError('Baseline provenance mismatch')
    if digest(reference_path)!=protocol['baseline_reference_sha256']:
        raise ValueError('Reference bytes changed after prospectively frozen protocol')
    if digest(Path(__file__).with_name('p3b1d_perturbed_analyse.py'))!=PREVIOUS_READER_SHA256:
        raise ValueError('Prior validated HDF5 reader changed; abort prospective comparison')
    if manifest.get('source')!='PHARE_runtime' or manifest.get('phare_commit')!=PHARE_SHA:
        raise ValueError('Not a provenance-verified PHARE run')
    if manifest.get('protocol_id')!=protocol['protocol_id'] or manifest.get('protocol_sha256')!=digest(protocol_path):
        raise ValueError('Research run does not match frozen sensitivity protocol')
    if manifest.get('baseline_reference_sha256')!=protocol['baseline_reference_sha256'] or manifest.get('baseline_run_id')!=44796208433:
        raise ValueError('Run used wrong baseline reference')
    c=protocol['case']; b=reference['baseline_case']
    if manifest.get('case')!=c:
        raise ValueError('Run-time model does not equal frozen new case')
    allowed={'name','time_step','solver_steps','snapshot_every_steps'}
    if set(c)!=set(b) or {k for k in c if c[k]!=b[k]}!=allowed:
        raise ValueError('More than the allowed controlled time-step parameters changed')
    if not(c['time_step']==0.0005 and b['time_step']==0.001 and c['solver_steps']==2000 and b['solver_steps']==1000
           and c['snapshot_every_steps']==200 and b['snapshot_every_steps']==100
           and c['final_time']==b['final_time']==1 and c['snapshot_times']==b['snapshot_times']):
        raise ValueError('Paired time-step/duration/cadence invalid')
    if reference['baseline_protocol_sha256']!='4a9e52f7165fc23540dd5f743efe95dd38786c999824ddf829ffe48748291c53':
        raise ValueError('Wrong original v0.4 baseline protocol')
    if len(reference['time_series'])!=11:
        raise ValueError('Baseline physical observation count not 11')
    for i, (br,t) in enumerate(zip(reference['time_series'],c['snapshot_times'])):
        if not math.isclose(br['t'],t,abs_tol=1e-11):
            raise ValueError(f'Baseline time alignment invalid at index {i}')
    return c


def initial_digest_equality(diagdir, reference):
    """Compare logical initial datasets exactly, independent of HDF5 layout/metadata."""
    checks={}
    for filename,keys in FIELDS0.items():
        with h5py.File(diagdir/filename,'r') as h:
            p=next(iter(h['t']['0.0000000000']['pl0'].values()))
            for key in keys:
                ds=p[key]
                if filename=='EM_B.h5':
                    ghosts=int(ds.attrs['ghosts']); physical_nodes=64+(1 if key=='EM_B_x' else 0)
                    observed=np.ascontiguousarray(ds[ghosts:ghosts+physical_nodes])
                else:
                    observed=np.ascontiguousarray(ds[()])
                path=filename+':'+key
                expected=reference['initial_logical_dataset_digests'][path]
                actualhash=hashlib.sha256(observed.tobytes(order='C')).hexdigest()
                checks[path]={'match':(actualhash==expected['sha256'] and list(observed.shape)==expected['shape'] and
                                       str(observed.dtype)==expected['dtype']),
                              'expected_sha256':expected['sha256'],'observed_sha256':actualhash,
                              'dtype':str(observed.dtype),'shape':list(observed.shape)}
    return {'pass':all(v['match'] for v in checks.values()),'datasets':checks}


def pairwise_screen(reference_rows, fine_rows, case, caps=None):
    """Apply absolute perturbation-normalised engineering screens, not a convergence theorem."""
    if caps is None:caps={'wave':0.01,'energy_exchange':0.02,'momentum':0.01}
    if len(fine_rows)!=len(reference_rows):
        raise ValueError('Paired sample length mismatch')
    delta=case['field_amplitude']; SE=case['domain_length']*delta**2/2
    SP=case['domain_length']*case['density']*case['ion_mass']*delta/case['B0']
    wb=[]; ke=[]; ub=[]; mom=[]; output=[]
    p0=np.array(fine_rows[0]['P'],dtype=float);pb0=np.array(reference_rows[0]['P'],dtype=float)
    k0=fine_rows[0]['K']; kb0=reference_rows[0]['K']
    u0=fine_rows[0]['UB_transverse'];ub0=reference_rows[0]['UB_transverse']
    for a,b in zip(fine_rows,reference_rows):
        if not math.isclose(a['t'],b['t'],rel_tol=0,abs_tol=1e-10):
            raise ValueError('Non-identical observation time; no interpolation permitted')
        za=complex(a['initial_wave_mode_real'],a['initial_wave_mode_imag'])
        zb=complex(b['initial_wave_mode_real'],b['initial_wave_mode_imag'])
        dv=abs(za-zb); dK=abs((a['K']-k0)-(b['K']-kb0)); dUB=abs((a['UB_transverse']-u0)-(b['UB_transverse']-ub0))
        dp=float(np.linalg.norm((np.array(a['P'])-p0)-(np.array(b['P'])-pb0)))
        if not all(math.isfinite(x) for x in (dv,dK,dUB,dp)):
            raise ValueError('Non-finite paired metrics')
        wb.append(dv);ke.append(dK);ub.append(dUB);mom.append(dp)
        output.append({'time':a['t'],'complex_mode_abs_diff':dv,'kinetic_exchange_abs_diff':dK,
                       'transverse_field_exchange_abs_diff':dUB,'momentum_trajectory_abs_diff':dp,
                       'fine_wave_amplitude':a['wave_amplitude'],'base_wave_amplitude':b['wave_amplitude'],
                       'fine_K':a['K'],'base_K':b['K'],'fine_UBtrans':a['UB_transverse'],'base_UBtrans':b['UB_transverse']})
    g={
        'S1_first_harmonic':{'pass':max(wb)/delta<=caps['wave'],'max_abs_diff':max(wb),'scaled':max(wb)/delta,'relative_perturbation_limit':caps['wave']},
        'S2_kinetic_exchange':{'pass':max(ke)/SE<=caps['energy_exchange'],'max_abs_diff':max(ke),'scaled':max(ke)/SE,'relative_perturbation_energy_limit':caps['energy_exchange']},
        'S2_magnetic_exchange':{'pass':max(ub)/SE<=caps['energy_exchange'],'max_abs_diff':max(ub),'scaled':max(ub)/SE,'relative_perturbation_energy_limit':caps['energy_exchange']},
        'S3_momentum_trajectory':{'pass':max(mom)/SP<=caps['momentum'],'max_abs_diff':max(mom),'scaled':max(mom)/SP,'relative_perturbation_momentum_limit':caps['momentum']},
    }
    return {'gates':g,'all_pass':all(x['pass'] for x in g.values()),'rows':output,'reference_scales':{'delta':delta,'SE':SE,'SP':SP},'limits':caps}


def run_analysis(proto_path, ref_path, diagdir, outdir):
    outdir.mkdir(parents=True,exist_ok=True)
    report={'status':'INCONCLUSIVE_DATA_CONTRACT','strict_scope':'ONE time-step sensitivity comparison, not convergence order or Shield validation',
            'not_convergence_proof':True,'not_physical_validation':True}
    try:
        proto=json.loads(proto_path.read_text());ref=json.loads(ref_path.read_text());m=json.loads((diagdir/'p3b1d_dt_half_run.json').read_text())
        case=validate_contract(proto,ref,m,proto_path,ref_path)
        rows,sizes=read_rows(diagdir,proto)
        original_gates=screen(rows,proto)
        initial=initial_digest_equality(diagdir,ref)
        pair=pairwise_screen(ref['time_series'],rows,case)
        report.update({'case':case,'comparison_to_run':ref['baseline_run_id'],'baseline_raw_zip_sha256':BASELINE_ARTIFACT_SHA256,
                       'baseline_reference_sha256':digest(ref_path),'protocol_sha256':digest(proto_path),
                       'phare_commit':PHARE_SHA,'research_commit':m.get('research_commit'),'HDF5_sha256':{x:digest(diagdir/x) for x in sizes},
                       'total_hdf5_bytes':sum(sizes.values()),'S0_identical_initial_arrays':initial,
                       'fine_case_gates':original_gates['gates'],'fine_case_raw_rows':original_gates['rows'],
                       'pair_sensitivity':pair,'S4_original_single_case_gates_pass':original_gates['status']=='BOUNDED_SINGLE_CASE_SCREEN_PASS',
                       'trend_warning':'Two dt values with one grid and seed cannot establish convergence order or numerical method consistency in general.'})
        if not initial['pass']:
            report['status']='INCONCLUSIVE_INITIAL_STATE'
        else:
            report['status']='BOUNDED_DT_SENSITIVITY_SCREEN_PASS' if (pair['all_pass'] and report['S4_original_single_case_gates_pass']) else 'BOUNDED_DT_SENSITIVITY_SCREEN_FAIL'
        with (outdir/'paired_timeseries.csv').open('w',newline='') as f:
            fields=['time','complex_mode_abs_diff','kinetic_exchange_abs_diff','transverse_field_exchange_abs_diff','momentum_trajectory_abs_diff',
                    'fine_wave_amplitude','base_wave_amplitude','fine_K','base_K','fine_UBtrans','base_UBtrans']
            dw=csv.DictWriter(f,fieldnames=fields);dw.writeheader();dw.writerows(pair['rows'])
    except Exception as e:
        report['error']=type(e).__name__+': '+str(e)
    (outdir/'sensitivity_assessment.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({key:report[key] for key in ['status','error'] if key in report},indent=2))
    return 0 if report['status']=='BOUNDED_DT_SENSITIVITY_SCREEN_PASS' else 3


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--protocol',type=Path,required=True);ap.add_argument('--reference',type=Path,required=True)
    ap.add_argument('--diagdir',type=Path,required=True);ap.add_argument('--outdir',type=Path,required=True)
    a=ap.parse_args();return run_analysis(a.protocol,a.reference,a.diagdir,a.outdir)

if __name__=='__main__':raise SystemExit(main())
