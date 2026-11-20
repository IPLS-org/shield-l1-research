#!/usr/bin/env python3
"""Prospective 1D PHARE seed-vs-seed screen, not statistical variance or Shield validation.

Read only actual raw HDF5 using the previously qualified unchanged owned-node reader.
Do not regard workflow exit 0 or synthetic fixtures as independent physical proof.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import h5py
from p3b1d_perturbed_analyse import read_rows,screen,digest,one_patch,extract_component,EXPECTED

PHARE_SHA='7e5c72c1e5eec5452c8e3cfd000ce26c01e084ef'
READER_SHA='edbbed49b5ff737caa7bc0426852263ac87a117733d7baac30d526e5bac72898'
PREVIOUS_ARTIFACT_SHA='a66631911992cbaebc92e4d056901ea30bb5330afb9ff154df4c57640f62c9e7'
BASELINE=44819929641
FIELD_CENTRES={'EM_B_x':'primal','EM_B_y':'dual','EM_B_z':'dual'}


def sha_array(a):
    a=np.ascontiguousarray(a)
    return hashlib.sha256(str(a.dtype).encode()+str(a.shape).encode()+a.tobytes()).hexdigest()


def require_contract(p,r,manifest,proto_path,ref_path):
    if p.get('protocol_id')!='P3B1D-SEED-128-20261120-v0.7':raise ValueError('Wrong experiment identity')
    if p.get('status')!='REVIEW_DRAFT_NOT_AUTHORISATION_TO_DISPATCH':raise ValueError('Protocol status has changed')
    if digest(Path(__file__).with_name('p3b1d_perturbed_analyse.py'))!=READER_SHA or p['previous_reader_sha256']!=READER_SHA:
        raise ValueError('Old qualified reader not unchanged')
    if digest(ref_path)!=p['baseline_reference_sha256']:raise ValueError('Frozen reference SHA mismatch')
    if p['baseline_artifact_sha256']!=PREVIOUS_ARTIFACT_SHA or r['baseline_raw_artifact_sha256']!=PREVIOUS_ARTIFACT_SHA:
        raise ValueError('Prior raw artifact provenance changed')
    if p['phare_commit']!=PHARE_SHA or r['baseline_phare_commit']!=PHARE_SHA or manifest.get('phare_commit')!=PHARE_SHA:
        raise ValueError('PHARE version not qualified')
    if r['reference_status']!='RECALCULATED_FROM_REAL_PREVIOUS_PHARE_HDF5':raise ValueError('Reference origin not real PHARE HDF5')
    if r['baseline_run_id']!=p['baseline_run_id'] or r['baseline_run_id']!=BASELINE:
        raise ValueError('Wrong original 128-cell seed run')
    c,b=p['case'],r['baseline_case']
    if set(c)!=set(b) or {k for k in c if c[k]!=b[k]}!={'seed','name'}:
        raise ValueError('Only seed and case-name may change')
    if set(p['changed_case_parameters_only'])!={'seed','name'}:raise ValueError('Improper permitted difference declaration')
    if not (c['seed']==271828 and b['seed']==314159 and c['name']=='P1_grid128_second_seed_271828'):
        raise ValueError('Seed pairing incorrect')
    if not (c['cells']==b['cells']==128 and c['dx']==b['dx']==.05 and c['domain_length']==b['domain_length']==6.4
        and c['particles_per_cell']==b['particles_per_cell']==400 and c['time_step']==b['time_step']==.0005
        and c['solver_steps']==b['solver_steps']==2000 and c['snapshot_every_steps']==b['snapshot_every_steps']==200
        and c['snapshot_times']==b['snapshot_times'] and c['final_time']==b['final_time']==1.0):
        raise ValueError('Not equal numerical grid or physical duration')
    if not(c['electron_temperature']==c['resistivity']==c['hyper_resistivity']==0):raise ValueError('Ideal physical settings changed')
    patch=p['grid_control']
    if not(patch['smallest_patch_size']==patch['largest_patch_size']==128 and patch['max_nbr_levels']==1):
        raise ValueError('Grid patch source declaration mismatch')
    if not(manifest.get('source')=='PHARE_runtime' and manifest.get('protocol_id')==p['protocol_id']
        and manifest.get('case')==c and manifest.get('protocol_sha256')==digest(proto_path)
        and manifest.get('baseline_reference_sha256')==digest(ref_path)
        and manifest.get('baseline_run_id')==BASELINE and manifest.get('patch_control')==patch):
        raise ValueError('PHARE runtime provenance contract failed')
    if len(r['time_series'])!=11 or any(abs(a['t']-t)>1e-11 for a,t in zip(r['time_series'],c['snapshot_times'])):
        raise ValueError('Reference time series missing/misaligned')
    assert p['baseline_grid_pair_observed_differences']==r['prior_grid_difference_same_seed_64to128']
    return c


def initial_arrays(diagdir,c):
    """Exact baseline/new owned magnetic field hashes and initial particles; ghost fields excluded."""
    b={}
    with h5py.File(diagdir/'EM_B.h5','r') as f:
        patch=one_patch(f,0.,c,'EM_B.h5')
        for field,centre in FIELD_CENTRES.items():
            arr=extract_component(patch[field],centre,c['cells'])
            b[field]={'sha256':sha_array(arr),'shape':list(arr.shape),'dtype':str(arr.dtype)}
    particles={}
    with h5py.File(diagdir/'ions_pop_protons_domain.h5','r') as f:
        patch=one_patch(f,0.,c,'ions_pop_protons_domain.h5')
        for field in ('weight','iCell','delta','v'):
            a=np.asarray(patch[field][...]);particles[field]={'sha256':sha_array(a),'shape':list(a.shape),'dtype':str(a.dtype)}
    with h5py.File(diagdir/'ions_charge_density.h5','r') as f:
        patch=one_patch(f,0.,c,'ions_charge_density.h5')
        rho=extract_component(patch['charge_density'],'primal',c['cells'])[:-1]/c['ion_charge']
    return b,particles,{'density_rms':float(np.std(rho-c['density'])),'density_mean':float(np.mean(rho))}


def paired_screen(ref_rows,now_rows,case,caps,grid_comparison):
    if len(ref_rows)!=11 or len(now_rows)!=11:raise ValueError('Missing one of eleven observation times')
    scale_wave=case['field_amplitude']; scale_energy=.5*case['domain_length']*scale_wave**2
    scale_momentum=case['domain_length']*case['density']*case['ion_mass']*scale_wave/case['B0']
    a=ref_rows[0];b=now_rows[0]
    comparison=[]
    for x,y in zip(ref_rows,now_rows):
        if abs(float(x['t'])-float(y['t']))>1e-11:raise ValueError('Time misalignment; no interpolation permitted')
        wave=abs(complex(y['initial_wave_mode_real'],y['initial_wave_mode_imag'])-
                 complex(x['initial_wave_mode_real'],x['initial_wave_mode_imag']))
        ki=abs((y['K']-b['K'])-(x['K']-a['K']))
        ub=abs((y['UB_transverse']-b['UB_transverse'])-(x['UB_transverse']-a['UB_transverse']))
        mom=float(np.linalg.norm((np.asarray(y['P'])-np.asarray(b['P']))-
                                 (np.asarray(x['P'])-np.asarray(a['P']))))
        if not all(math.isfinite(z) for z in (wave,ki,ub,mom)):
            raise ValueError('Nonfinite paired values')
        comparison.append({'t':float(y['t']),'wave_abs_complex_diff':wave,'wave_scaled':wave/scale_wave,
             'ion_exchange_abs_diff':ki,'ion_exchange_scaled':ki/scale_energy,
             'magnetic_exchange_abs_diff':ub,'magnetic_exchange_scaled':ub/scale_energy,
             'momentum_evolution_abs_diff':mom,'momentum_evolution_scaled':mom/scale_momentum,
             'seed314159_wave_real':x['initial_wave_mode_real'],'seed314159_wave_imag':x['initial_wave_mode_imag'],
             'seed271828_wave_real':y['initial_wave_mode_real'],'seed271828_wave_imag':y['initial_wave_mode_imag'],
             'seed314159_K':x['K'],'seed271828_K':y['K'],
             'seed314159_UB_transverse':x['UB_transverse'],'seed271828_UB_transverse':y['UB_transverse']})
    keys=[('S1_complex_mode','wave_scaled','complex_mode_scaled','wave_abs_complex_diff'),
          ('S2_ion_exchange','ion_exchange_scaled','kinetic_exchange_scaled','ion_exchange_abs_diff'),
          ('S2_magnetic_exchange','magnetic_exchange_scaled','transverse_exchange_scaled','magnetic_exchange_abs_diff'),
          ('S3_momentum_trajectory','momentum_evolution_scaled','momentum_trajectory_scaled','momentum_evolution_abs_diff')]
    gates={name:{'pass':max(row[scaled] for row in comparison)<=caps[cap],
                 'observed_max_scaled':max(row[scaled] for row in comparison),
                 'registered_max_scaled':caps[cap],
                 'observed_max_absolute':max(row[ab] for row in comparison)}
          for name,scaled,cap,ab in keys}
    ratios={}
    link={'S1_complex_mode':'wave_abs_complex_diff','S2_ion_exchange':'ion_exchange_abs_diff',
          'S2_magnetic_exchange':'magnetic_exchange_abs_diff','S3_momentum_trajectory':'momentum_evolution_abs_diff'}
    for name,field in link.items():
        previous=grid_comparison[field]
        ratios[name]={'new_seed_pair_observed_difference':gates[name]['observed_max_absolute'],
            'previous_grid_pair_observed_difference':previous,
            'new_seed_over_prior_grid_observed_ratio':gates[name]['observed_max_absolute']/previous if previous>0 else None}
    return {'gates':gates,'all_pass':all(z['pass'] for z in gates.values()),'paired_rows':comparison,
        'observed_seed_over_grid_ratios_DESCRIPTIVE_ONLY':ratios,
        'statistical_limit':'Only 2 seeds and 1 pair at one grid: NOT a standard deviation, confidence interval or partition of grid error.'}


def analyse(proto_path,ref_path,diagdir,outdir):
    outdir.mkdir(parents=True,exist_ok=True)
    report={'status':'INCONCLUSIVE_DATA_CONTRACT',
        'scope':'Exactly one 128-cell seed perturbation. A screen and observed difference, NOT a variance estimate, convergence order, or Project Shield validation.'}
    try:
        p=json.loads(proto_path.read_text());r=json.loads(ref_path.read_text())
        manifest=json.loads((diagdir/'p3b1d_seed128_run.json').read_text())
        c=require_contract(p,r,manifest,proto_path,ref_path)
        rows,sizes=read_rows(diagdir,p)
        physical=screen(rows,p)
        B,particles,noise=initial_arrays(diagdir,c)
        if B != r['initial_owned_B_sha256']:
            report['status']='INCONCLUSIVE_INITIAL_FIELD'
            raise RuntimeError('Initial physical B nodes differ; cannot attribute change exclusively to particle seed')
        if not any(particles[k]['sha256']!=r['initial_particles_sha256'][k]['sha256'] for k in ('delta','v')):
            report['status']='INVALID_CASE'
            raise RuntimeError('Particle seed has not altered any initial positions or velocities')
        # Changing a random seed can alter particle ordering and permissible weights;
        # physical number and weighting are verified by the existing HDF5/G2 gates.
        pair=paired_screen(r['time_series'],rows,c,p['numeric_caps'],r['prior_grid_difference_same_seed_64to128'])
        # A separate hard limit prevents classifying absurdly shifted initial scalar moments as legitimate.
        initial={'baseline_K':r['time_series'][0]['K'],'new_K':rows[0]['K'],
                 'initial_K_abs_diff':abs(r['time_series'][0]['K']-rows[0]['K']),
                 'baseline_P':r['time_series'][0]['P'],'new_P':rows[0]['P'],
                 'initial_P_difference_norm':float(np.linalg.norm(np.asarray(r['time_series'][0]['P'])-np.asarray(rows[0]['P']))),
                 'baseline_density_rms':r['initial_density']['density_rms'],'new_density_rms':noise['density_rms'],
                 'initial_B_exact_hash_match':True,'initial_particle_positions_or_velocities_changed':True,
                 'initial_particle_count_unchanged':True}
        report.update({'experiment_id':p['protocol_id'],'case':c,'reference_run':BASELINE,
           'phare_commit':PHARE_SHA,'research_commit':manifest.get('research_commit'),
           'protocol_sha256':digest(proto_path),'baseline_reference_sha256':digest(ref_path),
           'baseline_artifact_sha256':PREVIOUS_ARTIFACT_SHA,
           'raw_hdf5_sha256':{name:digest(diagdir/name) for name in sizes},
           'total_hdf5_bytes':sum(sizes.values()),'single_case_gates':physical['gates'],
           'single_case_pass':physical['status']=='BOUNDED_SINGLE_CASE_SCREEN_PASS',
           'paired_seed_comparison':pair,'initial_sampling_observations':initial,
           'final_case_observations':physical['rows'][-1],
           'post_run_decision':'Hold separate GO/NO-GO assessment for model progression; do not launch further seed sweeps automatically.',
           'limitations':['two seeds at one grid do not estimate uncertainty variance',
             'unchanged deterministic initial magnetic field does not imply identical particle moments',
             'seed differences cannot be attributed to mesh or physical validation']})
        report['status']='BOUNDED_SEED_SENSITIVITY_SCREEN_PASS' if (report['single_case_pass'] and pair['all_pass']) else 'BOUNDED_SEED_SENSITIVITY_SCREEN_FAIL'
        with (outdir/'seed_pair_timeseries.csv').open('w',newline='') as file:
            writer=csv.DictWriter(file,fieldnames=list(pair['paired_rows'][0]));writer.writeheader();writer.writerows(pair['paired_rows'])
    except Exception as exc:
        if report['status']=='INCONCLUSIVE_DATA_CONTRACT' or report['status'] in ('INCONCLUSIVE_INITIAL_FIELD','INVALID_CASE'):
            report['error']=type(exc).__name__+': '+str(exc)
        else:
            report['error']=type(exc).__name__+': '+str(exc)
            report['status']='INCONCLUSIVE_DATA_CONTRACT'
    (outdir/'seed_sensitivity_assessment.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','error') if k in report},indent=2))
    return 0 if report['status']=='BOUNDED_SEED_SENSITIVITY_SCREEN_PASS' else 3


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--protocol',type=Path,required=True);p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--diagdir',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args();return analyse(a.protocol,a.reference,a.diagdir,a.outdir)

if __name__=='__main__':raise SystemExit(main())
