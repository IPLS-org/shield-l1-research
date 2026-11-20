#!/usr/bin/env python3
"""Read actual PHARE HDF5 for ONE 128-v-64 1D spatial sensitivity screen.

Predeclared engineering gates. Not pure truncation error, statistical qualification,
multidimensional solver validation or Project Shield physical feasibility evidence.
Previously qualified raw HDF5 reader is imported unchanged, checked by SHA256.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import h5py
from p3b1d_perturbed_analyse import read_rows, screen, digest, extract_component, one_patch

PHARE_SHA = '7e5c72c1e5eec5452c8e3cfd000ce26c01e084ef'
BASELINE_ARCHIVE_SHA = '920659a199661afd765072877385cf87ba9237afc7b4105ee9b80bdf0caafdbf'
READER_SHA = 'edbbed49b5ff737caa7bc0426852263ac87a117733d7baac30d526e5bac72898'


def require_contract(protocol, reference, manifest, protocol_file, reference_file):
    """Fail closed before computing scientific comparisons."""
    if protocol['protocol_id'] != 'P3B1D-GRID-128-20261120-v0.6':
        raise ValueError('Wrong experiment identity')
    if protocol['status'] != 'REVIEW_DRAFT_NOT_AUTHORISATION_TO_DISPATCH':
        raise ValueError('Predeclared preregistration status changed')
    if digest(Path(__file__).with_name('p3b1d_perturbed_analyse.py')) != READER_SHA:
        raise ValueError('Previously qualified reader has been modified')
    if protocol['previous_reader_sha256'] != READER_SHA:
        raise ValueError('Reader version absent from protocol')
    if digest(reference_file) != protocol['baseline_reference_sha256']:
        raise ValueError('Frozen 64-cell baseline reference bytes changed')
    if (protocol['phare_commit'] != PHARE_SHA or reference['baseline_phare_commit'] != PHARE_SHA
        or manifest.get('phare_commit') != PHARE_SHA):
        raise ValueError('PHARE source mismatch')
    if (protocol['baseline_artifact_sha256'] != BASELINE_ARCHIVE_SHA
        or reference['baseline_raw_artifact_sha256'] != BASELINE_ARCHIVE_SHA):
        raise ValueError('Original evidence provenance mismatch')
    if reference['baseline_run_id'] != protocol['baseline_run_id'] or reference['baseline_run_id'] != 44805400036:
        raise ValueError('Wrong reference run')
    if reference['reference_status'] != 'INDEPENDENTLY_RECALCULATED_FROM_ACTUAL_PHARE_HDF5':
        raise ValueError('Reference is not from previous PHARE raw HDF5')
    c, b = protocol['case'], reference['baseline_case']
    permitted = {'name','cells','dx'}
    if set(c) != set(b) or {k for k in c if c[k] != b[k]} != permitted:
        raise ValueError('Unexpected physics or sampling configuration change')
    if set(protocol['changed_case_parameters_only']) != permitted:
        raise ValueError('Changed parameter declaration mismatch')
    if not (c['cells']==128 and b['cells']==64 and c['dx']==0.05 and b['dx']==0.1
        and c['domain_length']==b['domain_length']==6.4
        and c['particles_per_cell']==b['particles_per_cell']==400
        and c['time_step']==b['time_step']==0.0005 and c['solver_steps']==b['solver_steps']==2000
        and c['snapshot_every_steps']==b['snapshot_every_steps']==200
        and c['snapshot_times']==b['snapshot_times'] and c['final_time']==b['final_time']==1.0):
        raise ValueError('Not an unchanged-duration / controlled-grid pair')
    if c['seed'] != b['seed'] or c['seed'] != 314159:
        raise ValueError('Seed changed unexpectedly')
    if not(c['electron_temperature']==c['resistivity']==c['hyper_resistivity']==0):
        raise ValueError('Non-ideal physics switches disallow conservation gate')
    pc=protocol['grid_control']
    if not(pc['smallest_patch_size']==pc['largest_patch_size']==128 and pc['max_nbr_levels']==1):
        raise ValueError('Unreviewed refinement patch policy')
    if manifest.get('source') != 'PHARE_runtime' or manifest.get('protocol_id') != protocol['protocol_id']:
        raise ValueError('Not the intended real PHARE execution')
    if manifest.get('protocol_sha256') != digest(protocol_file) or manifest.get('case') != c:
        raise ValueError('Run input differs from frozen prospective case')
    if manifest.get('baseline_reference_sha256') != digest(reference_file) or manifest.get('baseline_run_id') != 44805400036:
        raise ValueError('Run uses different old reference')
    if manifest.get('patch_control') != pc:
        raise ValueError('Run patch control differs from frozen experiment')
    if len(reference['time_series']) != 11 or any(abs(z['t']-q)>1e-10 for z,q in zip(reference['time_series'],c['snapshot_times'])):
        raise ValueError('Previous reference missing aligned physical timestamps')
    return c


def grid_pair_screen(reference_rows, trial_rows, case, caps=None):
    """Metric-based paired comparison at matched physical wavenumber and times.

    Not cell-to-cell matching or exact initial particle identity. At fixed PPC,
    128 cells has double macro-particle count, thus a changed noise realisation.
    """
    if caps is None:
        caps = {'complex_mode_scaled':0.02,'kinetic_exchange_scaled':0.03,
                'transverse_exchange_scaled':0.03,'momentum_trajectory_scaled':0.02}
    if len(reference_rows) != 11 or len(trial_rows) != 11:
        raise ValueError('Missing paired observation time')
    delta=case['field_amplitude']; L=case['domain_length']; SE=L*delta*delta/2
    SP=L*case['density']*case['ion_mass']*delta/case['B0']
    previous=reference_rows[0]; new=trial_rows[0]
    out=[]
    for br,r in zip(reference_rows,trial_rows):
        if abs(float(r['t'])-float(br['t']))>1e-10:
            raise ValueError('Different observation times; no interpolation allowed')
        z0=complex(br['initial_wave_mode_real'],br['initial_wave_mode_imag'])
        z1=complex(r['initial_wave_mode_real'],r['initial_wave_mode_imag'])
        dw=abs(z1-z0)
        dk=abs((r['K']-new['K'])-(br['K']-previous['K']))
        dm=abs((r['UB_transverse']-new['UB_transverse'])-(br['UB_transverse']-previous['UB_transverse']))
        dp=float(np.linalg.norm((np.asarray(r['P'])-np.asarray(new['P']))-(np.asarray(br['P'])-np.asarray(previous['P']))))
        vals=[dw,dk,dm,dp]
        if not np.all(np.isfinite(vals)):
            raise ValueError('Nonfinite computed paired metric')
        out.append({'t':float(r['t']), 'wave_abs_complex_diff':dw, 'wave_diff_scaled':dw/delta,
                    'ion_kinetic_exchange_abs_diff':dk,'ion_kinetic_exchange_scaled':dk/SE,
                    'magnetic_transverse_exchange_abs_diff':dm,'magnetic_exchange_scaled':dm/SE,
                    'momentum_trajectory_abs_diff':dp,'momentum_trajectory_scaled':dp/SP,
                    'coarse_complex_mode_real':br['initial_wave_mode_real'],'coarse_complex_mode_imag':br['initial_wave_mode_imag'],
                    'fine_complex_mode_real':r['initial_wave_mode_real'],'fine_complex_mode_imag':r['initial_wave_mode_imag'],
                    'coarse_ion_K':br['K'],'fine_ion_K':r['K'],
                    'coarse_UBtrans':br['UB_transverse'],'fine_UBtrans':r['UB_transverse']})
    keys=[('S1_complex_wave','wave_diff_scaled','complex_mode_scaled'),
          ('S2_kinetic_exchange','ion_kinetic_exchange_scaled','kinetic_exchange_scaled'),
          ('S2_magnetic_exchange','magnetic_exchange_scaled','transverse_exchange_scaled'),
          ('S3_momentum_trajectory','momentum_trajectory_scaled','momentum_trajectory_scaled')]
    gates={name:{'pass':max(o[field] for o in out)<=caps[cap],
                 'observed_max_scaled':max(o[field] for o in out),
                 'registered_max_scaled':caps[cap],
                 'reference_scale':delta if name=='S1_complex_wave' else SP if name=='S3_momentum_trajectory' else SE}
           for name,field,cap in keys}
    initial={'coarse_particle_count':previous['particle_count'],'fine_particle_count':new['particle_count'],
             'same_macro_particle_population':False,'same_random_seed':True,'same_particles':False,
             'coarse_initial_K':previous['K'],'fine_initial_K':new['K'],'initial_K_difference':new['K']-previous['K'],
             'coarse_initial_P':previous['P'],'fine_initial_P':new['P'],
             'initial_momentum_difference_norm':float(np.linalg.norm(np.asarray(new['P'])-np.asarray(previous['P']))),
             'coarse_initial_wave':previous['wave_amplitude'],'fine_initial_wave':new['wave_amplitude'],
             'initial_wave_complex_difference':out[0]['wave_abs_complex_diff'],
             'stochastic_caveat':'Particle count doubles at fixed 400 ppc; same seed does not preserve individual particle positions or velocities.'}
    return {'gates':gates,'all_pass':all(g['pass'] for g in gates.values()),'paired_rows':out,
            'initial_sampling_diagnostics':initial,'caps':caps,'reference_scales':{'delta':delta,'SE':SE,'SP':SP}}


def initial_density_rms(diagdir,case):
    """Observation-only initial deposited charge fluctuations as PIC noise indicator."""
    with h5py.File(diagdir/'ions_charge_density.h5','r') as h:
        patch=one_patch(h,0.,case,'ions_charge_density.h5')
        rho=extract_component(patch['charge_density'],'primal',case['cells'])
        # A primal periodic endpoint duplicates the first physical location.
        owned=rho[:-1]/case['ion_charge']
        return float(np.std(owned-case['density'])),float(np.mean(owned))


def analyse(proto_path, ref_path, diagdir, outdir):
    outdir.mkdir(parents=True,exist_ok=True)
    report={'status':'INCONCLUSIVE_DATA_CONTRACT',
            'scientific_scope':'ONE 1D grid sensitivity screen at fixed PPC, not an asymptotic convergence order or Project Shield validation',
            'new_grid_128_particles_51200_vs_coarse_64_particles_25600':True}
    try:
        proto=json.loads(proto_path.read_text()); ref=json.loads(ref_path.read_text())
        manifest=json.loads((diagdir/'p3b1d_grid128_run.json').read_text())
        case=require_contract(proto,ref,manifest,proto_path,ref_path)
        rows,sizes=read_rows(diagdir,proto)
        old=screen(rows,proto)
        pair=grid_pair_screen(ref['time_series'],rows,case,proto['numeric_caps'])
        initial_rho_sd,initial_rho_mean=initial_density_rms(diagdir,case)
        prev_gates=old['gates']; all_original=old['status']=='BOUNDED_SINGLE_CASE_SCREEN_PASS'
        report.update({'case':case,'reference_run':44805400036,'phare_commit':PHARE_SHA,
                       'research_commit':manifest.get('research_commit'),
                       'protocol_sha256':digest(proto_path),'baseline_reference_sha256':digest(ref_path),
                       'previous_raw_artifact_sha256':BASELINE_ARCHIVE_SHA,
                       'raw_hdf5_sha256':{k:digest(diagdir/k) for k in sizes},
                       'total_hdf5_bytes':sum(sizes.values()),
                       'new_case_gates':prev_gates,'original_screen_all_pass':all_original,
                       'paired_comparison':pair,'fine_initial_deposited_density_cell_rms':initial_rho_sd,
                       'fine_initial_deposited_density_cell_mean':initial_rho_mean,
                       'final_single_case_observations':old['rows'][-1],
                       'reasons_this_is_not_grid_convergence':['two grid resolutions only',
                             'fixed PPC doubles the macro-particles and changes sampling noise',
                             'no third grid or independent seeds', 'no exact continuum truth']})
        report['status']='BOUNDED_GRID_SENSITIVITY_SCREEN_PASS' if (all_original and pair['all_pass']) else 'BOUNDED_GRID_SENSITIVITY_SCREEN_FAIL'
        with (outdir/'grid_pair_timeseries.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(pair['paired_rows'][0]));writer.writeheader();writer.writerows(pair['paired_rows'])
    except Exception as exc:
        report['error']=type(exc).__name__+': '+str(exc)
    (outdir/'grid_sensitivity_assessment.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({x:report[x] for x in ('status','error') if x in report},indent=2))
    return 0 if report['status']=='BOUNDED_GRID_SENSITIVITY_SCREEN_PASS' else 3


def main():
    p=argparse.ArgumentParser();p.add_argument('--protocol',type=Path,required=True);p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--diagdir',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args();return analyse(a.protocol,a.reference,a.diagdir,a.outdir)

if __name__=='__main__':raise SystemExit(main())
