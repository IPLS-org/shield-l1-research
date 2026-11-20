#!/usr/bin/env python3
"""Independent raw-HDF5 one-case conservation screening; never Shield efficacy proof.

Strictly observation-only. Do not import the PHARE runtime or trust its self-reported gate results.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import h5py
import numpy as np

EXPECTED = {
    'EM_B.h5': ('EM_B_x','EM_B_y','EM_B_z'),
    'EM_E.h5': ('EM_E_x','EM_E_y','EM_E_z'),
    'ions_charge_density.h5': ('charge_density',),
    'ions_charge_flux.h5': ('charge_flux_x','charge_flux_y','charge_flux_z'),
    'ions_pop_protons_domain.h5': ('weight','iCell','delta','v'),
}
CENTERING = {
    'EM_B_x': 'primal', 'EM_B_y':'dual', 'EM_B_z':'dual',
    'EM_E_x': 'dual', 'EM_E_y':'primal', 'EM_E_z':'primal',
    'charge_density':'primal',
    'charge_flux_x':'primal','charge_flux_y':'primal','charge_flux_z':'primal',
}
PHARE_SHA = '7e5c72c1e5eec5452c8e3cfd000ce26c01e084ef'


def digest(path):
    hashobj=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):hashobj.update(b)
    return hashobj.hexdigest()


def validate_config(protocol, manifest):
    c=protocol['case']
    assert protocol['protocol_id']=='P3B1D-PERTURBED-ONE-20261120-v0.4', 'unexpected protocol ID'
    assert protocol['phare_commit']==PHARE_SHA, 'PHARE source not approved'
    assert manifest['source']=='PHARE_runtime', 'not a PHARE runtime manifest'
    assert manifest['phare_commit']==PHARE_SHA, 'PHARE source mismatch'
    assert manifest['protocol_id']==protocol['protocol_id'], 'protocol ID mismatch'
    assert manifest['case']==c, 'recorded case does not equal frozen protocol'
    assert c['cells']==64 and c['particles_per_cell']==400 and c['seed']==314159
    assert c['electron_temperature']==c['resistivity']==c['hyper_resistivity']==0
    assert c['time_step']*c['solver_steps']==c['final_time']==1.0
    assert c['snapshot_every_steps']==100 and c['solver_steps']==1000
    assert c['field_amplitude']==0.01 and c['domain_length']==6.4
    assert c['interpolation_order']==2
    return c


def extract_component(ds,centre,cells):
    if ds.dtype!=np.dtype('float64'):
        raise ValueError(f'{ds.name} is {ds.dtype}, expected float64')
    if 'ghosts' not in ds.attrs:
        raise ValueError(f'{ds.name} ghost attr missing')
    ghosts=int(np.asarray(ds.attrs['ghosts']).item())
    a=np.asarray(ds[...])
    n=cells+(1 if centre=='primal' else 0)
    if a.ndim!=1 or a.size!=n+2*ghosts:
        raise ValueError(f'{ds.name} shape {a.shape} incompatible with {n} physical nodes and {ghosts} ghosts')
    part=a[ghosts:ghosts+n]
    if not np.all(np.isfinite(part)):
        raise ValueError(f'{ds.name} nonfinite')
    return part


def verify_hdf5_root(h,c,name):
    for required in ('dimension','cell_width','origin','periodicities','domain_box','interpOrder','layoutType'):
        if required not in h.attrs:raise ValueError(f'{name} missing {required} HDF5 metadata')
    if int(np.asarray(h.attrs['dimension']).item())!=1:raise ValueError('dimension not 1')
    if int(np.asarray(h.attrs['interpOrder']).item())!=c['interpolation_order']:raise ValueError('interpolation order mismatch')
    if str(h.attrs['layoutType'])!='yee':raise ValueError('layoutType must be yee')
    if np.asarray(h.attrs['cell_width']).reshape(-1).size!=1 or abs(float(h.attrs['cell_width'][0])-c['dx'])>1e-12:raise ValueError('dx mismatch')
    if int(np.asarray(h.attrs['domain_box']).reshape(-1)[0])!=c['cells']-1:raise ValueError('domain_box mismatch')
    if not bool(np.asarray(h.attrs['periodicities']).reshape(-1)[0]):raise ValueError('nonperiodic')
    if not np.isclose(float(np.asarray(h.attrs['origin']).reshape(-1)[0]),0.0,atol=1e-12):raise ValueError('origin mismatch')
    if 't' not in h:raise ValueError('missing t group')
    actual=sorted(float(x) for x in h['t'].keys())
    times=c['snapshot_times']
    if len(actual)!=len(times) or not np.allclose(actual,times,atol=1e-8,rtol=0):
        raise ValueError(f'{name} time groups do not match prescribed {times}: got {actual}')
    return actual


def one_patch(h,t,c,name):
    group=h['t'][f'{t:.10f}']
    if list(group.keys())!=['pl0']:
        raise ValueError(f'{name}: expected one level 0 only')
    patches=group['pl0']
    if len(patches)!=1:
        raise ValueError(f'{name}: this bounded first case requires exactly one physical patch; got {len(patches)}')
    p=next(iter(patches.values()))
    def at(key):
        return int(np.asarray(p.attrs[key]).reshape(-1)[0])
    if (at('lower'),at('upper'),at('nbrCells'))!=(0,c['cells']-1,c['cells']):
        raise ValueError(f'{name}: incomplete or duplicate physical patch')
    return p


def read_rows(diagdir, protocol):
    c=protocol['case']; n=c['cells']; dx=c['dx']; expectedmacro=n*c['particles_per_cell']
    if not diagdir.is_dir():raise ValueError('diagnostic directory missing')
    observed=set(p.name for p in diagdir.glob('*.h5'))
    if observed != set(EXPECTED):raise ValueError(f'HDF5 file set mismatch missing={set(EXPECTED)-observed}, extra={observed-set(EXPECTED)}')
    sizes={name:(diagdir/name).stat().st_size for name in EXPECTED}
    cap=protocol['resource_limits']
    if any(s>cap['max_single_file_bytes'] for s in sizes.values()) or sum(sizes.values())>cap['max_hdf5_bytes']:
        raise ValueError('HDF5 output exceeded resource bounds')
    with (h5py.File(diagdir/'EM_B.h5') as bf,
          h5py.File(diagdir/'EM_E.h5') as ef,
          h5py.File(diagdir/'ions_charge_density.h5') as nf,
          h5py.File(diagdir/'ions_charge_flux.h5') as jf,
          h5py.File(diagdir/'ions_pop_protons_domain.h5') as pf):
        handles={'EM_B.h5':bf,'EM_E.h5':ef,'ions_charge_density.h5':nf,'ions_charge_flux.h5':jf,'ions_pop_protons_domain.h5':pf}
        for name,h in handles.items(): verify_hdf5_root(h,c,name)
        for attr,v in [('pop_charge',c['ion_charge']),('pop_mass',c['ion_mass'])]:
            if attr not in pf.attrs or abs(float(np.asarray(pf.attrs[attr]).item())-v)>1e-14:
                raise ValueError('population mass/charge mismatch')
        rows=[]
        for t in c['snapshot_times']:
            b=one_patch(bf,t,c,'EM_B.h5')
            e=one_patch(ef,t,c,'EM_E.h5')
            den=one_patch(nf,t,c,'ions_charge_density.h5')
            flu=one_patch(jf,t,c,'ions_charge_flux.h5')
            p=one_patch(pf,t,c,'ions_pop_protons_domain.h5')
            fields={key:extract_component(ds,CENTERING[key],n) for patch in (b,e,den,flu) for key,ds in patch.items()}
            for name in ('EM_B_x','EM_B_y','EM_B_z','EM_E_x','EM_E_y','EM_E_z','charge_density','charge_flux_x','charge_flux_y','charge_flux_z'):
                if name not in fields:raise ValueError(f'missing required dataset {name}')
            for key in EXPECTED['ions_pop_protons_domain.h5']:
                if key not in p:raise ValueError(f'missing particle dataset {key}')
            w=np.asarray(p['weight'][...]);vel=np.asarray(p['v'][...]);idx=np.asarray(p['iCell'][...]);deltas=np.asarray(p['delta'][...])
            if w.dtype!=np.dtype('float64') or vel.dtype!=np.dtype('float64') or deltas.dtype!=np.dtype('float64') or idx.dtype.kind not in 'iu':
                raise ValueError('unexpected particle dtype')
            if w.shape!=(expectedmacro,1) or vel.shape!=(expectedmacro,3) or idx.shape!=(expectedmacro,1) or deltas.shape!=(expectedmacro,1):
                raise ValueError('particle count or array dimensions invalid')
            w=w[:,0];idx=idx[:,0];deltas=deltas[:,0]
            if not np.all(np.isfinite(w)) or np.min(w)<=0 or not np.all(np.isfinite(vel)) or not np.all(np.isfinite(deltas)):
                raise ValueError('invalid particle array')
            if np.any(idx<0) or np.any(idx>=n) or np.any(deltas<0) or np.any(deltas>1):
                raise ValueError('particle ownership or position out of bounds')
            Bx,By,Bz=(fields[f'EM_B_{key}'] for key in 'xyz')
            # Primal periodic boundary endpoint duplicated in HDF5, count exactly one period.
            ubx=dx*.5*(.5*Bx[0]**2+math.fsum(map(float,Bx[1:-1]**2))+.5*Bx[-1]**2)
            uby=dx*.5*math.fsum(map(float,By**2));ubz=dx*.5*math.fsum(map(float,Bz**2))
            ki=dx*.5*c['ion_mass']*math.fsum(map(float,w*np.sum(vel*vel,axis=1)))
            m=[dx*c['ion_mass']*math.fsum(map(float,w*vel[:,j])) for j in range(3)]
            totalN=dx*math.fsum(map(float,w))
            rho=fields['charge_density'];nrho=dx*(.5*float(rho[0])+math.fsum(map(float,rho[1:-1]))+.5*float(rho[-1]))/c['ion_charge']
            x=(np.arange(n)+.5)*dx
            mode=complex(np.mean((By+1j*Bz)*np.exp(-1j*2*math.pi*x/c['domain_length'])))
            inverse=complex(np.mean((By+1j*Bz)*np.exp(+1j*2*math.pi*x/c['domain_length'])))
            rows.append({'t':t,'particle_count':len(w),'N':totalN,'N_deposited':nrho,'P':m,'K':ki,'UBx':ubx,'UBy':uby,'UBz':ubz,'UB_transverse':uby+ubz,'U':math.fsum((ki,ubx,uby,ubz)),
                         'Bx_max_deviation':float(np.max(np.abs(Bx-c['B0']))),'Bx_periodic_mismatch':abs(float(Bx[0]-Bx[-1])),
                         'initial_wave_mode_real':mode.real,'initial_wave_mode_imag':mode.imag,'wave_amplitude':abs(mode),'wrong_helicity_amplitude':abs(inverse),
                         'max_abs_E':float(max(np.max(np.abs(fields[f'EM_E_{q}'])) for q in 'xyz'))})
    return rows,sizes


def screen(rows,protocol):
    c=protocol['case'];sc=protocol['registered_scales'];g={}
    ref=rows[0]; eps=np.finfo(np.float64).eps
    tauN=max(1e-12,8*eps*c['cells']*c['particles_per_cell'])
    g['G1_initial_perturbation']={'pass':abs(ref['wave_amplitude']-c['field_amplitude'])<=.02*c['field_amplitude'] and abs(ref['initial_wave_mode_imag'])<=.02*c['field_amplitude'] and ref['wrong_helicity_amplitude']<=.02*c['field_amplitude'],
                                  'observed_amplitude':ref['wave_amplitude'],'expected':c['field_amplitude'],'relative_tolerance':.02,'wrong_helicity_amplitude':ref['wrong_helicity_amplitude']}
    dN=[abs(r['N']-ref['N'])/sc['particle_equivalent'] for r in rows]
    g['G2_weighted_number']={'pass':max(dN)<=tauN and abs(ref['N']-sc['particle_equivalent'])<=1e-10 and max(abs(r['N']-r['N_deposited']) for r in rows)<=1e-8,
                             'max_relative_drift':max(dN),'registered_tolerance':tauN,'max_deposition_discrepancy':max(abs(r['N']-r['N_deposited']) for r in rows)}
    bx=max(max(r['Bx_max_deviation'],r['Bx_periodic_mismatch']) for r in rows)
    g['G3_longitudinal_Bx_1D']={'pass':bx<=1e-12,'max_deviation':bx,'registered_tolerance':1e-12}
    pm=[math.sqrt(math.fsum((r['P'][i]-ref['P'][i])**2 for i in range(3))) for r in rows]
    g['G4_momentum_screen']={'pass':max(pm)<=.05*sc['momentum'],'max_absolute_drift':max(pm),'max_scaled_drift':max(pm)/sc['momentum'],'registered_limit':.05*sc['momentum']}
    energy_delta=[math.fsum((r['K']-ref['K'],r['UBx']-ref['UBx'],r['UBy']-ref['UBy'],r['UBz']-ref['UBz'])) for r in rows]
    maxeng=max(abs(d) for d in energy_delta)
    g['G5_ideal_energy_screen']={'pass':maxeng<=.05*sc['transverse_magnetic_energy'],'max_absolute_excursion':maxeng,'max_scaled_excursion':maxeng/sc['transverse_magnetic_energy'],'registered_limit':.05*sc['transverse_magnetic_energy']}
    for row,p,ee in zip(rows,pm,energy_delta):
        row['deltaN']=row['N']-ref['N'];row['deltaP_norm']=p;row['deltaU_compensated']=ee
        row['deltaK']=row['K']-ref['K'];row['deltaUB_transverse']=row['UB_transverse']-ref['UB_transverse']
    gatepass=all(q['pass'] for q in g.values())
    return {'status':'BOUNDED_SINGLE_CASE_SCREEN_PASS' if gatepass else 'BOUNDED_SINGLE_CASE_SCREEN_FAIL',
            'not_full_conservation_qualification':True,'not_shielding_validation':True,
            'gates':g,'G6_exchange_exploratory':{'deltaK_end':rows[-1]['deltaK'],'deltaUBtrans_end':rows[-1]['deltaUB_transverse']},
            'G7_discrete_work':'NOT_TESTED_internal_Eavg_not_saved','rows':rows}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--protocol',required=True,type=Path)
    ap.add_argument('--diagdir',required=True,type=Path)
    ap.add_argument('--outdir',required=True,type=Path)
    args=ap.parse_args();args.outdir.mkdir(parents=True,exist_ok=True)
    result={'status':'INCONCLUSIVE_DATA_CONTRACT','not_full_conservation_qualification':True,'not_shielding_validation':True}
    try:
        doc=json.loads(args.protocol.read_text()); manifest=json.loads((args.diagdir/'p3b1d_perturbed_run.json').read_text())
        cfg=validate_config(doc,manifest)
        if manifest['protocol_sha256']!=digest(args.protocol):raise ValueError('protocol file sha differs from run provenance')
        rows,sizes=read_rows(args.diagdir,doc)
        result=screen(rows,doc)
        result['experiment_id']=doc['protocol_id'];result['source_sha']=PHARE_SHA
        result['research_sha']=manifest['research_commit'];result['protocol_sha256']=digest(args.protocol)
        result['total_hdf5_bytes']=sum(sizes.values());result['hdf5_sha256']={x:digest(args.diagdir/x) for x in EXPECTED}
    except Exception as e:
        result['error']=f'{type(e).__name__}: {e}'
    (args.outdir/'assessment.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    if 'rows' in result:
        with open(args.outdir/'timeseries.csv','w',newline='') as f:
            writer=csv.writer(f);writer.writerow(['time','N','N_deposited','P_x','P_y','P_z','K','UBx','UBy','UBz','U','deltaK','deltaUB_transverse','deltaU_compensated','Bx_dev','wave_amplitude'])
            for r in result['rows']:writer.writerow([r['t'],r['N'],r['N_deposited'],*r['P'],r['K'],r['UBx'],r['UBy'],r['UBz'],r['U'],r['deltaK'],r['deltaUB_transverse'],r['deltaU_compensated'],r['Bx_max_deviation'],r['wave_amplitude']])
    print(json.dumps({key:result[key] for key in ('status','gates','error') if key in result},indent=2))
    return 0 if result['status']=='BOUNDED_SINGLE_CASE_SCREEN_PASS' else 3

if __name__=='__main__':raise SystemExit(main())
