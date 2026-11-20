#!/usr/bin/env python3
"""Project Shield L1 Phase 3B.1D: ONE diagnostic qualification case only.

No physical conservation PASS/FAIL is inferred here. A diagnostic-contract
exception means INCONCLUSIVE pending static review of raw HDF5 artifacts.
Pinned PHARE source: 7e5c72c1e5eec5452c8e3cfd000ce26c01e084ef.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import traceback

PHARE_SHA = '7e5c72c1e5eec5452c8e3cfd000ce26c01e084ef'
EXPECTED = ['EM_B.h5','EM_E.h5','ions_charge_density.h5',
            'ions_charge_flux.h5','ions_pop_protons_domain.h5']
FIELDS = {
    'EM_B.h5': ('EM_B_x', 'EM_B_y', 'EM_B_z'),
    'EM_E.h5': ('EM_E_x', 'EM_E_y', 'EM_E_z'),
    'ions_charge_density.h5': ('charge_density',),
    'ions_charge_flux.h5': ('charge_flux_x','charge_flux_y','charge_flux_z'),
    'ions_pop_protons_domain.h5': ('weight','iCell','delta','v'),
}
TIMES = (0.0, 0.025, 0.05)
CELLS = 64
DL = 0.1
PPC = 100
DT = 0.001
T_FINAL = 0.05
SEED = 314159


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for buf in iter(lambda: f.read(1 << 20), b''):
            h.update(buf)
    return h.hexdigest()


def run(outdir: Path):
    # No dry-run false positives. Only run in the explicitly authorised job.
    if os.environ.get('PHARE_DRY_RUN') not in (None, '', '0', 'false', 'False'):
        raise RuntimeError('PHARE_DRY_RUN is set; cannot qualify real diagnostic output')
    if outdir.exists() and any(outdir.iterdir()):
        raise RuntimeError(f'Refusing nonempty diagnostic output directory: {outdir}')

    import numpy as np
    import pyphare.pharein as ph
    from pyphare.simulator.simulator import Simulator, startMPI
    import subprocess

    phare_root = Path(os.environ['PHARE_SOURCE_ROOT'])
    actual = subprocess.check_output(['git', '-C', str(phare_root), 'rev-parse', 'HEAD'], text=True).strip()
    if actual != PHARE_SHA:
        raise RuntimeError(f'Wrong PHARE SHA: {actual}, expected {PHARE_SHA}')

    startMPI()
    sim = ph.Simulation(
        time_step=DT, final_time=T_FINAL, cells=CELLS, dl=DL,
        interp_order=2, resistivity=0.0, hyper_resistivity=0.0,
        diag_options={'format':'phareh5','options':{'dir':str(outdir),'mode':'overwrite'}}
    )
    ph.MaxwellianFluidModel(
        bx=lambda x: 1.0 + 0.0*x,
        by=lambda x: 0.0*x,
        bz=lambda x: 0.0*x,
        protons={
            'charge':1.0, 'mass':1.0, 'nbr_part_per_cell':PPC,
            'density':lambda x: 1.0 + 0.0*x,
            'vbulkx':lambda x: 0.0*x,
            'vbulky':lambda x: 0.0*x,
            'vbulkz':lambda x: 0.0*x,
            'vthx':lambda x: 0.001 + 0.0*x,
            'vthy':lambda x: 0.001 + 0.0*x,
            'vthz':lambda x: 0.001 + 0.0*x,
            'init':{'seed':SEED},
        },
    )
    ph.ElectronModel(closure='isothermal',Te=0.0)
    for quantity in ('B','E'):
        ph.ElectromagDiagnostics(quantity=quantity,write_timestamps=np.array(TIMES))
    for quantity in ('charge_density','charge_flux'):
        ph.FluidDiagnostics(quantity=quantity,write_timestamps=np.array(TIMES))
    ph.ParticleDiagnostics(quantity='domain',population_name='protons',
                           write_timestamps=np.array(TIMES))
    Simulator(sim).run().reset()
    ph.global_vars.sim = None
    if sim.dry_run:
        raise RuntimeError('PHARE unexpectedly ran in dry-run mode')
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / 'p3b1d_smoke_run.json').write_text(json.dumps({
        'source':'PHARE_runtime', 'phare_commit': PHARE_SHA,
        'research_commit':os.environ.get('GITHUB_SHA'),
        'config':{'cells':CELLS,'dx':DL,'ppc':PPC,'dt':DT,'final_time':T_FINAL,
                  'seed':SEED, 'te':0.0, 'eta':0.0, 'hyper_eta':0.0,
                  'interpolation_order':2,'expected_snapshot_times':TIMES,
                  'species':['protons']},
        'scope':'DIAGNOSTIC ONLY - NO CONSERVATION PASS'
    },indent=2,sort_keys=True)+'\n')


def inspect(diagdir: Path, report: Path, times=TIMES, cells=CELLS,
            max_total_bytes=64*1024*1024, require_phare_provenance=False):
    import h5py
    import numpy as np

    summary = {'status':'INCONCLUSIVE_DIAGNOSTIC_CONTRACT',
               'not_conservation_qualification':True,
               'input_provenance':'UNVERIFIED_INPUT',
               'expected_times':list(times),
               'expected_files':EXPECTED,
               'files':{}, 'warnings':[], 'errors':[]}
    errors = summary['errors']
    try:
        if require_phare_provenance:
            manifest = diagdir / 'p3b1d_smoke_run.json'
            if not manifest.is_file():
                raise ValueError('No runtime provenance manifest: no claim of actual PHARE execution')
            provenance = json.loads(manifest.read_text())
            if provenance.get('source') != 'PHARE_runtime' or provenance.get('phare_commit') != PHARE_SHA:
                raise ValueError('Incorrect PHARE runtime provenance')
            if provenance.get('config',{}).get('cells') != cells:
                raise ValueError('PHARE runtime manifest case mismatch')
            summary['input_provenance']='PHARE_RUNTIME_MANIFEST_PRESENT_NOT_INDEPENDENTLY_ATTESTED'
        actual_files = sorted(p.name for p in diagdir.glob('*.h5'))
        missing = sorted(set(EXPECTED) - set(actual_files))
        if missing:
            raise ValueError('Required HDF5 files absent: '+', '.join(missing))
        total_bytes = sum((diagdir / f).stat().st_size for f in actual_files)
        summary['total_hdf5_bytes'] = total_bytes
        if total_bytes > max_total_bytes:
            raise ValueError('HDF5 output exceeds hard post-run contract limit')
        for name in EXPECTED:
            path = diagdir / name
            info = {'sha256':sha256(path),'size':path.stat().st_size,
                    'times':[], 'dtypes':{},'patches':{},'snapshots':[]}
            with h5py.File(path,'r') as h:
                for attr in ('dimension','cell_width','origin','periodicities',
                             'domain_box','interpOrder','layoutType'):
                    if attr not in h.attrs:
                        raise ValueError(f'{name}: missing root attribute {attr}')
                if int(np.asarray(h.attrs['dimension']).item()) != 1:
                    raise ValueError(f'{name}: not one-dimensional')
                dx = np.asarray(h.attrs['cell_width']).reshape(-1)
                if dx.size != 1 or not np.isclose(dx[0],DL,atol=1e-12,rtol=0):
                    raise ValueError(f'{name}: inconsistent dx')
                if not bool(np.asarray(h.attrs['periodicities']).reshape(-1)[0]):
                    raise ValueError(f'{name}: not periodic')
                if int(np.asarray(h.attrs['interpOrder']).item()) != 2:
                    raise ValueError(f'{name}: incorrect interpolation order')
                if name.endswith('_domain.h5'):
                    for attr in ('pop_mass','pop_charge'):
                        if attr not in h.attrs:
                            raise ValueError(f'{name}: missing {attr}')
                    summary['particle_species']={
                        'mass':float(np.asarray(h.attrs['pop_mass']).item()),
                        'charge':float(np.asarray(h.attrs['pop_charge']).item())}
                if 't' not in h:
                    raise ValueError(f'{name}: missing time group')
                recorded = sorted([float(v) for v in h['t'].keys()])
                info['times']=recorded
                if len(recorded) != len(times) or any(
                    not np.isclose(x,y,atol=5e-9,rtol=0) for x,y in zip(recorded,times)):
                    raise ValueError(f'{name}: expected {times}, got {recorded}')
                for tkey in sorted(h['t'].keys(), key=float):
                    group=h['t'][tkey]
                    if list(group.keys()) != ['pl0']:
                        raise ValueError(f'{name} at {tkey}: unexpected levels {list(group.keys())}')
                    patches=group['pl0']; owned=[]; partcounts=0; weight_sum=0.0
                    for patch_name,pdata in patches.items():
                        if 'lower' not in pdata.attrs or 'upper' not in pdata.attrs:
                            raise ValueError(f'{name}: patch {patch_name} missing lower/upper')
                        lower=int(np.asarray(pdata.attrs['lower']).reshape(-1)[0]);upper=int(np.asarray(pdata.attrs['upper']).reshape(-1)[0])
                        if upper < lower:
                            raise ValueError(f'{name}: malformed patch bounds')
                        owned.extend(range(lower,upper+1))
                        for field in FIELDS[name]:
                            if field not in pdata:
                                raise ValueError(f'{name}: {patch_name} missing dataset {field}')
                            ds=pdata[field]
                            if field=='iCell':
                                if ds.dtype.kind not in ('i','u'):
                                    raise ValueError(f'{name}: integer cells not stored as integers')
                            elif ds.dtype != np.dtype('float64'):
                                raise ValueError(f'{name}: {field} has {ds.dtype}, expected float64')
                            info['dtypes'][field]=str(ds.dtype)
                            arr=ds[()]
                            if not np.all(np.isfinite(arr)):
                                raise ValueError(f'{name}: nonfinite {field}')
                            if field=='weight':
                                partcounts += arr.size
                                weight_sum += float(np.sum(arr,dtype=np.float64))
                                if np.any(arr <= 0):
                                    raise ValueError(f'{name}: nonpositive particle weight')
                        if name.endswith('_domain.h5'):
                            w=pdata['weight'][:].reshape(-1)
                            v=pdata['v'][:].reshape(-1)
                            ii=pdata['iCell'][:].reshape(-1)
                            if v.size != 3*w.size or ii.size != w.size:
                                raise ValueError(f'{name}: inconsistent particle arrays')
                            if np.any(ii < lower) or np.any(ii > upper):
                                raise ValueError(f'{name}: particle outside owned patch')
                    if sorted(owned) != list(range(cells)):
                        raise ValueError(f'{name} at {tkey}: incomplete or duplicated physical cells')
                    item={'timestamp':float(tkey),'owned_cells':len(owned)}
                    if name.endswith('_domain.h5'):
                        item['macro_particles']=partcounts
                        item['dx_weighted_sum']=DL*weight_sum
                        if partcounts != PPC*cells:
                            summary['warnings'].append(f'{name} at {tkey}: particle count changed from initial target')
                    info['snapshots'].append(item)
            summary['files'][name]=info
        # Diagnostic data structure and nominal timing pass only, not solver accuracy.
        summary['status']='DIAGNOSTIC_SCHEMA_OK_NOT_PHYSICS_VALIDATED'
    except Exception as err:
        errors.append(str(err))
        summary['exception_type']=type(err).__name__
    report.parent.mkdir(parents=True,exist_ok=True)
    report.write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':summary['status'],'errors':errors,
                      'total_bytes':summary.get('total_hdf5_bytes')},indent=2))
    return 0 if not errors else 3


def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest='cmd',required=True)
    r=sub.add_parser('run');r.add_argument('--outdir',type=Path,required=True)
    s=sub.add_parser('inspect');s.add_argument('--diagdir',type=Path,required=True)
    s.add_argument('--report',type=Path,required=True)
    s.add_argument('--times',type=float,nargs=3,default=TIMES)
    s.add_argument('--cells',type=int,default=CELLS)
    s.add_argument('--require-phare-provenance',action='store_true')
    a=p.parse_args()
    if a.cmd=='run':run(a.outdir)
    else:sys.exit(inspect(a.diagdir,a.report,times=a.times,cells=a.cells,
                          require_phare_provenance=a.require_phare_provenance))

if __name__=='__main__':main()
