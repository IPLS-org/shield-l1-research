#!/usr/bin/env python3
"""Execute the pre-registered 18 PHARE dispersion cases. DO NOT run until reviewed.

Requires compiled PHARE at pinned revision. --list-only does not import PHARE.
Outputs one uniquely named HDF5 diagnostic directory per case.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

PIN='7e5c72c1e5eec5452c8e3cfd000ce26c01e084ef'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--phare-root',type=Path)
    p.add_argument('--output-root',type=Path)
    p.add_argument('--list-only',action='store_true')
    a=p.parse_args()
    doc=json.loads(a.manifest.read_text())
    assert doc['phare_commit']==PIN
    cases=doc['cases']
    assert len(cases)==18 and len(set(x['id'] for x in cases))==18
    if doc.get('initialisation')!='continuum_eigenvector':
        raise ValueError('Expected predeclared continuum initialisation')
    for c in cases:
        assert all(ch.isalnum() or ch in '_-' for ch in c['id'])
        assert abs(c['cells']*c['dx_di']-6.4)<1e-12
        assert abs(c['dt_omega_ci_inv']*c['dump_every']-0.01)<1e-12
        assert c['modes']==[1,2,3,4] and c['final_time_omega_ci_inv']==10.0
        assert c['B0']==1.0 and c['density']==1.0 and c['electron_temperature']==0.0
        assert c['amplitude']==0.01 and c['vthermal_alfven']==0.001
    if a.list_only:
        for c in cases:
            print(f"{c['id']:<31} N={c['cells']:3} dx={c['dx_di']:.8f} dt={c['dt_omega_ci_inv']:.6f} ppc={c['ppc']:3} order={c['interp_order']} seed={c['seed']}")
        print(f'{len(cases)} cases; 8 frequencies per case; no simulations run')
        return
    if not a.phare_root or not a.output_root:
        p.error('--phare-root and --output-root are required to execute')
    root=a.phare_root.resolve()
    dest=a.output_root.resolve()
    source=root/'tests/functional/dispersion/dispersion.py'
    assert source.exists(), source
    dest.mkdir(parents=True,exist_ok=True)
    # Verify exact source by the entire Git revision before any simulation.
    import subprocess
    found=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()
    if found!=PIN: raise RuntimeError(f'PHARE revision mismatch: {found}')
    spec=importlib.util.spec_from_file_location('pinned_phare_dispersion',source)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.AMPLITUDE!=0.01: raise RuntimeError('Unexpected upstream wave amplitude')
    if tuple(module.RUNS['high_k']['modes'])!=(1,2,3,4): raise RuntimeError('Upstream modes changed')
    # Fixed-physics continuum initialisation: the upstream benchmark uses
    # k_grid and an empirical di_eff to construct a discrete-mode velocity.
    # That makes the *physical* initial bulk velocity vary with dx and order.
    # For convergence studies we retain upstream PHARE's configuration and
    # simulation machinery, but replace those two pure helper functions ONLY
    # in this loaded module. B(x), u(x), phase, temperature and density are
    # thereby identical continuum functions for every numerical setting.
    # No empirical correction is used to initialise the physical waves.
    import math
    def physical_wave_modes(lx, dl, interp_order, modes):
        k0=2*math.pi/lx
        aic=[(m,m*k0,module.omega_aic(m*k0)) for m in modes]
        whistler=[(m,-m*k0,-module.omega_whistler(m*k0)) for m in modes]
        return aic,whistler
    module.wave_modes=physical_wave_modes
    module.k_grid=lambda k, dl: k
    record_source={"initialisation":"continuum_eigenvector", "continuum_wave_functions":"u=-q/omega_cold*b", "phase_seed":0}
    module.startMPI()
    orig=os.getcwd()
    os.chdir(dest)
    try:
        for ix,c in enumerate(cases,1):
            name=c['id']
            outdir=dest / module.diag_dir(name,c['interp_order'])
            if outdir.exists():
                raise FileExistsError(f'Case directory must be fresh: {outdir}')
            print(f'RUN {ix}/{len(cases)} {name}',flush=True)
            started=time.monotonic()
            record={'case':c,'phare_commit':PIN,'upstream_dispersion_sha256':sha(source),
                    'status':'STARTED','case_index':ix,**record_source,'github_sha':os.getenv('GITHUB_SHA'),
                    'github_run_id':os.getenv('GITHUB_RUN_ID')}
            (dest/f'{name}.json').write_text(json.dumps(record,indent=2)+'\n')
            sim=module.config(name=name,interp_order=c['interp_order'],cells=c['cells'],
                              dl=c['dx_di'],modes=tuple(c['modes']),ppc=c['ppc'],
                              time_step=c['dt_omega_ci_inv'],final_time=c['final_time_omega_ci_inv'],
                              dump_every=c['dump_every'],seed=c['seed'])
            module.Simulator(sim).run().reset()
            module.ph.global_vars.sim=None
            if sim.dry_run: raise RuntimeError('PHARE_DRY_RUN was enabled: no plasma execution')
            diagnostic=outdir/'EM_B.h5'
            if not diagnostic.is_file(): raise FileNotFoundError(diagnostic)
            record.update(status='COMPLETED',elapsed_seconds=time.monotonic()-started,
                          diagnostic_sha256=sha(diagnostic))
            (dest/f'{name}.json').write_text(json.dumps(record,indent=2)+'\n')
    finally:
        os.chdir(orig)

if __name__=='__main__':
    main()
