#!/usr/bin/env python3
"""Static/mock verification only; not a PHARE plasma simulation."""
import contextlib
import csv
import io
import importlib.util
import json
import math
from pathlib import Path
import tempfile
from unittest.mock import patch
import sys
import numpy as np
import h5py

HERE=Path(__file__).resolve().parent
MATRIX=HERE.parent/'protocol/phase3b1c-matrix.json'
manifest=json.loads(MATRIX.read_text())
assert len(manifest['cases'])==18
assert len(set(c['id'] for c in manifest['cases']))==18
assert manifest['initialisation']=='continuum_eigenvector'

sys.path.insert(0,str(HERE))
import p3b1c_execute as exe
import p3b1c_analyse as ana

# Pure reference calculations independently checked via analytic Hall splitting.
for k in [2*math.pi*m/6.4 for m in range(1,5)]:
    assert math.isclose(ana.cold(k,'whistler')-ana.cold(k,'aic'),k*k,rel_tol=1e-14)
for c in manifest['cases']:
    assert c['dt_omega_ci_inv']*c['dump_every']==0.01 or math.isclose(c['dt_omega_ci_inv']*c['dump_every'],0.01)
    assert math.isclose(c['cells']*c['dx_di'],6.4)

# A mock upstream dispersion.py confirms the driver uses identical continuum
# B(x) and u(x) initial functions independent of grid, order, ppc and seed.
# It deliberately does NOT pretend to emulate PHARE C++/MPI numerics.
fake='''
import numpy as np
from pathlib import Path
from types import SimpleNamespace
AMPLITUDE=0.01
RUNS={'high_k':{'modes':(1,2,3,4)}}
ph=SimpleNamespace(global_vars=SimpleNamespace(sim='ACTIVE'))
reference=None

def startMPI(): return None

def omega_aic(k):return k*(np.sqrt(1+k*k/4)-k/2)
def omega_whistler(k):return k*(np.sqrt(1+k*k/4)+k/2)
def wave_modes(lx,dl,order,modes):raise AssertionError('driver did not override wave_modes')
def k_grid(k,dl): raise AssertionError('driver did not override k_grid')
def diag_dir(name,order):return f'dispersion_{name}_o{order}'
def config(name,interp_order,cells,dl,modes,ppc,time_step,final_time,dump_every,seed):
    global reference
    assert abs(cells*dl-6.4)<1e-12
    aic,whistler=wave_modes(cells*dl,dl,interp_order,modes)
    phases=np.random.default_rng(0).uniform(0,2*np.pi,2*len(modes))
    x=np.linspace(0,6.4,25,endpoint=False)
    modes_all=[(q,w,p) for (_,q,w),p in zip(aic+whistler,phases)]
    b=np.array([sum(AMPLITUDE*np.exp(1j*(q*xx+p)) for q,w,p in modes_all) for xx in x])
    u=np.array([sum(-np.sign(q)*k_grid(abs(q),dl)/w*AMPLITUDE*np.exp(1j*(q*xx+p)) for q,w,p in modes_all) for xx in x])
    if reference is None:reference=(b,u)
    else:
        assert np.allclose(b,reference[0],atol=1e-12),name
        assert np.allclose(u,reference[1],atol=1e-12),name
    return SimpleNamespace(name=name,interp_order=interp_order,dry_run=False)
class Simulator:
    def __init__(self,sim):self.sim=sim
    def run(self):
        d=Path(diag_dir(self.sim.name,self.sim.interp_order));d.mkdir(exist_ok=True)
        (d/'EM_B.h5').write_bytes(b'SIMULATED MOCK FILE NOT PLASMA')
        return self
    def reset(self):return None
'''
with tempfile.TemporaryDirectory() as td:
    base=Path(td)
    root=base/'PHARE'
    source=root/'tests/functional/dispersion/dispersion.py'
    source.parent.mkdir(parents=True)
    source.write_text(fake)
    out=base/'output'
    args=['p3b1c_execute.py','--manifest',str(MATRIX),'--phare-root',str(root),'--output-root',str(out)]
    with patch.object(sys,'argv',args),patch('subprocess.check_output',return_value=exe.PIN+'\n'):
        with contextlib.redirect_stdout(io.StringIO()) as console:
            exe.main()
    assert len(list(out.glob('*.json')))==18
    for c in manifest['cases']:
        item=json.loads((out/(c['id']+'.json')).read_text())
        assert item['status']=='COMPLETED' and item['initialisation']=='continuum_eigenvector'
    assert 'RUN 18/18' in console.getvalue()

# Synthesize 1D magnetic fields with analytically known AIC and whistler
# frequencies. Test that the HDF5 parser and independent spectral phase-fit
# reconstruct those frequencies without accessing PHARE.
anchor=manifest['cases'][0]
n=anchor['cells'];length=6.4
with tempfile.TemporaryDirectory() as td:
    base=Path(td)
    d=base/f"dispersion_{anchor['id']}_o{anchor['interp_order']}"
    d.mkdir()
    path=d/'EM_B.h5'
    x=(np.arange(n)+.5)*length/n
    with h5py.File(path,'w') as h5:
        h5t=h5.create_group('t')
        for j,t in enumerate(np.arange(1001)*0.01):
            b=np.zeros(n,dtype=complex)
            for branch,sign in [('aic',1),('whistler',-1)]:
                for m in (1,2,3,4):
                    k=2*math.pi*m/length
                    w=ana.cold(k,branch)
                    # AIC waves +q and positive frequency; whistler waves -q and negative freq
                    b+=0.01*np.exp(1j*(sign*k*x-sign*w*t))
            tgrp=h5t.create_group(f'{t:.8f}')
            pl0=tgrp.create_group('pl0')
            patch=pl0.create_group('p0')
            patch.attrs['lower']=[0]
            patch.attrs['nbrCells']=[n]
            patch.create_dataset('EM_B_y',data=b.real)
            patch.create_dataset('EM_B_z',data=b.imag)
    sio=io.StringIO()
    writer=csv.DictWriter(sio,fieldnames=['case_id','branch','mode','t_omega_ci_inv','amplitude','unwrapped_phase'])
    rec, (times, bs)=ana.analyse_case(anchor,base,writer)
    assert len(times)==1001 and bs.shape==(1001,n) and len(rec)==8
    for r in rec:
        assert abs(r['omega_measured']-r['omega_cold'])<1e-12,(r['branch'],r['mode'],r['omega_measured'],r['omega_cold'])
        assert r['phase_rms_rad']<1e-12
    with h5py.File(path,'r+') as h5:
        h5['t']['0.00000000']['pl0']['p0']['EM_B_y'][0]=float('nan')
    try:
        ana.fields_for_case(path,n)
        raise AssertionError('nonfinite field data incorrectly accepted')
    except ValueError as err:
        assert 'non-finite' in str(err)

print('P3B1C v0.2 mock/static preflight PASS: 18 cases, fixed-physics initialisation, 8 analytic synthetic HDF5 modes, nonfinite rejection.')
print('No compiled PHARE solver or new physical plasma experiment was executed.')
