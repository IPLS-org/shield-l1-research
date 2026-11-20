"""SYNTHETIC-ONLY offline analysis tests. NO real PHARE simulation or measured physics.

These tests generate PHARE-shaped HDF5 only to validate our own reader.
"""
import json
import math
from pathlib import Path
import shutil
import sys

import h5py
import numpy as np
import pytest

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'scripts'))
import p3b1d_perturbed_analyse as A
DOC=json.loads((HERE/'protocol'/'phase3b1d-one-perturbed-v0.4.json').read_text())
C=DOC['case']
ROOT=HERE/'protocol'/'phase3b1d-one-perturbed-v0.4.json'


def create_synthetic(dirpath):
    n=C['cells'];ppc=C['particles_per_cell'];dx=C['dx'];k=2*math.pi/C['domain_length'];delta=C['field_amplitude']
    times=C['snapshot_times']
    dirpath.mkdir(exist_ok=True,parents=True)
    ghost=4
    for name,fields in A.EXPECTED.items():
        with h5py.File(dirpath/name,'w') as f:
            f.attrs.update(dict(dimension=1,cell_width=[dx],origin=[0.0],periodicities=[1],domain_box=[n-1],interpOrder=C['interpolation_order'],layoutType='yee'))
            if name.endswith('_domain.h5'):
                f.attrs['pop_mass']=C['ion_mass'];f.attrs['pop_charge']=C['ion_charge']
            ts=f.create_group('t')
            for t in times:
                p=ts.create_group(f'{t:.10f}').create_group('pl0').create_group('p0#0')
                p.attrs['lower']=[0];p.attrs['upper']=[n-1];p.attrs['nbrCells']=[n];p.attrs['origin']=[0.0]
                if name.endswith('_domain.h5'):
                    p.attrs['ghosts']=2
                    arrays={
                      'weight':np.full((n*ppc,1),1/ppc,dtype=np.float64),
                      'v':np.zeros((n*ppc,3),dtype=np.float64),
                      'delta':np.full((n*ppc,1),.5,dtype=np.float64),
                      'iCell':np.repeat(np.arange(n,dtype=np.int32),ppc).reshape(-1,1)}
                    for field in fields:p.create_dataset(field,data=arrays[field])
                else:
                    for field in fields:
                        cent=A.CENTERING[field]
                        idx=np.arange(n+(1 if cent=='primal' else 0)+2*ghost,dtype=float)-ghost
                        x=(idx+(.5 if cent=='dual' else 0.0))*dx
                        if field=='EM_B_x': values=np.ones(len(x))
                        elif field=='EM_B_y':values=delta*np.cos(k*x)
                        elif field=='EM_B_z':values=delta*np.sin(k*x)
                        elif field=='charge_density':values=np.ones(len(x))
                        else:values=np.zeros(len(x))
                        d=p.create_dataset(field,data=np.asarray(values,dtype=np.float64));d.attrs['ghosts']=ghost
    return dirpath


@pytest.fixture(scope='module')
def data(tmp_path_factory):
    return create_synthetic(tmp_path_factory.mktemp('synthetic_test')/'phare_shaped_not_phare')


def test_exact_registered_case():
    assert C['solver_steps']==1000
    assert C['snapshot_times']==[x/10 for x in range(11)]
    assert DOC['registered_scales']['transverse_magnetic_energy']==.00032


def test_hdf5_reader_and_independent_hand_solution(data):
    rows,sizes=A.read_rows(data,DOC)
    assert len(rows)==11 and sum(sizes.values())<DOC['resource_limits']['max_hdf5_bytes']
    for row in rows:
        assert row['particle_count']==25600
        assert abs(row['N']-6.4)<2e-14
        assert abs(row['N_deposited']-6.4)<2e-14
        assert abs(row['K'])<1e-30
        # Direct independent analytic circular transverse energy for each physical sample.
        assert abs(row['UB_transverse']-6.4*.01**2/2)<1e-16
        assert abs(row['UBx']-6.4/2)<1e-14


def test_synthetic_screen_would_pass(data):
    rows,_=A.read_rows(data,DOC)
    res=A.screen(rows,DOC)
    assert res['status']=='BOUNDED_SINGLE_CASE_SCREEN_PASS'
    assert all(v['pass'] for v in res['gates'].values())
    # This proves our own code processes hand-crafted arrays; NOT that PHARE passes.


def test_tiny_number_drift_rejected(data):
    rows,_=A.read_rows(data,DOC)
    rows[-1]['N']+=1e-3
    res=A.screen(rows,DOC)
    assert not res['gates']['G2_weighted_number']['pass']


def test_incorrect_initial_mode_rejected(data):
    rows,_=A.read_rows(data,DOC)
    rows[0]['wave_amplitude']=0.0
    res=A.screen(rows,DOC)
    assert not res['gates']['G1_initial_perturbation']['pass']


def test_wrong_helicity_rejected(data):
    rows,_=A.read_rows(data,DOC)
    rows[0]['wrong_helicity_amplitude']=0.01
    assert not A.screen(rows,DOC)['gates']['G1_initial_perturbation']['pass']


def test_momentum_drift_rejected(data):
    rows,_=A.read_rows(data,DOC)
    rows[3]['P'][2]+=0.01
    assert not A.screen(rows,DOC)['gates']['G4_momentum_screen']['pass']


def test_energy_drift_rejected(data):
    rows,_=A.read_rows(data,DOC)
    rows[-1]['UBz']+=0.00008
    assert not A.screen(rows,DOC)['gates']['G5_ideal_energy_screen']['pass']


def test_Bx_violation_rejected(data):
    rows,_=A.read_rows(data,DOC)
    rows[5]['Bx_max_deviation']=3e-4
    assert not A.screen(rows,DOC)['gates']['G3_longitudinal_Bx_1D']['pass']


def test_data_contract_rejects_float32(data,tmp_path):
    bad=tmp_path/'copy';shutil.copytree(data,bad)
    with h5py.File(bad/'EM_B.h5','r+') as h:
        p=h['t']['0.0000000000']['pl0']['p0#0']
        a=p['EM_B_y'][...].astype('float32');attrs=dict(p['EM_B_y'].attrs)
        del p['EM_B_y'];d=p.create_dataset('EM_B_y',data=a)
        for k,v in attrs.items():d.attrs[k]=v
    with pytest.raises(ValueError,match='float64'):A.read_rows(bad,DOC)


def test_data_contract_rejects_missing_timestamp(data,tmp_path):
    bad=tmp_path/'copy';shutil.copytree(data,bad)
    with h5py.File(bad/'EM_E.h5','r+') as h:del h['t']['0.5000000000']
    with pytest.raises(ValueError,match='time groups'):A.read_rows(bad,DOC)


def test_data_contract_rejects_wrong_particle_ownership(data,tmp_path):
    bad=tmp_path/'copy';shutil.copytree(data,bad)
    with h5py.File(bad/'ions_pop_protons_domain.h5','r+') as h:
        d=h['t']['0.0000000000']['pl0']['p0#0']['iCell'];d[0,0]=64
    with pytest.raises(ValueError,match='ownership'):A.read_rows(bad,DOC)


def test_no_runtime_provenance_is_rejected():
    fake={'source':'SYNTHETIC_FIXTURE','case':C,'phare_commit':DOC['phare_commit'],'protocol_id':DOC['protocol_id']}
    with pytest.raises(AssertionError,match='not a PHARE runtime manifest'):A.validate_config(DOC,fake)


def test_compensated_energy_deltas_do_not_use_background_total():
    rows=[dict(N=6.4,N_deposited=6.4,P=[0,0,0],K=1e-5,UBx=3.2,UBy=.0002,UBz=.00012,UB_transverse=.00032,
               wave_amplitude=.01,initial_wave_mode_imag=0,wrong_helicity_amplitude=0,Bx_max_deviation=0,Bx_periodic_mismatch=0)]
    rows2=[dict(rows[0],K=1e-5+1e-6,UBz=.00012-1e-6,UB_transverse=.00032-1e-6)]
    res=A.screen(rows+rows2,DOC)
    assert abs(res['gates']['G5_ideal_energy_screen']['max_absolute_excursion'])<1e-17


def test_older_real_unperturbed_case_is_not_spuriously_accepted():
    # An unperturbed case cannot be called a deliberately perturbed science PASS.
    fake=[dict(N=6.4,N_deposited=6.4,P=[0,0,0],K=0,UBx=3.2,UBy=0,UBz=0,UB_transverse=0,wave_amplitude=0,
               initial_wave_mode_imag=0,wrong_helicity_amplitude=0,Bx_max_deviation=0,Bx_periodic_mismatch=0)]
    assert not A.screen(fake,DOC)['gates']['G1_initial_perturbation']['pass']
