#!/usr/bin/env python3
"""Independent assessment of Phase 3B.1C raw B-field HDF5 diagnostics.

No PHARE import. It does not create simulation results. It reports INCONCLUSIVE
when no raw outputs are present. Bootstrap uses fixed RNG seed for reproducibility.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
import h5py

BOOTSTRAP_SEED=20261120
BOOTSTRAP_N=2500


def fields_for_case(path,n):
    times=[]; fields=[]
    with h5py.File(path,'r') as h5:
        for tkey in sorted(h5['t'].keys(),key=float):
            arr=np.full(n,np.nan+1j*np.nan,dtype=complex)
            counts=np.zeros(n,dtype=int)
            for p in h5['t'][tkey]['pl0'].values():
                off=int(p.attrs['lower'][0]);nc=int(p.attrs['nbrCells'][0]);
                y=p['EM_B_y'][...];z=p['EM_B_z'][...]
                ghost=(len(y)-nc)//2
                if (len(y)-nc)%2 or len(y)!=len(z): raise ValueError('bad staggered field lengths')
                arr[off:off+nc]=y[ghost:ghost+nc]+1j*z[ghost:ghost+nc]
                counts[off:off+nc]+=1
            if not np.all(counts==1):raise ValueError('incomplete physical cell coverage')
            times.append(float(tkey));fields.append(arr)
    result=np.stack(fields)
    if not np.isfinite(result).all():raise ValueError('non-finite raw magnetic field values')
    return np.asarray(times),result


def cold(k,branch):
    return k*(math.sqrt(1+k*k/4)+(1 if branch=='whistler' else -1)*k/2)

def empirical(k,branch,dx,order):
    kg=2/dx*math.sin(k*dx/2)
    di=1-(order+1)*dx*dx/12
    return kg*(math.sqrt(1+(kg*di)**2/4)+(1 if branch=='whistler' else -1)*kg*di/2)

def analyse_case(case,rawdir,seriesout):
    n=case['cells'];dx=case['dx_di'];order=case['interp_order']
    path=rawdir/f"dispersion_{case['id']}_o{order}"/'EM_B.h5'
    times,b=fields_for_case(path,n)
    if len(times)!=1001 or not np.allclose(times,np.arange(1001)*0.01,atol=1e-7):
        raise ValueError(f'{case["id"]}: incorrect diagnostic times')
    freq=np.fft.fft(b,axis=1)/n
    rows=[]
    for branch in ('aic','whistler'):
        for mode in case['modes']:
            k=2*math.pi*mode/6.4
            qidx=mode if branch=='aic' else -mode
            y=freq[:,qidx]; ph=np.unwrap(np.angle(y))
            slope, intercept=np.polyfit(times,ph,1)
            w=abs(float(slope));fitted=np.polyval([slope,intercept],times)
            if not math.isfinite(w) or w<=0 or not np.isfinite(y).all():
                raise ValueError('non-finite or non-positive mode measurement')
            mid=len(times)//2
            slope_a=np.polyfit(times[:mid+1],ph[:mid+1],1)[0]
            slope_b=np.polyfit(times[mid:],ph[mid:],1)[0]
            omega_theory=cold(k,branch)
            omega_empirical=empirical(k,branch,dx,order)
            end_ratio=abs(y[-1])/abs(y[0])
            rec=dict(case_id=case['id'],family=case['family'],N=n,dx_di=dx,
                 dt=case['dt_omega_ci_inv'],ppc=case['ppc'],order=order,seed=case['seed'],
                 branch=branch,mode=mode,k_di=k,omega_measured=w,omega_cold=omega_theory,
                 omega_empirical=omega_empirical,
                 signed_error_vs_cold=(w-omega_theory)/omega_theory,
                 signed_error_vs_empirical=(w-omega_empirical)/omega_empirical,
                 ratio_end_to_start=float(end_ratio),
                 ratio_min_to_start=float(abs(y).min()/abs(y[0])),
                 phase_rms_rad=float(np.sqrt(np.mean((ph-fitted)**2))),
                 half_frequency_disagreement_rel=float(abs(abs(slope_a)-abs(slope_b))/w))
            rows.append(rec)
            for j,t in enumerate(times):
                seriesout.writerow(dict(case_id=case['id'],branch=branch,mode=mode,
                                       t_omega_ci_inv=t,amplitude=abs(y[j]),
                                       unwrapped_phase=ph[j]))
    return rows, (times,b)

def bootstrap_grid(rows):
    rng=np.random.default_rng(BOOTSTRAP_SEED)
    out=[];fail=[];uncertain=[]
    for branch in ('aic','whistler'):
        for mode in (1,2,3,4):
            grp={n:np.array([r['omega_measured'] for r in rows if r['family']=='grid' and r['N']==n and r['branch']==branch and r['mode']==mode]) for n in (64,96,128)}
            if any(len(v)!=3 for v in grp.values()):raise ValueError('missing 3-seed grid observations')
            x=np.array([(6.4/n)**2 for n in (64,96,128)])
            means=np.array([grp[n].mean() for n in (64,96,128)])
            slope,intercept=np.polyfit(x,means,1)
            t=cold(2*math.pi*mode/6.4,branch)
            sims=np.empty(BOOTSTRAP_N)
            reversals=0
            for j in range(BOOTSTRAP_N):
                yy=np.array([rng.choice(grp[n],size=3,replace=True).mean() for n in (64,96,128)])
                sims[j]=np.polyfit(x,yy,1)[1]
                if abs(yy[-1]-t)>abs(yy[0]-t)+0.0025*t: reversals+=1
            lo,hi=np.quantile(sims,[0.025,0.975])
            limit=.01*t
            reversals_probability=reversals/BOOTSTRAP_N
            outcome='PASS' if lo >=t-limit and hi <=t+limit and reversals_probability<0.95 else 'FAIL' if hi<t-limit or lo>t+limit or reversals_probability>0.95 else 'INCONCLUSIVE'
            out.append(dict(branch=branch,mode=mode,k_di=2*math.pi*mode/6.4,
                            omega_cold=t,omega0_extrapolated=float(intercept),
                            bootstrap_95_low=float(lo),bootstrap_95_high=float(hi),
                            relative_bias=(intercept-t)/t,
                            reversal_bootstrap_probability=reversals_probability,
                            result=outcome,grid_N64=means[0],grid_N96=means[1],grid_N128=means[2]))
    return out

def save_csv(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def compare(raw_a,raw_b):
    ta,ba=raw_a;tb,bb=raw_b
    if not np.array_equal(ta,tb) or ba.shape!=bb.shape:return float('inf')
    return float(np.max(np.abs(ba-bb)))

def assess(cases,records,raw_grids,issues):
    byid={c['id']:c for c in cases}
    gates={}
    gate0=not issues and len(raw_grids)==18 and len(records)==144
    gates['G0_data_completeness']='PASS' if gate0 else 'INCONCLUSIVE'
    if not gate0:return gates,[],[]
    health=[r for r in records if r['ppc']==400]
    bad=[r for r in health if not (0.3<r['ratio_end_to_start']<3.0)
         or r['phase_rms_rad']>.25 or r['half_frequency_disagreement_rel']>.005]
    gates['G1_high_ppc_health']='PASS' if not bad else 'INCONCLUSIVE'
    anchor='grid_n64_s314159'
    diff=compare(raw_grids[anchor],raw_grids['repeat_anchor'])
    ars=[r for r in records if r['case_id']==anchor]
    brr={(r['branch'],r['mode']):r for r in records if r['case_id']=='repeat_anchor'}
    wdiff=max(abs(r['omega_measured']-brr[(r['branch'],r['mode'])]['omega_measured'])/r['omega_measured'] for r in ars)
    gates['G2_seeded_repeat']='PASS' if diff<=1e-8 and wdiff<=1e-6 else 'FAIL'
    grid=bootstrap_grid(records)
    gates['G3_spatial_convergence']='FAIL' if any(r['result']=='FAIL' for r in grid) else 'PASS' if all(r['result']=='PASS' for r in grid) else 'INCONCLUSIVE'
    timestep=[]
    for branch in ('aic','whistler'):
        for mode in (1,2,3,4):
            fit={r['dt']:r for r in records if (r['case_id']==anchor or r['family']=='timestep') and r['branch']==branch and r['mode']==mode}
            # Anchor provides dt=0.0005.
            gap=abs(fit[.0005]['omega_measured']-fit[.00025]['omega_measured'])/fit[.00025]['omega_measured']
            timestep.append(dict(branch=branch,mode=mode,relative_gap_dt0005_vs_dt00025=gap,
                                 result='PASS' if gap<=.0025 else 'INCONCLUSIVE_PENDING_REPLICATION'))
    gates['G4_timestep']='PASS' if all(r['result']=='PASS' for r in timestep) else 'INCONCLUSIVE'
    stats={}
    for ppc in (100,400):
        mode_var=[]
        for branch in ('aic','whistler'):
            for mode in (1,2,3,4):
                obs=[r['omega_measured']/r['omega_cold'] for r in records if r['N']==64 and r['order']==2 and r['dt']==.0005 and r['ppc']==ppc and r['branch']==branch and r['mode']==mode and r['case_id']!='repeat_anchor']
                if len(obs)!=3:raise ValueError(f'seed count unexpected, ppc={ppc}, {branch} {mode}: {len(obs)}')
                mode_var.append(float(np.std(obs,ddof=1)**2))
        stats[str(ppc)]=math.sqrt(float(np.mean(mode_var)))
    ratio=stats['400']/stats['100'] if stats['100']>0 else float('inf')
    independent_particle_seed_differences=[]
    for name in ('grid_n64_s271828','grid_n64_s161803'):
        independent_particle_seed_differences.append(compare(raw_grids['grid_n64_s314159'],raw_grids[name]))
    seed_variation_confirmed=all(x>1e-12 for x in independent_particle_seed_differences)
    gates['G5_noise_screen']='PASS' if seed_variation_confirmed and math.isfinite(ratio) and ratio<=.8 and gates['G1_high_ppc_health']=='PASS' else 'INCONCLUSIVE'
    return gates,grid,timestep,dict(grid_seeded_repeat_max_field_diff=diff,
              grid_seeded_repeat_max_rel_omega_diff=wdiff,
              seed_variation_confirmed=seed_variation_confirmed,seed_field_differences=independent_particle_seed_differences,
              pooled_noise_std_ppc100=stats['100'],pooled_noise_std_ppc400=stats['400'],pooled_noise_ratio_400_over_100=ratio,
              g1_unstable_modes=[r['case_id']+':'+r['branch']+str(r['mode']) for r in bad])

def charts(out,grid,timestep,records):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,5))
    for branch in ('aic','whistler'):
        for mode in (1,2,3,4):
            dat=[]
            for n in (64,96,128):
                sel=[r for r in records if r['family']=='grid' and r['N']==n and r['branch']==branch and r['mode']==mode]
                dat.append([100*np.mean([r['signed_error_vs_cold'] for r in sel]),100*np.std([r['signed_error_vs_cold'] for r in sel],ddof=1)])
            ax.errorbar([.1,6.4/96,.05],[v[0] for v in dat],yerr=[v[1] for v in dat],marker='o',label=f'{branch} m{mode}',capsize=2)
    ax.set(xlabel=r'Grid spacing $\Delta x/d_i$', ylabel='Signed error vs exact cold theory (%)',title='3B.1C seeded 3-level grid convergence (mean +/- 1 SD)')
    ax.grid(alpha=.3);ax.legend(ncol=2,fontsize=7);fig.tight_layout();fig.savefig(out/'spatial-convergence.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5))
    for branch in ('aic','whistler'):
        for mode in (1,2,3,4):
            sel=sorted([r for r in records if (r['family']=='timestep' or r['case_id']=='grid_n64_s314159') and r['branch']==branch and r['mode']==mode],key=lambda r:r['dt'])
            ax.plot([r['dt'] for r in sel],[100*r['signed_error_vs_cold'] for r in sel],marker='o',label=f'{branch} m{mode}')
    ax.set(xlabel=r'Time step $\Delta t\,\Omega_{ci}$',ylabel='Signed error vs cold theory (%)',title='3B.1C temporal sensitivity at fixed grid')
    ax.grid(alpha=.3);ax.legend(ncol=2,fontsize=7);fig.tight_layout();fig.savefig(out/'temporal-convergence.png',dpi=150);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True);p.add_argument('--raw-dir',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    cases=json.loads(a.manifest.read_text())['cases']
    a.output.mkdir(parents=True,exist_ok=True)
    records=[];raw={};issues=[]
    with (a.output/'per-mode-amplitude-phase-timeseries.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['case_id','branch','mode','t_omega_ci_inv','amplitude','unwrapped_phase']);writer.writeheader()
        for c in cases:
            try:
                path=a.raw_dir/f"dispersion_{c['id']}_o{c['interp_order']}"/'EM_B.h5'
                if not path.exists():raise FileNotFoundError(path)
                rr,bb=analyse_case(c,a.raw_dir,writer)
                records.extend(rr);raw[c['id']]=bb
            except Exception as exc:
                issues.append(f"{c['id']}: {type(exc).__name__}: {exc}")
    save_csv(a.output/'frequency-results.csv',records)
    got=assess(cases,records,raw,issues)
    if len(got)==4:
        gates,grid,timestep,extras=got
    else:
        gates,grid,timestep=got;extras={}
    save_csv(a.output/'grid-bootstrap.csv',grid);save_csv(a.output/'timestep-differences.csv',timestep)
    if all(v=='PASS' for v in gates.values()):verdict='PASS'
    elif any(v=='FAIL' for v in gates.values()):verdict='FAIL'
    else:verdict='INCONCLUSIVE'
    report=dict(status=verdict,meaning='Numerical qualification for tested 1D problem only; not Project Shield validation',
                gates=gates,extras=extras,issues=issues,measurements=len(records),
                expected_measurements=144,bootstrap_seed=BOOTSTRAP_SEED,bootstrap_draws=BOOTSTRAP_N)
    (a.output/'assessment.json').write_text(json.dumps(report,indent=2)+'\n')
    if len(got)==4:charts(a.output,grid,timestep,records)
    print(json.dumps(report,indent=2))
