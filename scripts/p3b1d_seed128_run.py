#!/usr/bin/env python3
"""Prospectively declared ONE 128-cell alternative-seed PHARE case. NOT Shield validation."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1<<20), b''):
            h.update(b)
    return h.hexdigest()


def run(protocol_path: Path, output: Path):
    import numpy as np
    import pyphare.pharein as ph
    from pyphare.simulator.simulator import Simulator, startMPI

    doc = json.loads(protocol_path.read_text())
    cfg = doc['case']
    assert doc['protocol_id'] == 'P3B1D-SEED-128-20261120-v0.7'
    assert doc['phare_commit'] == '7e5c72c1e5eec5452c8e3cfd000ce26c01e084ef'
    assert doc['status'] == 'REVIEW_DRAFT_NOT_AUTHORISATION_TO_DISPATCH'
    refpath = protocol_path.parent / 'p3b1d-seed-reference-run44819929641.json'
    if sha256(refpath) != doc['baseline_reference_sha256']:
        raise RuntimeError('The immutable baseline reference is not the prospectively registered input')
    ref = json.loads(refpath.read_text())
    if ref['baseline_raw_artifact_sha256'] != doc['baseline_artifact_sha256']:
        raise RuntimeError('Baseline source artifact identity has changed')
    if ref['baseline_phare_commit'] != doc['phare_commit']:
        raise RuntimeError('Baseline PHARE version does not match')
    allowed = set(doc['changed_case_parameters_only'])
    assert allowed == {'name','seed'} and cfg['seed'] == 271828 and ref['baseline_case']['seed'] == 314159
    differences = {k for k in cfg if cfg[k] != ref['baseline_case'].get(k)}
    if differences != allowed or set(cfg) != set(ref['baseline_case']):
        raise RuntimeError(f'Changed physics or unexpected parameter differences: {differences}')
    assert cfg['snapshot_times'] == ref['baseline_case']['snapshot_times']
    assert cfg['name'] == 'P1_grid128_second_seed_271828'

    assert cfg['cells'] == 128 and cfg['dx'] == 0.05 and cfg['particles_per_cell'] == 400
    assert cfg['solver_steps'] == 2000 and cfg['time_step'] == 0.0005 and cfg['snapshot_every_steps'] == 200
    assert cfg['final_time'] == 1.0 and cfg['field_amplitude'] == 0.01
    assert cfg['electron_temperature'] == cfg['resistivity'] == cfg['hyper_resistivity'] == 0
    assert cfg['interpolation_order'] == 2 and cfg['seed'] == 271828
    assert doc['grid_control']['smallest_patch_size'] == doc['grid_control']['largest_patch_size'] == cfg['cells']
    assert doc['grid_control']['max_nbr_levels'] == 1
    if os.getenv('PHARE_DRY_RUN', '').lower() not in ('', '0', 'false'):
        raise RuntimeError('dry-run is not an actual plasma simulation')
    if output.exists() and any(output.iterdir()):
        raise RuntimeError(f'Refusing nonempty output directory: {output}')
    root = Path(os.environ['PHARE_SOURCE_ROOT'])
    found = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
    if found != doc['phare_commit']:
        raise RuntimeError('Wrong pinned PHARE revision')
    if abs(cfg['cells'] * cfg['dx'] - cfg['domain_length']) > 1e-12:
        raise RuntimeError('Inconsistent domain length')
    steps = np.arange(0,cfg['solver_steps']+1,cfg['snapshot_every_steps'],dtype=int)
    times = steps * cfg['time_step']
    if len(times) != len(cfg['snapshot_times']) or not np.allclose(times, cfg['snapshot_times'],atol=1e-12,rtol=0):
        raise RuntimeError('Snapshot timestamps do not match frozen protocol')
    startMPI()
    L = cfg['domain_length']; delta = cfg['field_amplitude']; k = 2*np.pi/L
    sim = ph.Simulation(time_step=cfg['time_step'],final_time=cfg['final_time'],
                        cells=cfg['cells'],dl=cfg['dx'],interp_order=cfg['interpolation_order'],
                        smallest_patch_size=cfg['cells'],largest_patch_size=cfg['cells'],
                        resistivity=0.0,hyper_resistivity=0.0,
                        diag_options={'format':'phareh5','options':{'dir':str(output),'mode':'overwrite'}})
    ph.MaxwellianFluidModel(
        bx=lambda x: 1.0+0*x,
        by=lambda x: delta*np.cos(k*x),
        bz=lambda x: delta*np.sin(k*x),
        protons={'charge':1.0,'mass':1.0,'nbr_part_per_cell':cfg['particles_per_cell'],
                 'density':lambda x: 1.0+0*x,
                 'vbulkx':lambda x: 0.0*x,'vbulky':lambda x: 0.0*x,'vbulkz':lambda x: 0.0*x,
                 'vthx':lambda x: 0.001+0*x,'vthy':lambda x: 0.001+0*x,'vthz':lambda x: 0.001+0*x,
                 'init':{'seed':cfg['seed']}})
    ph.ElectronModel(closure='isothermal', Te=0.0)
    for name in ('B','E'):
        ph.ElectromagDiagnostics(quantity=name,write_timestamps=times)
    for name in ('charge_density','charge_flux'):
        ph.FluidDiagnostics(quantity=name,write_timestamps=times)
    ph.ParticleDiagnostics(quantity='domain',population_name='protons',write_timestamps=times)
    Simulator(sim).run().reset()
    ph.global_vars.sim = None
    if sim.dry_run:
        raise RuntimeError('PHARE returned dry_run: no real physics execution')
    output.mkdir(parents=True,exist_ok=True)
    (output/'p3b1d_seed128_run.json').write_text(json.dumps({
       'source':'PHARE_runtime', 'case':cfg,'protocol_id':doc['protocol_id'],
       'protocol_sha256':sha256(protocol_path),'phare_commit':found,
       'baseline_reference_sha256':sha256(refpath), 'baseline_run_id':ref['baseline_run_id'], 'patch_control':doc['grid_control'],
       'research_commit':os.getenv('GITHUB_SHA'),
       'scope':'ONE 128-CELL ALTERNATIVE SEED. IDENTICAL MESH/PHYSICS; PARTICLE SAMPLING CHANGES; NO VARIANCE ESTIMATION OR SHIELD VALIDATION.'},indent=2,sort_keys=True)+'\n')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--protocol',type=Path,required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    args=parser.parse_args()
    run(args.protocol,args.outdir)

if __name__=='__main__':main()
