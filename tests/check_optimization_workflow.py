"""Optimized/reference E1, E2 and variable-k training and audit comparisons."""
import json, os, sys, tempfile
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from run_adaptive import run
from checkpoint_io import inspect
from compare_adaptive import near
exe=Path(sys.argv[1]).resolve()
def records(p): return [json.loads(x) for x in p.read_text().splitlines()]
with tempfile.TemporaryDirectory() as temporary:
    root=Path(temporary)
    for problem,k in [('E1',16),('E2',8),('E2',32),('E2',128)]:
        c=json.loads((ROOT/'configs/adaptive'/f'{problem.lower()}_main.json').read_text())
        # This recovery comparison needs every state's checkpoint. The main
        # production presets retain restart points only every five cycles.
        c.update(wavenumber=k,level=4,gap=2,ell=1,maximum_ell=2,ell_mode='fixed',checkpoint_interval_cycles=0,
                 cycles=1,threads=2,member_ids=[0,1,2],training_ids=[0,1],radius=10.,rank_cap=3,extra_checks=[])
        before=root/f'{problem}-{k}-reference';after=root/f'{problem}-{k}-optimized'
        with patch.dict(os.environ,{'ALOD_REFERENCE_EXECUTION':'1'}):m=run(c,before,exe,audit_threads=2)
        n=run(c,after,exe,audit_threads=2)
        assert m['audit_complete'] and n['audit_complete']
        a=[r for r in records(before/'solver.jsonl') if r['kind']=='accepted']
        b=[r for r in records(after/'solver.jsonl') if r['kind']=='accepted']
        assert len(a)==len(b)==3
        for x,y in zip(a,b):
            for key in ['coarse_marks','reference_marks','coarse_fingerprint','reference_fingerprint','rank','ell','training_stop']:
                assert x.get(key)==y.get(key),(problem,k,key,x.get(key),y.get(key))
            for key in ['eta','targets','frozen_scales','coarse_mean_total','reference_mean_total']:
                if x.get(key) is not None:near(x[key],y[key],relative=1e-8,absolute=1e-11)
        def samples(folder):return sorted([r for f in (folder/'audits').glob('*/samples.jsonl') for r in records(f) if r['kind']=='sample'],key=lambda r:(r['state_id'],r['sample']))
        a,b=samples(before),samples(after);assert len(a)==len(b)==9
        for x,y in zip(a,b):
            for key in ['e','f','g','E','F','G','exact_norm','energy']:
                if x.get(key) is not None:near(x[key],y[key],relative=1e-8,absolute=1e-11)
        assert len(list((after/'checkpoints/meshes').glob('*.bin')))==3
        cp,meta=inspect(after,exe);assert meta['format']==3
        if problem=='E1':
            blob=after/'checkpoints/meshes'/meta['geometry_file']
            damaged=bytearray(blob.read_bytes());damaged[20]^=1;blob.write_bytes(damaged)
            try:inspect(after,exe)
            except ValueError:pass
            else:raise AssertionError('corrupt shared dependency accepted')
            repaired=run(c,after,exe,resume='auto',audit_threads=2)
            assert repaired['audit_complete']
            assert list((after/'checkpoints/meshes/recovery').iterdir())
        print(problem,k,'marking, training decisions and audit observables agree')
