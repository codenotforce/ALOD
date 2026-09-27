"""Target stopping, resume, exact-only audits and frozen scientific controls."""
import json,sys,tempfile,os
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from run_adaptive import run,DEFAULT
from run_audit import audit
from scientific_diagnostics import run as diagnose
from compare_adaptive import near
exe=Path(sys.argv[1]).resolve()
def samples(folder):return [json.loads(x) for x in (folder/'samples.jsonl').read_text().splitlines() if json.loads(x)['kind']=='sample']
with tempfile.TemporaryDirectory() as directory:
    root=Path(directory)
    c=dict(DEFAULT,problem='E2',level=3,gap=2,ell=1,maximum_ell=2,ell_mode='fixed',cycles=2,
        radius=10.,rank_cap=4,member_ids=list(range(6)),training_ids=list(range(4)),threads=2,audit=True,
        exact_target=10.,audit_mode='exact')
    m=run(c,root/'target',exe)
    assert m['target_reached'] and m['accepted_states']==2 and m['audit_complete'] and not m['complete_horizon']
    assert m['termination_reason']=='exact_target' and not m['paper_complete']
    original=(root/'target/solver.jsonl').read_bytes()
    resumed=run(dict(c,cycles=3),root/'target',exe,resume='auto')
    assert resumed['target_reached'] and (root/'target/solver.jsonl').read_bytes()==original
    with patch.dict(os.environ,{'ALOD_CHECKPOINT_SYNC':'1','ALOD_INTEGRATION_COLD':'1'}):
        cold=run(c,root/'cold',exe)
    a=[json.loads(x) for x in original.splitlines() if json.loads(x)['kind']=='accepted']
    b=[json.loads(x) for x in (root/'cold/solver.jsonl').read_text().splitlines() if json.loads(x)['kind']=='accepted']
    for x,y in zip(a,b):
        for key in ('coarse_marks','reference_marks','rank','ell','target_reached'):assert x[key]==y[key]
        near(x['exact_target_ratio'],y['exact_target_ratio'])
    checkpoint=root/'target/checkpoints/latest'
    full=audit(checkpoint,root/'full',exe,threads=2)
    light=audit(checkpoint,root/'light',exe,threads=2,audit_mode='exact')
    stored=audit(checkpoint,root/'stored',exe,threads=2,audit_mode='exact',member_ids=[0,1])
    assert stored['timing']['accepted_values_reused'] and stored['timing']['reference_factorizations']==0
    assert full['timing']['reference_factorizations']==1 and light['timing']['reference_factorizations']==0
    for x,y in zip(samples(root/'full'),samples(root/'light')):
        near(x['E'],y['E']);assert all(y[k] is None for k in ('f','g','F','G','E_ref','reference_residual'))
    d=diagnose(checkpoint,root/'diagnostics',exe,['all'],ells=[1,2],ranks=[1,2],radii=[.6],levels=[1],threads=2,member_ids=[0,1])
    assert d['complete'] and len(d['completed'])==16 and len(d['skipped'])==1
    assert all(x['rank']==0 for x in samples(root/'diagnostics/rank-zero-0'))
    assert samples(root/'diagnostics/reference-0')[0]['reference_mesh_hash']!=samples(root/'full')[0]['reference_mesh_hash']
    e1=dict(c,problem='E1',exact_target=1e6,exact_scope='training_max')
    e=run(e1,root/'e1',exe);assert e['accepted_states']==1 and e['target_reached']
    r=diagnose(root/'e1/checkpoints/latest',root/'e1-diagnostics',exe,['residual'],threads=2,member_ids=[0]);assert r['complete']
    print('Exact target/resume, zero-reference audit, cold parity and all applicable scientific controls passed.')
