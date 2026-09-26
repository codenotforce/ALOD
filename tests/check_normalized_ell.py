"""Exercise initial checks, repeated promotions, cap reporting and resume."""
import json, math, sys, tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from run_adaptive import DEFAULT, run
from checkpoint_io import inspect, validate_resume_config

exe=Path(sys.argv[1]).resolve()
def events(p):return [json.loads(x) for x in (p/'solver.jsonl').read_text().splitlines()]
def check(p,c):
    rows=events(p);states=[r for r in rows if r['kind']=='accepted']
    checks=[r for r in rows if r['kind']=='ell_check']
    assert checks[0]['state_id']==0 and checks[0]['ell']==2
    for s in states:
        same=[r for r in checks if r['state_id']==s['state_id']]
        if same:
            assert same[-1]['ell']==s['ell'] and same[-1]['action']!='promote'
            assert same[-1]['criterion_satisfied'] or same[-1]['action']=='ell_cap_reached'
            assert s['ell_last_check']==s['state_id'] and s['theta'] is not None
            assert math.isclose(same[-1]['solution_scale'],s['frozen_scales'][0],rel_tol=1e-12)
            for a,b in zip(same,same[1:]):
                assert a['action']=='promote' and b['ell']==a['ell']+1
                assert a['mesh_revision']==b['mesh_revision']
        else:assert s['theta'] is None and s['theta_seconds']==0
    for r in checks:
        if r['eta']>0:
            assert math.isclose(r['ratio'],r['solution_scale']*r['theta']/r['eta'],rel_tol=2e-14)
    return checks

report={}
with tempfile.TemporaryDirectory() as td:
    root=Path(td)
    for problem in ['E1','E2']:
        c=dict(DEFAULT,problem=problem,level=8 if problem=='E1' else 6,gap=2,ell=2,
               cycles=3,threads=2,maximum_ell=4,member_ids=[0],training_ids=[0],
               ell_mode='lazy',ell_ratio_mode='solution_scaled',ell_threshold=1e-12,
               reference_theta=[.1,.1],radius=.6,rank_cap=4,
               maximum_nodes=200000,maximum_patch_entries=100000000,
               maximum_dense_entries=10000000)
        p=root/problem;run(c,p,exe,timeout=360)
        rows=check(p,c);first=[r for r in rows if r['state_id']==0]
        assert len(first)>=2,(problem,first)
        if problem=='E1':assert [r['ell'] for r in first]==[2,3,4]
        paused=root/(problem+'-resume')
        run(c,paused,exe,timeout=360,pause_state=0,pause_phase='ell')
        cp,meta=inspect(paused,exe)
        assert meta['phase']==2  # AfterEll; pending flag is exercised by the replay below.
        run(c,paused,exe,timeout=360,resume=cp)
        check(paused,c)
        a,b=events(p),events(paused)
        assert len(a)==len(b)
        for x,y in zip(a,b):
            for k in ('state_id','kind','ell','action','coarse_fingerprint','reference_fingerprint','rank'):
                assert x.get(k)==y.get(k),(problem,k)
            for k in ('theta','ratio','solution_scale'):
                if x.get(k) is not None:assert math.isclose(x[k],y[k],rel_tol=1e-7,abs_tol=1e-12)
        try:validate_resume_config(c,dict(c,ell_ratio_mode='raw'))
        except ValueError:pass
        else:raise AssertionError('mathematical policy change accepted on resume')
        report[problem]={'initial_ells':[r['ell'] for r in first],
                         'initial_actions':[r['action'] for r in first],
                         'same_mesh_recheck':True,'resume_verified':True}
print(json.dumps(report,indent=2))
