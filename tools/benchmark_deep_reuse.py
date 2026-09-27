"""Alternating bounded E1/E2 comparisons against the prior optimized executable."""
import argparse, hashlib, json, os, statistics
from pathlib import Path
from unittest.mock import patch
from run_adaptive import run
from compare_adaptive import near

def records(path):
    return [json.loads(line) for line in path.read_text().splitlines()]

def phases(path):
    totals={};counters={}
    if not path.exists():return totals,counters
    for row in records(path):
        name=row['phase']
        if row['event']=='end':totals[name]=totals.get(name,0)+row['seconds']
        if row['event']=='counter':
            if name.startswith(('patch_cache_','riesz_cache_')):counters[name]=row['value']
            else:counters[name]=counters.get(name,0)+row['value']
    return totals,counters

def compare(a,b):
    a=[r for r in records(a/'solver.jsonl') if r['kind']=='accepted']
    b=[r for r in records(b/'solver.jsonl') if r['kind']=='accepted']
    assert len(a)==len(b)
    for x,y in zip(a,b):
        for key in ('coarse_marks','reference_marks','rank','ell','training_stop','working_rank'):
            assert x[key]==y[key],(key,x[key],y[key])
        for key in ('eta','targets','frozen_scales','coarse_mean_total','reference_mean_total'):
            near(x[key],y[key],relative=1e-8,absolute=1e-10)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--before',type=Path,required=True);p.add_argument('--after',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--repeats',type=int,default=2)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    if a.repeats<1:raise ValueError('positive repeats required')
    root=Path(__file__).resolve().parents[1]
    result=dict(binary_sha256={name:hashlib.sha256(getattr(a,name).read_bytes()).hexdigest() for name in ('before','after')},runs=[])
    for problem,level,radius in [('E1',10,.6),('E2',8,10.)]:
        c=json.loads((root/'configs/adaptive'/f'{problem.lower()}_main.json').read_text())
        c.update(level=level,gap=4,cycles=1,state_limit=3,threads=32,ell=2,maximum_ell=2,ell_mode='fixed',
                 member_ids=list(range(16)),training_ids=list(range(16)),radius=radius,extra_checks=[])
        for repeat in range(a.repeats+1):
            folders={}
            for name in (('before','after') if repeat%2==0 else ('after','before')):
                folder=a.output/f'{problem}-{repeat}-{name}';folders[name]=folder
                env={'ALOD_AUDIT_ONESHOT':'1','ALOD_PATCH_CACHE_BYTES':'0'} if name=='before' else {}
                with patch.dict(os.environ,env):m=run(c,folder,getattr(a,name),audit_threads=32,audit_drain_workers=2)
                if not m['audit_complete']:raise RuntimeError(m.get('audit_failures'))
                adaptive,counters=phases(folder/'timings.jsonl');audits={};memory=0;pids=set()
                for file in (folder/'audits').glob('*/run.json'):
                    audit=json.loads(file.read_text());memory+=audit.get('memory_snapshot',False)
                    if audit.get('native_worker_pid'):pids.add(audit['native_worker_pid'])
                    times,counts=phases(file.parent/'timings.jsonl')
                    for key,value in times.items():audits[key]=audits.get(key,0)+value
                row=dict(problem=problem,version=name,warmup=repeat==0,repeat=repeat,
                    wall_seconds=m['wall_seconds'],solver_seconds=m['solver_wall_seconds'],drain_seconds=m['audit_drain_seconds'],
                    adaptive_phases=adaptive,audit_phases=audits,counters=counters,memory_snapshots=memory,audit_processes=len(pids))
                result['runs'].append(row);print(json.dumps(row),flush=True)
            compare(folders['before'],folders['after'])
            def samples(folder):return sorted([r for f in (folder/'audits').glob('*/samples.jsonl') for r in records(f) if r['kind']=='sample'],key=lambda r:(r['state_id'],r['sample']))
            for x,y in zip(samples(folders['before']),samples(folders['after'])):
                for key in ('e','f','g','energy','exact_norm'):near(x[key],y[key],relative=1e-8,absolute=1e-10)
    result['medians']={}
    for problem in ('E1','E2'):
        result['medians'][problem]={}
        for version in ('before','after'):
            rows=[r for r in result['runs'] if r['problem']==problem and r['version']==version and not r['warmup']]
            result['medians'][problem][version]={key:statistics.median(r[key] for r in rows) for key in ('wall_seconds','solver_seconds','drain_seconds')}
            for kind in ('adaptive_phases','audit_phases'):
                result['medians'][problem][version][kind]={key:statistics.median(r[kind].get(key,0) for r in rows) for key in set().union(*(r[kind] for r in rows))}
    result['decisions_and_observables_match']=True
    (a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
