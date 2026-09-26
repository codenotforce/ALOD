"""Alternating workflow benchmark with numerical and marking equivalence gates."""
import argparse, json, math, statistics, os
from pathlib import Path
from run_adaptive import run
from checkpoint_io import digest

def records(path):return [json.loads(x) for x in path.read_text().splitlines()]

def compare(a,b,key=''):
    if key.endswith('_seconds') or 'hash' in key:return
    if isinstance(a,dict):
        assert set(a)==set(b),(key,'fields')
        for k in a:compare(a[k],b[k],k)
    elif isinstance(a,list):
        assert len(a)==len(b),key
        for x,y in zip(a,b):compare(x,y,key)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-10),(key,a,b)
    else:assert a==b,(key,a,b)

def phase_sums(paths):
    out={}
    for path in paths:
        for r in records(path):
            if r['event']=='end':out[r['phase']]=out.get(r['phase'],0)+r['seconds']
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--before',type=Path,required=True);p.add_argument('--after',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--problem',choices=['E1','E2'],default='E1')
    p.add_argument('--level',type=int,default=8);p.add_argument('--gap',type=int,default=4)
    p.add_argument('--threads',type=int,default=32);p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--wavenumber',type=float,default=16)
    p.add_argument('--radius',type=float)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    c=json.loads(Path(f'configs/adaptive/{args.problem.lower()}_main.json').read_text())
    c.update(level=args.level,gap=args.gap,threads=args.threads,cycles=1,member_ids=list(range(16)),training_ids=list(range(16)),wavenumber=args.wavenumber,extra_checks=[])
    if args.radius is not None:c['radius']=args.radius
    report=dict(config=c,operational_environment={key:os.environ.get(key) for key in ['ALOD_PATCH_CACHE_BYTES','ALOD_SOURCE_MOMENT_BYTES','ALOD_REFERENCE_EXECUTION']},runs=[],binary_sha256={name:digest(getattr(args,name)) for name in ['before','after']})
    reference=None
    for iteration in range(args.repeats+1):
        for version in (['before','after'] if iteration%2==0 else ['after','before']):
            folder=args.output/f'{iteration}-{version}'
            m=run(c,folder,getattr(args,version),timeout=3600,audit_threads=args.threads,audit_workers=1)
            assert m['audit_complete'] and m['solver_completed']
            accepted=[r for r in records(folder/'solver.jsonl') if r['kind']=='accepted']
            samples=sorted([r for f in (folder/'audits').glob('*/samples.jsonl') for r in records(f) if r['kind']=='sample'],key=lambda r:(r['state_id'],r['sample']))
            if reference is None:reference=(accepted,samples)
            compare(reference[0],accepted);compare(reference[1],samples)
            phases=phase_sums([folder/'timings.jsonl']);audits=phase_sums((folder/'audits').glob('*/timings.jsonl'))
            counters={r['phase']:r['value'] for r in records(folder/'timings.jsonl') if r['event']=='counter'}
            row=dict(version=version,warmup=iteration==0,wall_seconds=m['wall_seconds'],solver_seconds=m['solver_wall_seconds'],checkpoint_bytes=m['checkpoint_bytes'],adaptive_phases=phases,audit_phases=audits,counters=counters)
            report['runs'].append(row);print(json.dumps(row),flush=True)
    report['numerical_comparison']='Marking/training decisions exact; floating observables rtol=1e-8, atol=1e-10; hashes excluded.'
    report['medians']={}
    for version in ['before','after']:
        rows=[r for r in report['runs'] if r['version']==version and not r['warmup']]
        result={key:statistics.median(r[key] for r in rows) for key in ['wall_seconds','solver_seconds','checkpoint_bytes']}
        for group in ['adaptive_phases','audit_phases']:
            keys=set().union(*(r[group] for r in rows))
            result[group]={key:statistics.median(r[group].get(key,0) for r in rows) for key in sorted(keys)}
        report['medians'][version]=result
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
