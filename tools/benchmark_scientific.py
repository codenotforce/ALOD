"""Bounded alternating cold/incremental and full/exact workflow measurements."""
import argparse,json,os,statistics
from pathlib import Path
from unittest.mock import patch
from run_adaptive import run
from benchmark_deep_reuse import phases,compare,records
from compare_adaptive import near
from checkpoint_io import digest

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    config=json.loads((Path(__file__).resolve().parents[1]/'configs/adaptive/e2_main.json').read_text())
    config.update(level=8,gap=4,cycles=1,state_limit=3,threads=32,ell=2,maximum_ell=2,ell_mode='fixed',
        member_ids=list(range(16)),training_ids=list(range(16)),radius=10.,extra_checks=[],exact_target=1e-12,exact_scope='nominal')
    evidence=dict(schema=1,config=config,executable_sha256=digest(a.executable),runs=[])
    for repeat in range(3):
        folders={}
        for name in (('cold','incremental','exact') if repeat%2==0 else ('exact','incremental','cold')):
            folder=a.output/f'{repeat}-{name}';folders[name]=folder
            env={'ALOD_INTEGRATION_COLD':'1','ALOD_CHECKPOINT_SYNC':'1'} if name=='cold' else {}
            with patch.dict(os.environ,env):m=run(dict(config,audit_mode='exact' if name=='exact' else 'full'),folder,a.executable,audit_threads=32,audit_drain_workers=2)
            if not m['audit_complete']:raise RuntimeError('incomplete benchmark audit')
            timing,counters=phases(folder/'timings.jsonl');factors=0;reused=0
            for f in (folder/'audits').glob('*/run.json'):
                audit=json.loads(f.read_text());factors+=audit['timing']['reference_factorizations'];reused+=audit['timing'].get('accepted_values_reused',False)
            row=dict(mode=name,repeat=repeat,warmup=repeat==0,wall_seconds=m['wall_seconds'],solver_seconds=m['solver_wall_seconds'],
                drain_seconds=m['audit_drain_seconds'],phases=timing,counters=counters,reference_factorizations=factors,accepted_snapshots_reused=reused)
            evidence['runs'].append(row);print(json.dumps(row),flush=True)
        compare(folders['cold'],folders['incremental']);compare(folders['cold'],folders['exact'])
        def samples(folder):return sorted([r for f in (folder/'audits').glob('*/samples.jsonl') for r in records(f) if r['kind']=='sample'],key=lambda r:(r['state_id'],r['sample']))
        before=samples(folders['cold'])
        for name in ('incremental','exact'):
            after=samples(folders[name]);assert len(before)==len(after)
            for x,y in zip(before,after):near(x['E'],y['E'],relative=1e-8,absolute=1e-10)
    evidence['numerically_equivalent']=True;evidence['medians']={}
    for mode in ('cold','incremental','exact'):
        rows=[x for x in evidence['runs'] if x['mode']==mode and not x['warmup']]
        out={k:statistics.median(x[k] for x in rows) for k in ('wall_seconds','solver_seconds','drain_seconds')}
        out['phases']={k:statistics.median(x['phases'].get(k,0) for x in rows) for k in set().union(*(x['phases'] for x in rows))}
        evidence['medians'][mode]=out
    (a.output/'summary.json').write_text(json.dumps(evidence,indent=2)+'\n')

if __name__=='__main__':main()
