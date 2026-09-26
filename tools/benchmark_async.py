"""Compare bounded audit drain scheduling with numerical equivalence gates."""
import argparse, json, statistics
from pathlib import Path
from run_adaptive import run
from benchmark_reuse import records, compare, phase_sums
from checkpoint_io import digest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--repeats',type=int,default=2)
    args=parser.parse_args()
    if args.repeats<1:parser.error('repeats must be positive')
    args.output.mkdir(parents=True,exist_ok=False)
    c=json.loads((Path(__file__).resolve().parents[1]/'configs/adaptive/e2_server32.json').read_text())
    c.update(cycles=2,state_limit=5,audit=True)
    report=dict(config=c,binary_sha256=digest(args.executable),runs=[])
    reference=None
    for iteration in range(args.repeats+1):
        for workers in ([1,2] if iteration%2==0 else [2,1]):
            folder=args.output/f'{iteration}-drain{workers}'
            m=run(c,folder,args.executable,timeout=1800,audit_workers=1,audit_threads=32,audit_drain_workers=workers)
            assert m['solver_completed'] and m['audit_complete'],m['status']
            accepted=[r for r in records(folder/'solver.jsonl') if r['kind']=='accepted']
            samples=sorted([r for f in (folder/'audits').glob('*/samples.jsonl') for r in records(f) if r['kind']=='sample'],key=lambda r:(r['state_id'],r['sample']))
            if reference is None:reference=(accepted,samples)
            compare(reference[0],accepted);compare(reference[1],samples)
            checks=[r for r in records(folder/'solver.jsonl') if r['kind']=='ell_check']
            assert checks and max(r['ritz_residual'] for r in checks)<=1e-9
            counters={}
            for r in records(folder/'timings.jsonl'):
                if r['event']=='counter' and r['phase'].startswith('training_schur_'):
                    counters[r['phase']]=counters.get(r['phase'],0)+r['value']
            # Auditors may overlap only after the numerical solver releases its budget.
            jobs=json.loads((folder/'audit_queue.json').read_text())['jobs'].values()
            timeline=sorted([(j['started_at'],1) for j in jobs]+[(j['finished_at'],-1) for j in jobs])
            active=peak=0
            for when,change in timeline:
                active+=change;peak=max(peak,active)
                if when<m['solver_finished_at']-0.01:assert active<=1
            assert peak<=workers
            row=dict(drain_workers=workers,warmup=iteration==0,wall_seconds=m['wall_seconds'],solver_seconds=m['solver_wall_seconds'],drain_seconds=m['audit_drain_seconds'],peak_audits=peak,maximum_ritz_residual=max(r['ritz_residual'] for r in checks),counters=counters,adaptive_phases=phase_sums([folder/'timings.jsonl']))
            report['runs'].append(row);print(json.dumps(row),flush=True)
    report['medians']={str(w):{key:statistics.median(r[key] for r in report['runs'] if r['drain_workers']==w and not r['warmup']) for key in ['wall_seconds','solver_seconds','drain_seconds']} for w in [1,2]}
    report['validation']='Five states, 48 audits per state; marking/training decisions exact, floating observables rtol=1e-8/atol=1e-10; warmup discarded; no full-horizon claim.'
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
