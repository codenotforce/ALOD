"""Shared immutable snapshots, disk spill, persistent fallback and numerical parity."""
import json, os, sys, tempfile, threading
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from run_adaptive import run,DEFAULT
from run_audit import audit
from native_audit import NativeAuditWorker
from compare_adaptive import near
exe=Path(sys.argv[1]).resolve()

def samples(folder):
    return sorted([json.loads(line) for file in (folder/'audits').glob('*/samples.jsonl')
        for line in file.read_text().splitlines() if json.loads(line)['kind']=='sample'],
        key=lambda row:(row['state_id'],row['sample']))

with tempfile.TemporaryDirectory() as directory:
    root=Path(directory)
    c=dict(DEFAULT,problem='E2',level=4,gap=2,ell=1,maximum_ell=1,ell_mode='fixed',cycles=2,
           radius=10.,rank_cap=4,member_ids=list(range(12)),training_ids=list(range(4)),threads=4,audit=True)
    with patch.dict(os.environ,{'ALOD_AUDIT_ONESHOT':'1'}):
        before=run(c,root/'oneshot',exe,audit_threads=4)
    after=run(c,root/'shared',exe,audit_threads=4)
    with patch.dict(os.environ,{'ALOD_SHARED_SNAPSHOT_BYTES':'0'}):
        spilled=run(c,root/'spill',exe,audit_threads=4)
    assert all(m['audit_complete'] for m in (before,after,spilled))
    assert all(json.loads(f.read_text())['timing']['batch_pipeline']
               for f in (root/'shared/audits').glob('*/run.json'))
    expected=samples(root/'oneshot')
    for name in ('shared','spill'):
        actual=samples(root/name);assert len(actual)==len(expected)
        for a,b in zip(expected,actual):
            for key in ('e','f','g','E','F','G','energy','exact_norm'):
                near(a[key],b[key],relative=1e-8,absolute=1e-10)
    manifests=[json.loads(f.read_text()) for f in (root/'shared/audits').glob('*/run.json')]
    if sys.platform=='linux':
        assert any(m.get('memory_snapshot') for m in manifests),'no adaptive memory snapshot consumed'
        assert len({m['native_worker_pid'] for m in manifests})==1,'audits did not share the adaptive process'
    assert all(not json.loads(f.read_text()).get('memory_snapshot') for f in (root/'spill/audits').glob('*/run.json'))
    # A resident standalone worker serves two independent audits and reuses the
    # immutable hierarchy; neither job may inherit mutable solution state.
    checkpoint=sorted((root/'shared/checkpoints').glob('state-*-phase-0-*.bin'))[-1]
    worker=NativeAuditWorker(exe,2)
    try:
        a=audit(checkpoint,root/'repeat-a',exe,threads=2,native_worker=worker)
        b=audit(checkpoint,root/'repeat-b',exe,threads=2,native_worker=worker)
        assert a['native_worker_pid']==b['native_worker_pid'] and a['audit_complete'] and b['audit_complete']
        timings=[json.loads(line) for line in (root/'repeat-b/timings.jsonl').read_text().splitlines()]
        assert any(t.get('phase')=='audit_hierarchy_hits' and t.get('value')==1 for t in timings)
        stopped=threading.Event();stopped.set()
        try:worker.prepare(checkpoint,stopped)
        except RuntimeError:pass
        else:raise AssertionError('cancelled worker request continued')
        assert worker.process.poll() is not None
    finally:worker.close()
    print(json.dumps(dict(shared_snapshots=sum(m.get('memory_snapshot',False) for m in manifests),
        disk_spill_validated=True,persistent_hierarchy_reuse=True,cancellation_verified=True)))
