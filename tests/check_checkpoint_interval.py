"""Periodic restarts preserve all-state audits and replay superseded tails."""
import json
import shutil
import sys
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from run_adaptive import DEFAULT, run
from checkpoint_io import inspect

exe=Path(sys.argv[1]).resolve()


def states(folder):
    return [json.loads(line) for line in (folder/'solver.jsonl').read_text().splitlines()
            if json.loads(line)['kind']=='accepted']


def equivalent(a,b):
    x,y=states(a),states(b)
    assert len(x)==len(y)
    for p,q in zip(x,y):
        for key in ('state_id','coarse_cycle','reference_sweep','ell',
                    'coarse_fingerprint','reference_fingerprint','coarse_marks','reference_marks'):
            assert p.get(key)==q.get(key),(key,p.get(key),q.get(key))


with tempfile.TemporaryDirectory() as td:
    root=Path(td)
    config=dict(DEFAULT,problem='E1',level=2,gap=2,ell=1,maximum_ell=1,
                ell_mode='fixed',cycles=6,threads=2,member_ids=[0],training_ids=[0],
                checkpoint_interval_cycles=5,audit=True,audit_mode='exact')
    full=root/'full'; result=run(config,full,exe,audit_threads=2,audit_drain_workers=1)
    assert result['audit_complete'] and result['accepted_states']==13
    points=sorted((full/'checkpoints').glob('state-*.bin'))
    assert [int(p.name.split('-')[1]) for p in points]==[10,12]
    assert all('-phase-0-' in p.name for p in points)
    assert len(list((full/'audit_snapshots').glob('state-*.bin')))==13
    for point in points:
        _,meta=inspect(point,exe)
        assert meta['phase']==0
    # Simulate an interruption after a committed checkpoint and later audits.
    replay=root/'replay';shutil.copytree(full,replay)
    (replay/'checkpoints'/points[-1].name).rename(replay/'saved-final.bin')
    result=run(config,replay,exe,resume='auto',audit_threads=2,audit_drain_workers=1)
    assert result['resumed_from']['state_id']==10 and result['audit_complete']
    assert len(list((replay/'audit_snapshots/recovery').glob('*.bin')))==2
    equivalent(full,replay)
    # Before the first scheduled point, restart from the beginning, preserving
    # the superseded journal and audit files rather than treating them as restarts.
    startup=root/'startup';shutil.copytree(full,startup)
    for point in (startup/'checkpoints').glob('state-*.bin'):
        point.rename(startup/('saved-'+point.name))
    result=run(config,startup,exe,resume='auto',audit_threads=2,audit_drain_workers=1)
    assert result['startup_recovery'] and result['audit_complete']
    assert list((startup/'recovery').glob('solver-before-first-checkpoint-*.jsonl'))
    equivalent(full,startup)
    # Without audit, only scheduled/final geometry and transactions are written.
    lean=root/'lean';run(dict(config,audit=False),lean,exe)
    assert not (lean/'audit_snapshots').exists()
    assert len(list((lean/'checkpoints').glob('state-*.bin')))==2
    assert len(list((lean/'checkpoints/meshes').glob('*.bin')))==2
    equivalent(full,lean)
    # A target reached before cycle five still creates a final restart point.
    target=root/'target';result=run(dict(config,exact_target=1e6),target,exe)
    assert result['target_reached'] and result['accepted_states']==1
    assert inspect(target,exe)[1]['state_id']==0
print('Five-cycle checkpoints, all-state audits, rollback, startup replay and final target passed.')
