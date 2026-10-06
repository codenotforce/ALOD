"""Real asynchronous audits: reclamation, numerical parity and receipt recovery."""
import json
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from checkpoint_io import inspect
from compare_adaptive import near
from run_adaptive import DEFAULT, run

exe = Path(sys.argv[1]).resolve()


def states(folder):
    return [json.loads(line) for line in (folder/'events.jsonl').read_text().splitlines()
            if json.loads(line)['kind'] == 'accepted']


def parity(actual, baseline):
    a, b = states(actual), states(baseline)
    assert len(a) == len(b)
    for x, y in zip(a, b):
        for field in ('state_id', 'coarse_cycle', 'reference_sweep', 'ell', 'rank',
                      'coarse_fingerprint', 'reference_fingerprint', 'coarse_marks', 'reference_marks'):
            assert x[field] == y[field], (field, x[field], y[field])
        for field in ('eta', 'frozen_scales'):
            near(x[field], y[field], relative=1e-10, absolute=1e-12)
        assert len(x['audit']) == len(y['audit'])
        for p, q in zip(x['audit'], y['audit']):
            for field in ('e', 'E', 'energy', 'exact_norm'):
                near(p[field], q[field], relative=1e-10, absolute=1e-12)


with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    config = dict(DEFAULT, problem='E1', level=2, gap=2, ell=1, maximum_ell=1,
                  ell_mode='fixed', cycles=6, threads=2, member_ids=[0], training_ids=[0],
                  checkpoint_interval_cycles=5, audit=True, audit_mode='exact')
    baseline = root/'all'
    assert run(config, baseline, exe, audit_workers=2, audit_threads=2, audit_drain_workers=2)['audit_complete']
    assert len(list((baseline/'audit_snapshots').glob('state-*.bin'))) == 13
    for mode, keep_final, expected in [('delete_after_audit', False, []),
                                       ('delete_after_audit', True, [12]),
                                       ('periodic', False, [10, 12])]:
        folder = root/(mode+str(keep_final))
        current = dict(config, checkpoint_retention=mode, keep_final=keep_final)
        manifest = run(current, folder, exe, audit_workers=2, audit_threads=2, audit_drain_workers=2)
        assert manifest['audit_complete'] and manifest['validation_passed']
        for name in ('checkpoints', 'audit_snapshots'):
            retained = sorted(int(p.name.split('-')[1]) for p in (folder/name).glob('state-*.bin'))
            assert retained == expected, (mode, name, retained)
        if not expected:
            assert not (folder/'checkpoints/latest').exists()
            assert not list(folder.rglob('mesh-*.bin'))
            assert manifest['checkpoint'] is None
        else:
            assert inspect(folder, exe)[1]['state_id'] == expected[-1]
        parity(folder, baseline)
        assert len(json.loads((folder/'audit_queue.json').read_text())['jobs']) == 13
        if mode == 'periodic' or keep_final:
            with patch('run_audit.audit', side_effect=AssertionError('completed audit was repeated')):
                assert run(current, folder, exe, resume='auto', audit_threads=2, audit_drain_workers=1)['audit_complete']
            parity(folder, baseline)
    # Interval zero also publishes training/ell restart files; successful audits
    # must reclaim those superseded phases, rather than leaving most bytes behind.
    lean = root/'zero-interval'
    assert run(dict(config, cycles=1, checkpoint_interval_cycles=0,
                    checkpoint_retention='delete_after_audit'), lean, exe,
               audit_threads=2, audit_drain_workers=1)['audit_complete']
    assert not list(lean.rglob('state-*.bin')) and not list(lean.rglob('mesh-*.bin'))
    replay = run(dict(config, cycles=1, checkpoint_interval_cycles=0,
                      checkpoint_retention='delete_after_audit'), lean, exe,
                 resume='auto', audit_threads=2, audit_drain_workers=1)
    assert replay['startup_recovery'] and replay['audit_complete']
    assert not list(lean.glob('checkpoints/state-*.bin'))
    # Failed audits retain their input and can be retried. Completed, reclaimed
    # prefix states survive an explicit pause/resume through durable receipts.
    failed = root/'failed'
    with patch('run_audit.audit', side_effect=RuntimeError('injected audit failure')):
        manifest = run(dict(config, cycles=0, checkpoint_retention='delete_after_audit'),
                       failed, exe, audit_threads=2, audit_drain_workers=1)
    assert not manifest['audit_complete'] and list((failed/'checkpoints').glob('state-*.bin'))
    assert run(dict(config, cycles=0, checkpoint_retention='delete_after_audit'), failed, exe,
               resume='auto', audit_threads=2, audit_drain_workers=1)['audit_complete']
    assert not list(failed.rglob('state-*.bin'))
    paused = root/'paused'
    periodic = dict(config, checkpoint_retention='periodic')
    assert run(periodic, paused, exe, pause_state=11, audit_threads=2,
               audit_drain_workers=1)['status'] == 'paused'
    assert inspect(paused, exe)[1]['state_id'] == 11
    assert run(periodic, paused, exe, resume='auto', audit_threads=2,
               audit_drain_workers=1)['audit_complete']
    parity(paused, baseline)
    target = root/'target'
    assert run(dict(periodic, exact_target=1e6), target, exe,
               audit_threads=2, audit_drain_workers=1)['target_reached']
    assert inspect(target, exe)[1]['state_id'] == 0
    capped = root/'partial-cycle'
    assert run(dict(periodic, state_limit=4), capped, exe,
               audit_threads=2, audit_drain_workers=1)['audit_complete']
    assert inspect(capped, exe)[1]['state_id'] == 3
    assert not states(capped)[-1]['cycle_complete']
    # Exercise enrichment, full reference errors and disk-only native workers,
    # in addition to the exact-only E1 shared-memory tests above.
    e2 = dict(config, problem='E2', cycles=1, radius=10., rank_cap=3,
              member_ids=[0, 1], training_ids=[0, 1], audit_mode='full')
    e2_baseline = root/'e2-all'
    assert run(e2, e2_baseline, exe, audit_threads=2, audit_drain_workers=1)['audit_complete']
    from os import environ
    with patch.dict(environ, {'ALOD_AUDIT_DISK_ONLY': '1'}):
        e2_deleted = root/'e2-deleted'
        assert run(dict(e2, checkpoint_retention='delete_after_audit'), e2_deleted, exe,
                   audit_workers=2, audit_threads=2, audit_drain_workers=2)['audit_complete']
    parity(e2_deleted, e2_baseline)
    for a, b in zip(states(e2_deleted), states(e2_baseline)):
        for p, q in zip(a['audit'], b['audit']):
            for field in ('f', 'F', 'g', 'G'):
                near(p[field], q[field], relative=1e-10, absolute=1e-12)
    assert not list(e2_deleted.rglob('state-*.bin'))
print('All/delete/keep-final/periodic retention, multi-worker parity, audit retry, receipt recovery and final target passed.')
