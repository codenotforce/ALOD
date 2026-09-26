"""Actual asynchronous overlap, absolute lazy checks, and crash recovery."""
import json, os, subprocess, sys, tempfile, threading, time
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from run_adaptive import DEFAULT,run
from run_audit import audit
from checkpoint_io import inspect,validate_resume_config,digest
from execution import run_lease
exe=Path(sys.argv[1]).resolve()
def events(p):return [json.loads(x) for x in (p/'solver.jsonl').read_text().splitlines()]
def equivalent(a,b):
    x,y=events(a),events(b);assert len(x)==len(y)
    for i,j in zip(x,y):
        for key in ('kind','state_id','ell','coarse_fingerprint','reference_fingerprint','coarse_marks','reference_marks','dictionary_hash','kernel_hash'):
            assert i.get(key)==j.get(key),(key,i.get(key),j.get(key))
with tempfile.TemporaryDirectory() as td:
    root=Path(td)
    c=dict(DEFAULT,problem='E2',level=3,gap=2,ell=1,maximum_ell=2,cycles=2,
           radius=10.,rank_cap=3,member_ids=[0,1,2],training_ids=[0,1],threads=2,audit=True,
           ell_ratio_mode='solution_scaled',ell_threshold=1e20,ell_absolute_threshold=1e-20)
    reference=root/'reference';run(dict(c,audit=False),reference,exe)
    startup=root/'startup';startup.mkdir()
    (startup/'run.json').write_text(json.dumps(dict(config=dict(c,audit=False))))
    (startup/'members.txt').write_bytes((reference/'members.txt').read_bytes())
    m=run(dict(c,audit=False),startup,exe,resume='auto')
    assert m['startup_recovery'];equivalent(reference,startup)
    gates=[r for r in events(reference) if r['kind']=='ell_check']
    assert gates[0]['absolute_exceeded'] and not gates[0]['ratio_exceeded'] and gates[0]['action']=='promote'
    assert gates[1]['state_id']==0 and gates[1]['ell']==2
    held=root/'overlap';started=threading.Event();release=threading.Event();errors=[]
    def slow(*args,**kwargs):
        started.set()
        if not release.wait(90):raise RuntimeError('solver waited for audit')
        return audit(*args,**kwargs)
    def launch():
        try:run(c,held,exe)
        except BaseException as e:errors.append(e)
    with patch('run_audit.audit',side_effect=slow):
        worker=threading.Thread(target=launch);worker.start()
        try:
            assert started.wait(90)
            deadline=time.monotonic()+90
            final=held/'checkpoints/state-000004-phase-0-ell-2.bin'
            while not final.exists() and time.monotonic()<deadline:time.sleep(.05)
            assert final.exists(),'adaptive solver blocked on held audit'
            try:
                with run_lease(held):raise AssertionError('double writer lease accepted')
            except RuntimeError:pass
        finally:release.set();worker.join(120)
    assert not worker.is_alive() and not errors,errors
    equivalent(reference,held)
    # Completed audits must not be recomputed when only postprocessing resumes.
    with patch('run_audit.audit',side_effect=AssertionError('completed audit repeated')):
        m=run(c,held,exe,resume='auto')
        assert m['audit_complete']
    # Recover an orphan accepted transaction despite a missing pointer and torn tail.
    resumed=root/'recovered';run(dict(c,audit=False),resumed,exe,pause_state=1)
    cp,metadata=inspect(resumed,exe)
    (resumed/'checkpoints/latest').unlink()
    (resumed/'solver.jsonl').write_text(''.join((resumed/'solver.jsonl').read_text().splitlines(True)[:metadata['committed_lines']])+'{"torn":')
    run(dict(c,audit=False),resumed,exe,resume='auto');equivalent(reference,resumed)
    bad,_=inspect(resumed,exe)
    damaged=bytearray(bad.read_bytes());damaged[-20]^=1;bad.write_bytes(damaged)
    run(dict(c,audit=False),resumed,exe,resume='auto');equivalent(reference,resumed)
    assert list((resumed/'checkpoints/recovery').glob('*.bin'))
    # An interrupted audit is retried under a new directory; math is unchanged.
    failed=root/'audit-failed'
    with patch('run_audit.audit',side_effect=RuntimeError('injected interruption')):
        m=run(c,failed,exe)
    assert m['solver_completed'] and not m['audit_complete']
    run(c,failed,exe,resume='auto');equivalent(reference,failed)
    # Real supervisor termination: its child cannot continue mutating the run.
    killed=root/'killed';cfg=root/'config.json';cfg.write_text(json.dumps(dict(c,audit=False,cycles=30)))
    log=(root/'killed.log').open('w')
    child=subprocess.Popen([sys.executable,str(ROOT/'tools/run_adaptive.py'),'--config',str(cfg),'--output',str(killed),'--executable',str(exe)],stdout=log,stderr=log)
    deadline=time.monotonic()+90
    while not (killed/'checkpoints/latest').exists() and child.poll() is None and time.monotonic()<deadline:time.sleep(.01)
    assert (killed/'checkpoints/latest').exists()
    child.kill();child.wait();log.close();time.sleep(.5)
    before={p.name:p.stat().st_size for p in (killed/'checkpoints').glob('*.bin')};time.sleep(.5)
    assert before=={p.name:p.stat().st_size for p in (killed/'checkpoints').glob('*.bin')},'orphan solver continued after supervisor death'
    # Keep the original declared horizon; pause after recovery for a bounded test.
    cp,meta=inspect(killed,exe)
    target=meta['state_id']+1
    m=run(dict(c,audit=False,cycles=30),killed,exe,resume='auto',pause_state=target)
    assert m['status']=='paused'
    for tau in (.2,-2,float('nan')):
        try:validate_resume_config(c,dict(c,ell_absolute_threshold=tau))
        except ValueError:pass
        else:raise AssertionError('changed tau accepted on resume')
    print(json.dumps(dict(overlap_verified=True,absolute_only_promotion=True,audit_reuse=True,torn_journal_recovered=True,real_supervisor_kill_recovered=True)))
