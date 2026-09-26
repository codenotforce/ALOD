"""Checkpoint recovery, immutable audits and numerical trajectory equivalence."""
import json
from pathlib import Path
import shutil
import sys
import tempfile
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from run_adaptive import DEFAULT, run
from run_audit import audit
from checkpoint_io import digest, inspect, recover_journal, validate_resume_config
from compare_adaptive import near

exe=Path(sys.argv[1]).resolve()
audit_exe=exe.with_name('alod_audit')
maximum_resume_difference=0.
def difference(a,b):
    return max((difference(x,y) for x,y in zip(a,b)),default=0.) if isinstance(b,list) else abs(a-b)
with tempfile.TemporaryDirectory() as temp:
    root=Path(temp)
    config=dict(DEFAULT, problem='E2', level=3, gap=2, ell=1, maximum_ell=2,
                cycles=1, m_ref=2, reference_theta=[.1,.1], theta=.1,
                radius=10., rank_cap=3, member_ids=list(range(48)), training_ids=[0,1],
                force_promotions=[0], emit_solution=True)
    baseline=root/'baseline'
    run(config,baseline,exe)
    def events(folder):
        return [json.loads(line) for line in (folder/'solver.jsonl').read_text().splitlines()]
    assert all(row['rank']>0 for row in events(baseline) if row['kind']=='accepted')
    def equivalent(folder):
        global maximum_resume_difference
        actual,expected=events(folder),events(baseline)
        assert len(actual)==len(expected)
        for a,b in zip(actual,expected):
            for key in ('kind','state_id','ell','coarse_cycle','reference_sweep','mesh_revision','next_phase',
                        'coarse_fingerprint','reference_fingerprint','rank','coarse_marks','reference_marks','ell_last_check'):
                assert a.get(key)==b.get(key),(key,a.get(key),b.get(key))
            for key in ('solution','eta','targets','frozen_scales','theta'):
                if a.get(key) is not None:
                    near(a[key],b[key],relative=1e-8,absolute=1e-10)
                    maximum_resume_difference=max(maximum_resume_difference,difference(a[key],b[key]))
    for phase,state in [('accepted',1),('training',1),('ell',0)]:
        folder=root/phase
        manifest=run(config,folder,exe,pause_state=state,pause_phase=phase)
        assert manifest['status']=='paused' and not manifest['solver_completed']
        checkpoint,metadata=inspect(folder,exe)
        if phase=='accepted':
            # Simulate interruption after checkpoint rename but before log flush.
            lines=(folder/'solver.jsonl').read_text().splitlines(True)
            (folder/'solver.jsonl').write_text(''.join(lines[:metadata['committed_lines']])+'{"partial":')
        run(config,folder,exe,resume=checkpoint)
        equivalent(folder)
    checkpoint,metadata=inspect(baseline,exe)
    original=digest(checkpoint)
    changed=dict(config,theta=.2)
    try: validate_resume_config(config,changed)
    except ValueError: pass
    else: raise AssertionError('changed math accepted')
    validate_resume_config(config,dict(config,cycles=2,threads=2,audit=True))
    damaged=root/'damaged.bin'
    data=bytearray(checkpoint.read_bytes());data[-20]^=1;damaged.write_bytes(data)
    try: inspect(damaged,exe)
    except ValueError: pass
    else: raise AssertionError('corrupt checkpoint accepted')
    journal=root/'corrupt.jsonl'
    journal.write_text((baseline/'solver.jsonl').read_text().replace('"state_id":0','"state_id":9',1))
    try: recover_journal(journal,metadata)
    except ValueError: pass
    else: raise AssertionError('corrupt committed prefix accepted')
    samples=[]
    for batch in (1,8):
        folder=root/f'audit-{batch}'
        manifest=audit(checkpoint,folder,audit_exe,batch_size=batch,pure=True)
        assert manifest['audit_complete'] and manifest['timing']['reference_factorizations']==1
        samples.append([json.loads(line) for line in (folder/'samples.jsonl').read_text().splitlines() if json.loads(line)['kind']=='sample'])
    assert len(samples[0])==50
    for a,b in zip(*samples):
        for key in ('e','f','g','E','F','G','exact_norm','energy'):
            near(a[key],b[key],relative=1e-10)
    assert digest(checkpoint)==original
    # Repository coefficient edits cannot change a checkpoint's audit inputs.
    import run_audit
    previous=run_audit.validate
    def changed_table(settings):
        rows=previous(settings)
        rows[0]=dict(rows[0],wave_x=99.)
        return rows
    with patch('run_audit.validate',side_effect=changed_table):
        frozen=root/'frozen-audit'
        audit(checkpoint,frozen,audit_exe,member_ids=[0])
    row=json.loads((frozen/'samples.jsonl').read_text().splitlines()[0])
    near(row['e'],samples[0][0]['e'],relative=1e-12)
    # Audit enabled and disabled must leave the mathematical journal unchanged.
    audited=root/'audited'
    manifest=run(dict(config,audit=True),audited,exe)
    assert manifest['solver_completed'] and manifest['audit_complete'] and not manifest['paper_complete']
    equivalent(audited)
    # A postprocessor failure leaves the numerical result available for retry.
    with patch('run_audit.audit',side_effect=RuntimeError('injected audit failure')):
        manifest=run(dict(config,audit=True,cycles=0),root/'audit-failed',exe)
    assert manifest['solver_completed'] and manifest['validation_passed'] and not manifest['audit_complete']
    with patch('export_results.export',side_effect=RuntimeError('injected export failure')):
        manifest=run(dict(config,cycles=0),root/'export-failed',exe)
    assert manifest['solver_completed'] and manifest['validation_passed'] and manifest['status']=='postprocess_failed'
    # A real accepted state 16 is transferred alone; the audit has no H0..15 files.
    h16=dict(config,cycles=16,m_ref=1,reference_theta=[.03],theta=.03,force_promotions=[],
             ell_mode='fixed',maximum_ell=1,member_ids=[0,1],state_limit=17,emit_solution=False)
    h16_manifest=run(h16,root/'h16',exe)
    source,metadata=inspect(root/'h16',exe)
    assert metadata['state_id']==16
    from checkpoint_io import pack
    portable=pack(source,root/'portable',exe)
    isolated=root/'isolated';isolated.mkdir();shutil.copyfile(portable,isolated/'state16.bin')
    report=audit(isolated/'state16.bin',isolated/'fresh',audit_exe,member_ids=[1],fresh=True,batch_size=1)
    assert report['audit_complete'] and report['timing']['state_id']==16
    print(json.dumps(dict(test='P5 checkpoint acceptance',pause_boundaries=3,maximum_resume_absolute_difference=maximum_resume_difference,
                          config=config,h16_config=h16,accepted_ranks=[row['rank'] for row in events(baseline) if row['kind']=='accepted'],
                          audited_members=50,batch_sizes=[1,8],direct_fresh_state=16,h16_checkpoint_bytes=source.stat().st_size,
                          h16_all_checkpoint_bytes=h16_manifest['checkpoint_bytes'],fresh_audit_seconds=report['wall_seconds'])))
