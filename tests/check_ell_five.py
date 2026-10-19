"""ell=5 must survive publication, audit, checkpoint load and restart."""
import json, subprocess, sys, tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from run_adaptive import DEFAULT, run
from run_audit import audit
from checkpoint_io import inspect
exe=Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory() as td:
    root=Path(td)
    c=dict(DEFAULT,level=3,gap=2,ell=5,maximum_ell=5,ell_mode='fixed',cycles=1,
        member_ids=[0],training_ids=[0],threads=2,audit=True,audit_mode='exact',
        checkpoint_retention='delete_after_audit',keep_final=True)
    paused=run(c,root/'run',exe,pause_state=0,audit_threads=2)
    assert paused['status']=='paused'
    cp,meta=inspect(root/'run',exe)
    assert meta['ell']==5
    resumed=run(c,root/'run',exe,resume='auto',audit_threads=2)
    assert resumed['status']=='complete' and resumed['audit_complete'] and resumed['validation_passed']
    cp,meta=inspect(root/'run',exe)
    assert meta['ell']==5 and meta['state_id']==2
    assert len(list((root/'run/checkpoints').glob('state-*.bin')))==1
    a=audit(cp,root/'audit',exe,member_ids=[0],threads=2,ell_override=5,audit_mode='exact')
    assert a['audit_complete'] and a['validation_passed']
    row=json.loads((root/'audit/samples.jsonl').read_text().splitlines()[0])
    assert row['ell']==5 and row['PG_residual']<1e-8
print('ell=5 checkpoint/audit/restart passed')
