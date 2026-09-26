"""An unpinned campaign worker can run and resume without duplicate work."""
import json, subprocess, sys, tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from run_adaptive import DEFAULT
from checkpoint_io import digest
exe=Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory() as td:
    root=Path(td);(root/'configs').mkdir()
    c=dict(DEFAULT,problem='E2',level=3,gap=2,ell=1,maximum_ell=1,ell_mode='fixed',
           cycles=0,radius=10.,rank_cap=3,member_ids=[0,1],training_ids=[0,1],threads=2,audit=True)
    (root/'configs/job.json').write_text(json.dumps(c))
    spec=dict(schema=1,jobs=[dict(name='job',config='configs/job.json',kind='adaptive',memory_gib=2)],
              maximum_active=1,reserve_gib=0,emergency_available_gib=0,disk_reserve_gib=0,ell_diagnostic_kappa=[])
    (root/'campaign.json').write_text(json.dumps(spec))
    cmd=[sys.executable,str(ROOT/'tools/run_campaign.py'),'--campaign',str(root),'--build',str(exe.parent)]
    subprocess.run(cmd,check=True,timeout=90)
    status=json.loads((root/'status.json').read_text());assert status['status']=='complete'
    assert status['finished'][0]['cpu_binding'] is None
    assert 'taskset' not in status['finished'][0]['command']
    checkpoint=root/'runs/job/checkpoints';before={p.name:digest(p) for p in checkpoint.glob('*.bin')}
    subprocess.run(cmd+['--resume'],check=True,timeout=30)
    assert before=={p.name:digest(p) for p in checkpoint.glob('*.bin')}
    assert len(json.loads((root/'status.json').read_text())['finished'])==1
print('Campaign launch, unpinned execution and completed-job resume passed')
