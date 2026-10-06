"""Exercise preset execution, artifact status and preflight errors end to end."""
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

root=Path(__file__).resolve().parents[1]
exe=Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory() as temp:
    out=Path(temp)/'run'
    cmd=[sys.executable,str(root/'tools/run_baseline.py'),'--config',str(root/'configs/smoke/e1_afem.json'),
         '--executable',str(exe),'--output',str(out)]
    env={**os.environ,'OPENBLAS_NUM_THREADS':'2','MKL_NUM_THREADS':'3'}
    subprocess.run(cmd,check=True,capture_output=True,text=True,env=env)
    info=json.loads((out/'run.json').read_text());states=[json.loads(x) for x in (out/'states.jsonl').read_text().splitlines()]
    assert info['status']=='completed' and info['states_written']==3 and info['stop_reason']=='state_limit'
    assert info['BLAS_threads'] is None
    assert info['thread_environment']['OPENBLAS_NUM_THREADS']=='2'
    assert info['thread_environment']['MKL_NUM_THREADS']=='3'
    assert len(states)==3 and all('solution' not in x for x in states)
    assert subprocess.run(cmd,capture_output=True).returncode!=0
    for args in [['E1','ALOD'],['E1','AFEM','--target=nan'],['E1','AFEM','--states=1','--states=2'],
                 ['E2','SLOD','--initial-level=6','--maximum-nodes=50']]:
        result=subprocess.run([str(exe),*args],capture_output=True,text=True)
        assert result.returncode!=0 and result.stdout==''
    result=subprocess.run([str(exe),'E1','AFEM','--states=4','--target=10'],check=True,capture_output=True,text=True)
    stopped=[json.loads(x) for x in result.stdout.splitlines()]
    assert len(stopped)==1 and stopped[0]['target_reached'] and stopped[0]['stop_reason']=='exact_error_target'
print('Preset artifacts, preflight errors and exact stopping passed.')
