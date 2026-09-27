"""Physical 16-load E2 training where changing local ordering changes POD rank."""
import json,os,sys,tempfile
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from run_adaptive import run
from compare_adaptive import near
exe=Path(sys.argv[1]).resolve()
c=json.loads((ROOT/'configs/adaptive/e2_main.json').read_text())
c.update(level=8,gap=4,ell=2,maximum_ell=2,ell_mode='fixed',cycles=0,state_limit=1,
         radius=10.,member_ids=list(range(16)),training_ids=list(range(16)),threads=4,audit=False,extra_checks=[])
with tempfile.TemporaryDirectory() as directory:
    states=[]
    for name,env in [('direct',{'ALOD_REFERENCE_EXECUTION':'1'}),('coordinates',{})]:
        folder=Path(directory)/name
        with patch.dict(os.environ,env):run(c,folder,exe)
        states.append(next(json.loads(line) for line in (folder/'solver.jsonl').read_text().splitlines() if json.loads(line)['kind']=='accepted'))
    a,b=states
    assert a['rank']==b['rank']==6 and a['working_rank']==b['working_rank']==12
    for key in ('coarse_marks','reference_marks','training_stop','training_evaluations','compression_trials'):assert a[key]==b[key],key
    for key in ('eta','targets','frozen_scales','coarse_mean_total','reference_mean_total'):near(a[key],b[key],relative=1e-8,absolute=1e-10)
    print('Physical E2: rank 12 -> 6, training decisions and estimators preserved.')
