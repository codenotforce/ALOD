"""Compare bounded initial states with the supplied E1/E2 production evidence."""
import json, sys, tempfile, subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from run_adaptive import run
from compare_adaptive import near
exe=Path(sys.argv[1]).resolve()
fixtures=json.loads((ROOT/'tests/fixtures/v4/initial_states.json').read_text())
with tempfile.TemporaryDirectory() as td:
    for problem,expected in fixtures.items():
        config=json.loads((ROOT/f'configs/adaptive/{problem.lower()}_main.json').read_text())
        config.update(cycles=0,state_limit=1,threads=2,audit=True,maximum_nodes=20000)
        folder=Path(td)/problem
        run(config,folder,exe,timeout=240)
        rows=[json.loads(x) for x in (folder/'solver.jsonl').read_text().splitlines()]
        actual={'state':next(x for x in rows if x['kind']=='accepted'),
                'check':next(x for x in rows if x['kind']=='ell_check')}
        for kind,fields in expected.items():
            for key,value in fields.items():
                if isinstance(value,str):assert actual[kind][key]==value
                else:near(actual[kind][key],value,relative=2e-7,absolute=1e-10)
        assert actual['state']['aot_factorizations']==0
        samples=[json.loads(x) for f in (folder/'audits').glob('*/samples.jsonl') for x in f.read_text().splitlines()]
        for row in samples:
            if row['kind']=='sample' and row['accepted_solution_distance'] is not None:
                assert row['accepted_solution_distance']<1e-8
        if problem=='E2':
            latest=folder/'checkpoints'/(folder/'checkpoints/latest').read_text().strip()
            wrong=subprocess.run([str(exe.with_name('alod_audit')),problem,f'--members={folder/"members.txt"}',f'--checkpoint={latest}','--enrichment-tests=adjoint'],capture_output=True,text=True)
            assert wrong.returncode!=0 and 'enrichment policy differs' in wrong.stderr
        print(problem,'initial observables and independent audit passed')
