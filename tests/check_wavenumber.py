"""Variable-kappa operators, immutable audits, geometry and restart isolation."""
import json, math, subprocess, sys, tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from run_adaptive import DEFAULT,run
from run_audit import audit
from checkpoint_io import inspect,validate_resume_config
exe=Path(sys.argv[1]).resolve()
geometries=[json.loads(x) for x in subprocess.check_output([str(exe.with_name('alod_geometry'))],text=True).splitlines()]
assert [r['wavenumber'] for r in geometries]==[8,16,32,64]
for r in geometries:
    assert r['kH0']<=math.pi and abs(r['H0_over_h0']-4)<1e-12
with tempfile.TemporaryDirectory() as temp:
    root=Path(temp);identities=set();errors=[]
    for k in [8,16,32,64]:
        c=dict(DEFAULT,wavenumber=k,level=3,gap=2,ell=2,maximum_ell=3,
               cycles=1,state_limit=2,ell_mode='fixed',member_ids=[0],training_ids=[0],threads=2,audit=True)
        folder=root/f'k{k}';m=run(c,folder,exe)
        assert m['solver_completed'] and m['audit_complete']
        cp,meta=inspect(folder,exe);identities.add(meta['space_identity']) if 'space_identity' in meta else None
        baseline=next(r for p in (folder/'audits').rglob('samples.jsonl') for r in map(json.loads,p.read_text().splitlines()) if r['kind']=='sample')
        assert baseline['wavenumber']==k and baseline['PG_residual']<1e-10
        errors.append(baseline['E'])
        boosted=root/f'quad{k}';audit(cp,boosted,exe.with_name('alod_audit'),quadrature_boost=4,refinement_steps=2)
        elevated=root/f'ell{k}';audit(cp,elevated,exe.with_name('alod_audit'),ell_override=3)
        row=json.loads((elevated/'samples.jsonl').read_text().splitlines()[0])
        assert row['ell']==3 and row['source_ell']==2 and row['accepted_solution_distance']>=0
        try:validate_resume_config(c,dict(c,wavenumber=2*k))
        except ValueError:pass
        else:raise AssertionError('changed kappa accepted on resume')
        for method in ['AFEM','UFEM','SLOD']:
            row=json.loads(subprocess.check_output([str(exe),'E1',method,f'--wavenumber={k}',
                '--initial-level=3','--states=1','--threads=2'],text=True))
            assert row['wavenumber']==k and row['residual']<1e-10
    assert max(errors)-min(errors)>.01
    for invalid in ['0','-1','nan','inf']:
        assert subprocess.run([str(exe),'E1','UFEM',f'--wavenumber={invalid}'],capture_output=True).returncode!=0
print(json.dumps(dict(status='passed',geometries=geometries,small_mesh_errors=errors)))
