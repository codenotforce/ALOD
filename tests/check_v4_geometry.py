"""Check actual NVB diameters and configurable LOD reference refinement."""
import json, math, subprocess, sys
from pathlib import Path
exe=Path(sys.argv[1]).resolve()
rows=[json.loads(x) for x in subprocess.check_output([str(exe.with_name('alod_geometry')),'--paper-v4'],text=True).splitlines()]
assert [r['wavenumber'] for r in rows]==[8,16,32,64,128]
for r in rows:
    assert math.isclose(r['kH0'],2*math.sqrt(2),rel_tol=1e-12)
    assert math.isclose(r['wavenumber']**3*r['h0']**2,8,rel_tol=1e-12)
    assert r['reference_level']==r['level']+int(math.log2(r['wavenumber']))
def baseline(gap):
    return json.loads(subprocess.check_output([str(exe),'E1','SLOD','--initial-level=3','--states=1',f'--reference-gap={gap}'],text=True))
a,b=baseline(2),baseline(3)
assert a!=b
assert a['residual']<1e-10 and b['residual']<1e-10
for gap in [0,9]:
    assert subprocess.run([str(exe),'E1','SLOD',f'--reference-gap={gap}'],capture_output=True).returncode!=0
print(json.dumps(dict(geometry=rows,baseline_gaps=[2,3],status='passed')))
