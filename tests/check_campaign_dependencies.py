"""Dependency barriers include successful audit completion through job exit."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from campaign_schedule import ready, validate_dependencies

jobs=[dict(name='e1'),dict(name='e2',depends_on=['e1']),
      dict(name='k8',depends_on=['e2']),dict(name='k16',depends_on=['e2']),
      dict(name='k32',depends_on=['k8','k16']),dict(name='k64',depends_on=['k32']),
      dict(name='k128',depends_on=['k64'])]
validate_dependencies(jobs)
finished=[]
for expected in (['e1'],['e2'],['k8','k16'],['k32'],['k64'],['k128']):
    complete={j['name'] for j in finished}
    assert [j['name'] for j in jobs if j['name'] not in complete and ready(j,finished)]==expected
    finished += [dict(name=name,status='complete') for name in expected]
assert not ready(jobs[1],[dict(name='e1',status='failed')])
try:
    validate_dependencies([dict(name='a',depends_on=['b']),dict(name='b',depends_on=['a'])])
except ValueError:
    pass
else:
    raise AssertionError('cycle accepted')
print('Campaign dependency barriers passed.')
