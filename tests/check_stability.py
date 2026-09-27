"""Cross-check the inf-sup diagnostic against dense SVD on both domains."""
import json
import math
import subprocess
import sys

executable = sys.argv[1]
for problem, k in [('E1', 8), ('E2', 16)]:
    values = []
    for seed in (17, 29):
        row = json.loads(subprocess.check_output(
            [executable, problem, str(k), '5', '4', str(seed), '750', '1e-9'],
            text=True))
        assert 0 < row['dense_gamma'] and math.isfinite(row['gamma'])
        assert row['dense_relative_difference'] < 1e-7
        assert row['inverse_residual'] <= 1e-9
        assert row['original_residual'] < 1e-7
        values.append(row['gamma'])
    assert abs(values[0]/values[1]-1) < 1e-7
invalid = subprocess.run([executable, 'E1', '0', '5', '4', '17', '750', '1e-9'],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
assert invalid.returncode != 0
print('Inf-sup SVD and independent-start checks passed.')
