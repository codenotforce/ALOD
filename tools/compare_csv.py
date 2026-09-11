"""Compare selected numeric and identity columns of two archived CSV files."""
import argparse
import csv
import json
import math
from pathlib import Path


def compare(expected, actual, columns, *, atol=1e-12, rtol=1e-8, audit_only_missing=(), extra_diagnostics=()):
    csv.field_size_limit(64 * 1024 * 1024)
    with Path(expected).open(newline='', encoding='utf-8-sig') as f:
        left = list(csv.DictReader(f))
    with Path(actual).open(newline='', encoding='utf-8-sig') as f:
        right = list(csv.DictReader(f))
    if len(left) != len(right):
        raise ValueError(f'row count: {len(left)} != {len(right)}')
    maximum = 0.0
    extra_count = 0
    identities = {'H_step','N_H','N_reference','sample','seed','ell','rank','step',
                  'marked_element','shared_rank','retrained_rank','regional_budget','inherit'}
    for row, (a, b) in enumerate(zip(left, right)):
        for key in columns:
            x, y = a[key], b[key]
            if key in identities:
                if int(x) != int(y):
                    raise ValueError(f'row {row}, {key}: integer identity differs')
                continue
            try:
                u, v = float(x), float(y)
            except ValueError:
                if x != y:
                    raise ValueError(f'row {row}, {key}: identity differs')
                continue
            if (key in audit_only_missing and a.get('split') != 'train'
                    and b.get('split') == a.get('split') and math.isnan(u) and math.isnan(v)):
                continue
            if key in extra_diagnostics and math.isnan(u) and (math.isnan(v) or math.isfinite(v)):
                extra_count += int(math.isfinite(v))
                continue
            if not math.isfinite(u) or not math.isfinite(v):
                raise ValueError(f'row {row}, {key}: nonfinite number')
            d = abs(u-v)
            maximum = max(maximum, d)
            if d > atol + rtol*abs(u):
                raise ValueError(f'row {row}, {key}: {u} != {v}')
    return {'rows': len(left), 'columns': columns, 'atol': atol, 'rtol': rtol,
            'maximum_absolute_difference': maximum, 'extra_diagnostic_values': extra_count, 'passed': True}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('expected', type=Path)
    p.add_argument('actual', type=Path)
    p.add_argument('--columns', required=True, help='Comma-separated column names')
    a = p.parse_args()
    print(json.dumps(compare(a.expected, a.actual, a.columns.split(',')), indent=2))
