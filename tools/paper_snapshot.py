"""Render current ALOD runs with existing controls and a separate full-run tail fit."""
import argparse
import copy
import csv
import gzip
import json
import shutil
from pathlib import Path
import numpy as np
from paper.data import read, clean, sha
from paper.figures import render_figures, STYLE_PROVENANCE
from paper.selection import completed_cycles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--historical', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cycles', type=int, default=25)
    args = parser.parse_args()
    current = read(args.native)['runs']
    runs = clean(read(args.historical)['runs'])
    normalized = clean(current)
    for name, run in normalized.items():
        # Preserve explicit cycle numbers: state IDs alone do not define cycles.
        for row, original in zip(run['rows'], current[name]['rows']):
            row['coarse_cycle'] = original['coarse_cycle']
        run['source_directory'] = current[name]['source_directory']
        if name in ('e1-ALOD', 'e2-ALOD'):
            run['mesh_cycle'] = current[name]['mesh_cycle']
            run['mesh_state'] = current[name]['mesh_state']
            if run['mesh_cycle'] != args.cycles:
                raise ValueError('mesh and curve endpoints must use the same cycle')
        runs[name] = run
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    for folder in ('figures', 'data', 'tables', 'tools/paper'):
        (output / folder).mkdir(parents=True)
    metrics = {}
    plotted = copy.deepcopy(runs)
    for name in ('e1-ALOD', 'e2-ALOD'):
        full = runs[name]
        endpoints = completed_cycles(full)
        if endpoints[-1]['coarse_cycle'] != 28:
            raise ValueError('the requested fit requires a completed 28-cycle run')
        tail = [r for r in endpoints if 19 <= r['coarse_cycle'] <= 28]
        if [r['coarse_cycle'] for r in tail] != list(range(19, 29)):
            raise ValueError('missing or duplicated cycle endpoints')
        x = np.log([r['N_on'] for r in tail])
        y = np.log([r['E'] for r in tail])
        slope, intercept = np.polyfit(x, y, 1)
        r2 = 1 - np.sum((y - slope * x - intercept)**2) / np.sum((y-y.mean())**2)
        metrics[name] = dict(p=-float(slope), log_intercept=float(intercept),
                             R_squared=float(r2), cycles=list(range(19, 29)),
                             points=tail, final=endpoints[-1])
        with (output / 'tables' / (name + '-fit-points.csv')).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['coarse_cycle', 'state', 'N_on', 'E'])
            writer.writeheader()
            writer.writerows({k: r[k] for k in writer.fieldnames} for r in tail)
        display = plotted[name]
        display['rows'] = [r for r in full['rows'] if r['coarse_cycle'] <= args.cycles]
        last = completed_cycles(display)[-1]
        if last['coarse_cycle'] != args.cycles or last['state'] != full['mesh_state']:
            raise ValueError('display endpoint does not match the exported mesh')
        display['promotions'] = [r for r in full.get('promotions', []) if r['state_id'] <= last['state']]
        display['audits'] = [r for r in full.get('audits', []) if r['state_id'] <= last['state']]
        metrics[name]['plotted_endpoint'] = last
    issues = []
    render_figures(plotted, output / 'figures', issues)
    if issues:
        raise ValueError('; '.join(issues))
    for filename, data in [('full-runs.json.gz', runs), ('plotted-runs.json.gz', plotted)]:
        (output / 'data' / filename).write_bytes(gzip.compress(
            json.dumps(dict(schema=1, runs=data), separators=(',', ':')).encode(), mtime=0))
    (output / 'convergence.json').write_text(json.dumps(metrics, indent=2)+'\n')
    lines = ['# Convergence fit: cycles 19 through 28', '',
             'Fit log(E) = a - p log(N_on) by unweighted least squares on ten completed cycle endpoints.',
             'E is the nominal member (sample 0) relative energy error; N_on = coarse free DOFs + enrichment rank.',
             'The fit uses the full 28-cycle runs, independently of the 25-cycle figure cutoff.', '',
             '| Experiment | p | R squared | Final N_on | Final relative error |',
             '|---|---:|---:|---:|---:|']
    for name, m in metrics.items():
        lines.append(f'| {name} | {m["p"]:.8f} | {m["R_squared"]:.8f} | {m["final"]["N_on"]} | {m["final"]["E"]:.10g} |')
    (output / 'convergence.md').write_text('\n'.join(lines)+'\n')
    (output / 'README.md').write_text('''# Paper figures from current ALOD experiments

Five PDF figures and PNG previews use the existing tuned plotting modules.
E1/E2 ALOD curves show the initial state and completed cycles 1 through 25.
Their ALOD mesh panels use the same cycle-25 endpoint. AFEM mesh panels and
AFEM/UFEM/LOD controls retain the previously supplied historical results and
original display limits. E3 uses current ALOD runs and historical controls;
ALOD displays completed cycles within the first 1% error crossing prefix.

The separate convergence report uses cycles 19 through 28 of the full current
E1/E2 runs. It is not fitted to the shortened figure curves. Full and plotted
canonical data are saved separately. Historical controls are not new runs of
the optimized executable. No experiments or snapshot files were removed.
''')
    (output / 'reproduce.py').write_text('''"""Regenerate the five figures from the saved, explicitly selected data."""
from pathlib import Path
import sys
root = Path(__file__).resolve().parent
sys.path.insert(0, str(root / 'tools'))
from paper.data import read
from paper.figures import render_figures
issues = []
render_figures(read(root / 'data/plotted-runs.json.gz')['runs'], root / 'figures', issues)
if issues:
    raise RuntimeError('; '.join(issues))
''')
    with (output / 'README.md').open('a') as f:
        f.write('\nRegenerate the figures with `python reproduce.py` (NumPy and Matplotlib required).\n')
    shutil.copy2(__file__, output / 'tools' / Path(__file__).name)
    for p in (Path(__file__).parent / 'paper').glob('*.py'):
        shutil.copy2(p, output / 'tools/paper' / p.name)
    report = dict(cycles=args.cycles, fit_cycles=[19,28], style=STYLE_PROVENANCE,
                  source_hashes=dict(native=sha(args.native), historical=sha(args.historical)),
                  origins={k:r.get('origin') for k,r in runs.items()},
                  files={str(p.relative_to(output)).replace('\\','/'):sha(p)
                         for p in output.rglob('*') if p.is_file()})
    (output / 'manifest.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:{q:v[q] for q in ['p','R_squared','final','plotted_endpoint']} for k,v in metrics.items()},indent=2))


if __name__ == '__main__':
    main()
