"""Measure fine-grid inf-sup constants and anchor tau(16)=0.2.

Uses the E3 initial NVB levels, two random starts, and two extra NVB
levels (halving h) for mesh sensitivity. Results are numerical estimates,
not certified inf-sup lower bounds. Existing successful jobs are resumed.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess


LEVELS = {8: 7, 16: 10, 32: 13, 64: 16, 128: 19}


def summarize(output, problem):
    rows = []
    for path in sorted(output.glob(f'{problem}-k*-l*-s*.json')):
        row = json.loads(path.read_text())
        rows.append(row)
    anchors = {r['level'] - LEVELS[16]: r['gamma'] for r in rows
               if r['k'] == 16 and r['seed'] == 17}
    for row in rows:
        offset = row['level'] - LEVELS[int(row['k'])]
        row['mesh_offset'] = offset
        if offset in anchors:
            row['tau_anchor_16'] = .2 * row['gamma'] / anchors[offset]
    fits = []
    for offset in (0, 2):
        selected = sorted((r for r in rows if r['mesh_offset'] == offset
                           and r['seed'] == 17), key=lambda r: r['k'])
        if len(selected) < 2:
            continue
        xs = [math.log(r['k']) for r in selected]
        ys = [math.log(r['gamma']) for r in selected]
        xm, ym = sum(xs)/len(xs), sum(ys)/len(ys)
        slope = sum((x-xm)*(y-ym) for x,y in zip(xs,ys))/sum((x-xm)**2 for x in xs)
        fits.append(dict(mesh_offset=offset, points=len(selected), power=-slope,
                         prefactor=math.exp(ym-slope*xm)))
        if offset in anchors:
            ax = [math.log(r['k']/16) for r in selected]
            ay = [math.log(r['gamma']/anchors[offset]) for r in selected]
            power = -sum(x*y for x,y in zip(ax,ay))/sum(x*x for x in ax)
            fits[-1]['anchored_power'] = power
            fits[-1]['anchored_max_relative_error'] = max(
                abs(anchors[offset]*(16/r['k'])**power/r['gamma']-1)
                for r in selected)
    checks = []
    for k in LEVELS:
        by_run = {(r['mesh_offset'], r['seed']): r for r in rows if r['k'] == k}
        base = by_run.get((0,17))
        if base is None:
            continue
        check = dict(k=k)
        if (0,29) in by_run:
            check['independent_start_relative_difference'] = abs(by_run[0,29]['gamma']/base['gamma']-1)
        if (2,17) in by_run:
            check['refinement_relative_change'] = by_run[2,17]['gamma']/base['gamma']-1
        checks.append(check)
    payload = dict(problem=problem, anchor_k=16, anchor_tau=.2, rows=rows,
                   fits=fits, checks=checks, expected_jobs=15, completed_jobs=len(rows),
                   interpretation='Numerical estimates on uniform meshes; not certified lower bounds or adaptive-state guarantees.')
    temp = output / 'summary.tmp'
    temp.write_text(json.dumps(payload, indent=2)+'\n')
    temp.replace(output / 'summary.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--problem', choices=['E1', 'E2'], default='E1')
    parser.add_argument('--timeout', type=float, default=86400)
    args = parser.parse_args()
    executable = args.executable.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    # Obtain the anchor early; finish all initial meshes before refinement checks.
    jobs = [(k, offset, seed) for offset, seed in ((0,17),(0,29),(2,17))
            for k in (16,8,32,64,128)]
    for k, offset, seed in jobs:
        level = LEVELS[k]+offset
        name = f'{args.problem}-k{k}-l{level}-s{seed}'
        result = args.output / f'{name}.json'
        if result.exists():
            continue
        command = [str(executable), args.problem, str(k), str(level), '4',
                   str(seed), '750', '1e-9']
        print(json.dumps(dict(event='start', job=name)), flush=True)
        with (args.output / f'{name}.log').open('w') as log:
            run = subprocess.run(command, stdout=subprocess.PIPE, stderr=log,
                                 text=True, timeout=args.timeout)
            if run.returncode:
                log.write(run.stdout)
        if run.returncode:
            raise RuntimeError(f'{name} failed; see its log')
        row = json.loads(run.stdout)
        if not math.isfinite(row['gamma']) or row['gamma'] <= 0:
            raise RuntimeError(f'{name}: invalid gamma')
        temp = result.with_suffix('.tmp')
        temp.write_text(json.dumps(row, indent=2)+'\n')
        temp.replace(result)
        summarize(args.output, args.problem)
        print(json.dumps(dict(event='complete', job=name, gamma=row['gamma'],
                              wall_seconds=row['wall_seconds'])), flush=True)
    summarize(args.output, args.problem)


if __name__ == '__main__':
    main()
