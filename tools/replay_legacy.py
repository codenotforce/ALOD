"""Rebuild an explicitly supplied archived production snapshot and replay P0.

Run on Linux/WSL with the old build dependencies installed. No server launcher
is invoked. Legacy source and generated binaries stay in the chosen work area.
"""
import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path
from compare_csv import compare

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--problem',choices=['E1','E2'],required=True)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True)
    a=p.parse_args();key=a.problem.lower();proof=ROOT/'docs/provenance/p0'
    patch=json.loads((proof/f'{key}-resource-patch.json').read_text())
    source=a.source.resolve();work=a.work.resolve()
    if work==source or source in work.parents or work in source.parents:
        raise ValueError('work and original source must be separate directory trees')
    work.mkdir(parents=True,exist_ok=False)
    def ignore(directory,names):
        if Path(directory)==source:
            return [n for n in names if n=='experiments' or n=='.git' or n.startswith('build')]
        return []
    shutil.copytree(source,work/'source',ignore=ignore)
    baseline=Path('experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE.sha256')
    (work/'source'/baseline).parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(source/baseline,work/'source'/baseline)
    runner=work/'source/benchmarks/reference_epoch_paper_runner.cpp'
    digest=hashlib.sha256(runner.read_bytes()).hexdigest()
    if digest!=patch['runner_original_sha256']:raise ValueError('runner does not match the frozen production snapshot')
    content=runner.read_text()
    for before,after in patch['resource_only_changes']:
        if before not in content:raise ValueError('resource patch target is missing')
        content=content.replace(before,after)
    runner.write_text(content)
    shutil.copyfile(proof/f'{key}-replay-config.json',work/'config.json')
    subprocess.run(['cmake','-S',str(work/'source'),'-B',str(work/'build'),'-DCMAKE_BUILD_TYPE=Release',
                    '-DLOD2D_BUILD_TESTS=OFF','-DLOD2D_BUILD_BENCHMARKS=ON','-DCMAKE_EXPORT_COMPILE_COMMANDS=ON'],check=True)
    subprocess.run(['cmake','--build',str(work/'build'),'--target','bench_helmholtz_adaptive_paper','-j2'],check=True)
    env={**os.environ,'OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','OMP_PROC_BIND':'false'}
    exe=work/'build/benchmarks/bench_helmholtz_adaptive_paper'
    with (work/'replay.log').open('w') as log:
        subprocess.run([str(exe),'--config='+str(work/'config.json'),'--output-dir='+str(work/'output'),
                        '--reference-cache-dir='+str(work/'cache'),'--manuscript-baseline='+str(work/'source'/baseline)],
                       env=env,stdout=log,stderr=subprocess.STDOUT,timeout=360,check=True)
    run=next((work/'output').iterdir())
    if key=='e1':
        report=json.loads((proof/'e1-replay-comparison.json').read_text())
        reports={'family':compare(proof/'e1-sc_lod_rhs_family.csv',run/'sc_lod_rhs_family.csv',report['family']['columns']),
                 'marks':compare(proof/'e1-coarse_mark_recipe.csv',run/'coarse_mark_recipe.csv',['step','marked_element'])}
    else:
        report=json.loads((proof/'e2-replay-comparison.json').read_text())
        reports=compare(proof/'e2-accepted-prefix.csv',run/'accepted_control_samples.csv',report['columns'],audit_only_missing=['eta_H','eta_D'],extra_diagnostics=report['extra_diagnostic_columns'])
    (work/'comparison.json').write_text(json.dumps(reports,indent=2)+'\n')
    print(json.dumps(reports,indent=2))


if __name__=='__main__':main()
