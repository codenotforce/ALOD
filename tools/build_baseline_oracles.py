"""Explicitly build a probe against a freshly built external legacy source tree.

The probe shares only CLI/refinement scheduling and serialization with alod_run.
All FEM, estimator, marking, manufactured solution and complete LOD model code
comes from the original source/library. Normal builds never call this tool.
"""
import argparse
import hashlib
import json
import os
import shlex
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--legacy-build', type=Path, required=True)
    p.add_argument('--work', type=Path, default=ROOT/'_local/baseline-oracle')
    p.add_argument('--fixtures', type=Path, default=ROOT/'tests/fixtures/baselines')
    p.add_argument('--native',action='store_true',help='Match a legacy library built with -march=native')
    a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True);a.fixtures.mkdir(parents=True,exist_ok=True)
    exe=a.work/'legacy_baseline'
    eigen_flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','eigen3'],text=True))
    # Match the library flags: Eigen allocation alignment is ABI-sensitive.
    subprocess.run(['c++','-O3','-DNDEBUG',*(['-march=native'] if a.native else []),'-std=c++20','-DALOD_LEGACY_ORACLE=1',
        '-I'+str(a.source/'include'),'-I'+str(ROOT/'include'),*eigen_flags,
        str(ROOT/'apps/alod_run.cpp'),str(ROOT/'src/problems/problems.cpp'),
        str(a.legacy_build/'liblod2d_core.a'),'-fopenmp','-lumfpack','-o',str(exe)],check=True)
    # The library is a development input, never copied into ALOD or a fixture.
    sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    provenance={'library_sha256':sha(a.legacy_build/'liblod2d_core.a'),
        'model_sha256':sha(a.source/'src/helmholtz/model.cpp'),
        'estimator_sha256':sha(a.source/'src/helmholtz/adaptive/estimator.cpp'),
        'probe_sha256':sha(ROOT/'apps/alod_run.cpp'),
        'architecture_flags':'-march=native' if a.native else 'portable default architecture, matching ALOD',
        'scope':'Old production kernels and full LOD model; shared deterministic scheduling/output probe. Not a production campaign.'}
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
    for problem in ['E1','E2']:
        for method in ['AFEM','UFEM','SLOD']:
            args=[problem,method,'--initial-level=2','--states=3','--emit-solution=1']
            r=subprocess.run([str(exe.resolve()),*args],capture_output=True,text=True,env=env,check=True,timeout=180)
            data={'provenance':provenance,'args':args,'output':[json.loads(x) for x in r.stdout.splitlines()]}
            (a.fixtures/f'{problem}_{method}.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
            print(f'{problem} {method}: captured {len(data["output"])} old states',flush=True)


if __name__=='__main__': main()
