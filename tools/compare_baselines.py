"""Compare a complete small baseline trajectory with an old-kernel fixture."""
import argparse
import json
import os
import subprocess
from pathlib import Path
from compare_oracle import compare


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--fixture',type=Path,required=True)
    p.add_argument('--threads',type=int,default=1)
    a=p.parse_args();f=json.loads(a.fixture.read_text(encoding='utf-8'))
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
    r=subprocess.run([str(a.executable.resolve()),*f['args'],f'--threads={a.threads}'],
        text=True,capture_output=True,env=env,timeout=180)
    if r.returncode: raise SystemExit(r.stderr)
    actual=[json.loads(x) for x in r.stdout.splitlines()]
    delta=compare(f['output'],actual)
    print(f'{a.fixture.stem}: {len(actual)} states, threads={a.threads}, max absolute delta={delta:.3g}')


if __name__=='__main__':main()
