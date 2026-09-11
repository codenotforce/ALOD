"""Path-free runtime and numerical build provenance for portable result bundles."""
import hashlib
import json
import os
import platform
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]

def provenance(executable):
    result=dict(platform=platform.system(),machine=platform.machine(),python=platform.python_version(),
                cpu_affinity=sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else None,
                blas_threads=1,omp_proc_bind='false',omp_places=None)
    build=subprocess.run([str(Path(executable).resolve()),'build-info'],capture_output=True,text=True,timeout=30)
    if build.returncode==0:
        result['build']=json.loads(build.stdout)
    git=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,capture_output=True,text=True)
    result['commit']=git.stdout.strip() if git.returncode==0 else None
    if git.returncode==0:
        status=subprocess.run(['git','status','--porcelain'],cwd=ROOT,capture_output=True,text=True,check=True)
        result['dirty']=bool(status.stdout.strip())
    h=hashlib.sha256()
    for directory in ('apps','include','src','tools','configs','data'):
        for path in sorted((ROOT/directory).rglob('*')):
            if path.is_file() and path.suffix in ('.cpp','.hpp','.h','.py','.json'):
                h.update(path.relative_to(ROOT).as_posix().encode()+b'\0')
                h.update(path.read_bytes().replace(b'\r\n',b'\n'))
    h.update((ROOT/'CMakeLists.txt').read_bytes().replace(b'\r\n',b'\n'))
    result['source_tree_sha256']=h.hexdigest()
    return result
