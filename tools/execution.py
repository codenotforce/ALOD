"""Unpinned runtime configuration and an OS-released exclusive run lease."""
import contextlib, os, subprocess, time
from pathlib import Path

def runtime_environment(threads=0, *, supervised=True):
    env=dict(os.environ)
    env.pop("ALOD_PARENT_PID",None)
    if supervised:env["ALOD_PARENT_PID"]=str(os.getpid())
    if threads:
        env['OMP_NUM_THREADS']=str(threads)
        env['OMP_DYNAMIC']='FALSE'
        env.pop('OMP_THREAD_LIMIT',None)
    for key in ('OMP_PLACES','OMP_PROC_BIND','GOMP_CPU_AFFINITY','KMP_AFFINITY'):
        env.pop(key,None)
    # BLAS and OpenMP runtime defaults/user environment are not forced to one.
    return env

@contextlib.contextmanager
def run_lease(output):
    output=Path(output).resolve()
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.with_name(output.name+'.run.lock').open('a+b') as lock:
        lock.seek(0);lock.write(b'0');lock.flush();lock.seek(0)
        try:
            if os.name=='nt':
                import msvcrt
                msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
            else:
                import fcntl
                fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except OSError as e:raise RuntimeError('run directory already has an active supervisor') from e
        try:yield
        finally:
            if os.name=='nt':
                lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1)
            else:fcntl.flock(lock,fcntl.LOCK_UN)

def cancellable_run(command, *, cancel=None, timeout=None, **kwargs):
    start=time.monotonic()
    with subprocess.Popen(command,**kwargs) as process:
        try:
            while process.poll() is None:
                if cancel is not None and cancel.is_set():raise RuntimeError('audit interrupted; checkpoint remains pending')
                if timeout is not None and time.monotonic()-start>timeout:raise subprocess.TimeoutExpired(command,timeout)
                time.sleep(.1)
            return subprocess.CompletedProcess(command,process.returncode)
        except BaseException:
            process.terminate()
            try:process.wait(timeout=5)
            except subprocess.TimeoutExpired:process.kill();process.wait()
            raise
