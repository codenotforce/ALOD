"""One bounded, cancellable native worker per active audit thread.

The immutable checkpoint is decoded once by prepare and ownership is transferred
inside the worker to the numerical audit. Queued states remain on disk.
"""
import json
import queue
import struct
import subprocess
import threading
import time
from pathlib import Path
from execution import runtime_environment

class NativeAuditWorker:
    def __init__(self, executable, threads=0):
        self.process=subprocess.Popen([str(Path(executable).resolve()), 'audit-worker'],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            env=runtime_environment(threads))
        self.responses=queue.Queue()
        def read():
            try:
                for line in self.process.stdout:
                    self.responses.put(json.loads(line))
            except BaseException as error:self.responses.put(error)
            finally:self.responses.put(RuntimeError('native audit worker exited'))
        self.reader=threading.Thread(target=read,daemon=True);self.reader.start()
        try:
            if self.receive(None,30).get('ready')!=1:raise RuntimeError('invalid audit worker handshake')
        except BaseException:self.close();raise

    def receive(self,cancel,timeout):
        deadline=time.monotonic()+timeout
        while True:
            if cancel is not None and cancel.is_set():raise RuntimeError('audit interrupted; checkpoint remains pending')
            remaining=deadline-time.monotonic()
            if remaining<=0:raise TimeoutError('native audit worker timed out')
            try:result=self.responses.get(timeout=min(.1,remaining))
            except queue.Empty:continue
            if isinstance(result,BaseException):raise result
            if result.get('ok') is False:raise RuntimeError(result.get('error','native audit failed'))
            return result

    def request(self,fields,cancel=None,timeout=604800):
        try:
            if cancel is not None and cancel.is_set():raise RuntimeError('audit interrupted; checkpoint remains pending')
            for field in fields:
                encoded=str(field).encode('utf-8')
                if len(encoded)>1024*1024:raise ValueError('audit worker field exceeds limit')
                self.process.stdin.write(struct.pack('<I',len(encoded))+encoded)
            self.process.stdin.flush()
            return self.receive(cancel,timeout)
        except BaseException:self.close();raise

    def prepare(self,checkpoint,cancel=None,timeout=604800):
        from checkpoint_io import resolve_checkpoint
        checkpoint=resolve_checkpoint(checkpoint).resolve()
        result=self.request(['prepare',checkpoint],cancel,timeout)
        self.memory_snapshot=result.get('memory_snapshot',False)
        return checkpoint,json.loads(result['metadata'])

    def run(self,arguments,output,cancel=None,timeout=604800):
        output=Path(output).resolve()
        result=self.request(['audit',output/'samples.jsonl',output/'stderr.log',
            output/'timings.jsonl',len(arguments),*arguments],cancel,timeout)
        return subprocess.CompletedProcess(arguments,result['returncode'])

    def close(self):
        if self.process.poll() is None:
            self.process.terminate()
            try:self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:self.process.kill();self.process.wait()
        for stream in (self.process.stdin,self.process.stdout):
            if stream is not None:stream.close()
        self.reader.join(timeout=2)


class SharedAuditWorker(NativeAuditWorker):
    """Connection to a bounded auditor inside the adaptive numerical process."""
    def __init__(self, endpoint):
        import socket
        from types import SimpleNamespace
        self.socket=socket.socket(socket.AF_UNIX,socket.SOCK_STREAM)
        deadline=time.monotonic()+5
        while True:
            try:self.socket.connect(str(endpoint));break
            except (FileNotFoundError,ConnectionRefusedError):
                if time.monotonic()>=deadline:self.socket.close();raise
                time.sleep(.02)
        self.responses=queue.Queue();self.closed=False
        self.process=SimpleNamespace(stdin=self.socket.makefile('wb'),stdout=self.socket.makefile('rb'),
                                     poll=lambda: 0 if self.closed else None,pid=None)
        def read():
            try:
                for line in self.process.stdout:self.responses.put(json.loads(line))
            except BaseException as error:self.responses.put(error)
            finally:self.responses.put(RuntimeError('shared audit connection closed'))
        self.reader=threading.Thread(target=read,daemon=True);self.reader.start()
        try:
            response=self.receive(None,30)
            if response.get('ready')!=1:raise RuntimeError('invalid shared audit handshake')
            self.process.pid=response['pid']
        except BaseException:self.close();raise

    def close(self):
        import socket
        if self.closed:return
        self.closed=True
        try:self.socket.shutdown(socket.SHUT_RDWR)
        except OSError:pass
        self.reader.join(timeout=2)
        self.process.stdin.close();self.process.stdout.close();self.socket.close()


def shutdown_shared(endpoint):
    worker=SharedAuditWorker(endpoint)
    try:worker.request(['shutdown'],timeout=30)
    finally:worker.close()
