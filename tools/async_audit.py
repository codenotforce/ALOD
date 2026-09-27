"""Durable checkpoint-backed audit queue; workers never mutate solver state."""
import json, threading, time, os
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path
from checkpoint_io import atomic_text, digest

class AuditQueue:
    def __init__(self, output, executable, members, workers=1, threads=None, *, drain_workers=None, shared_endpoint=None):
        if type(workers) is not int or workers<1:raise ValueError('audit workers must be positive')
        drain_workers=workers if drain_workers is None else drain_workers
        if type(drain_workers) is not int or drain_workers<workers:raise ValueError("invalid audit drain workers")
        self.limit=workers;self.drain_workers=drain_workers;self.pending={}
        self.output=Path(output);self.executable=Path(executable);self.members=members
        from run_audit import AuditExecutionContext
        self.execution_context=AuditExecutionContext(executable)
        self.threads=threads;self.cancel=threading.Event();self.jobs={};self.futures={}
        self.pool=ThreadPoolExecutor(max_workers=drain_workers,thread_name_prefix='alod-audit')
        self.lock=threading.Lock();self.stopped=False
        self.local=threading.local();self.native_workers=[]
        self.shared_endpoint=shared_endpoint
        self.persistent=not os.environ.get('ALOD_AUDIT_ONESHOT') and self.executable.stem=='alod_run'

    def native_worker(self):
        if not self.persistent:return None
        worker=getattr(self.local,'worker',None)
        from native_audit import SharedAuditWorker
        shared=self.shared_endpoint is not None and Path(self.shared_endpoint).exists()
        if shared and worker is not None and not isinstance(worker,SharedAuditWorker):
            worker.close();worker=None
        if worker is None or worker.process.poll() is not None:
            from native_audit import NativeAuditWorker
            worker=SharedAuditWorker(self.shared_endpoint) if shared else NativeAuditWorker(self.executable,self.threads or 0)
            self.local.worker=worker
            with self.lock:self.native_workers.append(worker)
        return worker

    def close_workers(self):
        for worker in self.native_workers:worker.close()
        self.native_workers.clear()

    def persist(self):
        atomic_text(self.output/'audit_queue.json',json.dumps(dict(schema=1,jobs=self.jobs),indent=2)+'\n')

    def discover(self):
        for checkpoint in sorted((self.output/'checkpoints').glob('state-*-phase-0-ell-*.bin')):
            name=checkpoint.name
            if name in self.jobs:continue
            with self.lock:
                self.jobs[name]=dict(status='pending',checkpoint=name,queued_at=time.time());self.persist()
            self.pending[name]=checkpoint
        self.pump()

    def pump(self):
        # Only running snapshots enter the executor: queued geometry remains on disk.
        available=self.limit-sum(not f.done() for f in self.futures.values())
        for name in list(self.pending)[:max(0,available)]:
            self.futures[name]=self.pool.submit(self.work,self.pending.pop(name))

    def begin_drain(self):
        self.limit=self.drain_workers
        self.pump()

    def work(self, checkpoint):
        from run_audit import audit
        name=checkpoint.name
        with self.lock:
            self.jobs[name].update(status='running',started_at=time.time());self.persist()
        try:
            sha=digest(checkpoint);base=self.output/'audits'/(sha[:16]+'-shared')
            # A killed worker's directory is immutable evidence; retries get new names.
            candidates=[base]+sorted(base.parent.glob(base.name+'-retry-*'))
            found=None
            for candidate in candidates:
                try:
                    m=json.loads((candidate/'run.json').read_text())
                    if m.get('audit_complete') and m['checkpoint_sha256']==sha and m['member_ids']==self.members and m['output_sha256']==digest(candidate/'samples.jsonl'):
                        found=candidate;break
                except (OSError,ValueError,KeyError):pass
            if found is None:
                found=base;attempt=0
                while found.exists():
                    attempt+=1;found=base.with_name(base.name+f'-retry-{attempt}')
                audit(checkpoint,found,self.executable,batch_size=8,threads=self.threads,cancel=self.cancel,execution_context=self.execution_context,native_worker=self.native_worker())
            samples=[json.loads(line) for line in (found/'samples.jsonl').read_text().splitlines() if json.loads(line)['kind']=='sample']
            state=int(name.split('-')[1])
            if len(samples)!=len(self.members) or {r['sample'] for r in samples}!=set(self.members) or any(r['state_id']!=state for r in samples):
                raise ValueError('audit sample/state coverage mismatch')
            with self.lock:
                self.jobs[name].update(status='complete',checkpoint_sha256=sha,output=str(found.relative_to(self.output)),finished_at=time.time());self.persist()
            return samples
        except BaseException as e:
            with self.lock:
                self.jobs[name].update(status='failed',reason=str(e),finished_at=time.time());self.persist()
            raise

    def finish(self, accepted):
        self.discover();self.begin_drain()
        while self.pending or any(not f.done() for f in self.futures.values()):
            self.pump()
            active=[f for f in self.futures.values() if not f.done()]
            if active:wait(active,timeout=.1,return_when=FIRST_COMPLETED)
        self.pool.shutdown(wait=True);self.close_workers();self.stopped=True
        failures=[]
        for state in accepted:
            candidates=[key for key in self.futures if key.startswith(f"state-{state['state_id']:06d}-phase-0-")]
            try:
                if len(candidates)!=1:raise ValueError('missing or ambiguous accepted checkpoint')
                state['audit']=self.futures[candidates[0]].result()
            except Exception as e:failures.append(dict(state_id=state['state_id'],reason=str(e)))
        return failures

    def close(self):
        if not self.stopped:
            self.cancel.set();self.pool.shutdown(wait=True,cancel_futures=True);self.close_workers();self.stopped=True
