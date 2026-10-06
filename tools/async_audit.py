"""Durable checkpoint-backed audit queue; workers never mutate solver state."""
import json, threading, time, os
from concurrent.futures import Future, ThreadPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path
from checkpoint_io import atomic_text, digest, sync_directory
from state_retention import StateRetention, accepted_hash, state_header

class AuditQueue:
    def __init__(self, output, executable, members, workers=1, threads=None, *, drain_workers=None, shared_endpoint=None, audit_mode="full", config=None):
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
        self.shared_endpoint=shared_endpoint;self.audit_mode=audit_mode
        self.persistent=not os.environ.get('ALOD_AUDIT_ONESHOT') and self.executable.stem=='alod_run'
        self.retention=StateRetention(output,config)
        self.previous={}
        queue_file=self.output/'audit_queue.json'
        if queue_file.exists():
            self.previous=json.loads(queue_file.read_text()).get('jobs',{})

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
        atomic_text(self.output/'audit_queue.json',json.dumps(dict(schema=2,jobs={**self.previous,**self.jobs}),indent=2)+'\n')

    def completed_samples(self, job):
        """Validate a durable receipt, even when its state input was reclaimed."""
        directory=(self.output/job['output']).resolve()
        if directory.parent != (self.output/'audits').resolve():
            raise ValueError('audit receipt points outside the run audit directory')
        manifest=json.loads((directory/'run.json').read_text())
        if (not manifest.get('audit_complete') or not manifest.get('validation_passed')
                or manifest.get('audit_mode','full')!=self.audit_mode
                or manifest['checkpoint_sha256']!=job['checkpoint_sha256']
                or manifest['member_ids']!=self.members
                or manifest['output_sha256']!=digest(directory/'samples.jsonl')):
            raise ValueError('completed audit receipt failed validation')
        records=[json.loads(line) for line in (directory/'samples.jsonl').read_text().splitlines()]
        samples=[row for row in records if row['kind']=='sample']
        sid=int(job['checkpoint'].split('-')[1])
        if (len(samples)!=len(self.members) or {r['sample'] for r in samples}!=set(self.members)
                or any(r['state_id']!=sid or not 0<=r['PG_residual']<=1e-8 for r in samples)
                or not records or records[-1].get('kind')!='audit_complete'
                or records[-1]['samples']!=len(samples)):
            raise ValueError('completed audit sample/state coverage mismatch')
        return samples

    def restore_completed(self):
        self.retention.update_journal()
        for name,job in list(self.previous.items()):
            if name in self.jobs or job.get('status')!='complete' or not job.get('accepted_hash'):
                continue
            if job.get('checkpoint')!=name:
                raise ValueError('audit receipt checkpoint identity mismatch')
            state=self.retention.accepted.get(int(name.split('-')[1]))
            if state is None or state['accepted_hash']!=job['accepted_hash']:
                continue  # Rollback/replay cannot reuse a superseded state's receipt.
            existing=next((p for p in self.retention.paths(name) if p.exists()),None)
            try:
                if existing is not None and digest(existing)!=job['checkpoint_sha256']:
                    continue
                samples=self.completed_samples(job)
            except (OSError,ValueError,KeyError):
                if existing is None:
                    raise ValueError(f'cannot recover reclaimed state {name}: audit results are missing or invalid')
                continue  # Retained input permits an ordinary audit retry.
            future=Future();future.set_result(samples)
            with self.lock:
                self.jobs[name]=dict(job);self.futures[name]=future;self.persist()
            del self.previous[name]

    def discover(self):
        self.restore_completed()
        # Periodic restart points and per-state audit inputs have separate
        # lifetimes. Prefer the audit path so shared-memory lookup still hits.
        snapshots={p.name:p for folder in ('checkpoints','audit_snapshots')
                   for p in (self.output/folder).glob('state-*-phase-0-ell-*.bin')}
        for checkpoint in sorted(snapshots.values()):
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
            metadata,_=state_header(checkpoint)
            accepted=[json.loads(line) for line in metadata['journal'].splitlines()
                      if json.loads(line).get('kind')=='accepted']
            if len(accepted)!=1 or accepted[0]['state_id']!=metadata['state_id']:
                raise ValueError('audit input has no unique accepted-state record')
            row_hash=accepted_hash(accepted[0])
            sha=digest(checkpoint);base=self.output/'audits'/(sha[:16]+'-'+self.audit_mode)
            # A killed worker's directory is immutable evidence; retries get new names.
            candidates=[base]+sorted(base.parent.glob(base.name+'-retry-*'))
            found=None
            for candidate in candidates:
                try:
                    m=json.loads((candidate/'run.json').read_text())
                    if m.get('audit_mode','full')==self.audit_mode and m.get('audit_complete') and m['checkpoint_sha256']==sha and m['member_ids']==self.members and m['output_sha256']==digest(candidate/'samples.jsonl'):
                        found=candidate;break
                except (OSError,ValueError,KeyError):pass
            if found is None:
                found=base;attempt=0
                while found.exists():
                    attempt+=1;found=base.with_name(base.name+f'-retry-{attempt}')
                audit(checkpoint,found,self.executable,batch_size=8,threads=self.threads,cancel=self.cancel,execution_context=self.execution_context,native_worker=self.native_worker(),audit_mode=self.audit_mode)
            samples=[json.loads(line) for line in (found/'samples.jsonl').read_text().splitlines() if json.loads(line)['kind']=='sample']
            state=int(name.split('-')[1])
            if len(samples)!=len(self.members) or {r['sample'] for r in samples}!=set(self.members) or any(r['state_id']!=state for r in samples):
                raise ValueError('audit sample/state coverage mismatch')
            if self.retention.mode!='all':
                # A receipt must not outlive results buffered only in the page
                # cache when its input is about to be removed.
                for file in (found/'samples.jsonl',found/'run.json'):
                    with file.open('rb') as stream:os.fsync(stream.fileno())
                sync_directory(found)
            with self.lock:
                self.jobs[name].update(status='complete',checkpoint_sha256=sha,accepted_hash=row_hash,output=str(found.relative_to(self.output)),finished_at=time.time());self.persist()
            return samples
        except BaseException as e:
            with self.lock:
                self.jobs[name].update(status='failed',reason=str(e),finished_at=time.time());self.persist()
            raise

    def reclaim(self, *, producer_done=False, protected=()):
        # Filesystem mutation stays on the supervisor thread. Native workers
        # release input ownership before their futures become complete.
        with self.lock:
            self.retention.reclaim(self.jobs,self.futures,producer_done=producer_done,protected=protected)

    def retention_summary(self):
        return self.retention.summary()

    def finish(self, accepted):
        self.discover();self.begin_drain()
        while self.pending or any(not f.done() for f in self.futures.values()):
            self.pump()
            self.reclaim()
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
