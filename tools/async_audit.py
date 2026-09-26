"""Durable checkpoint-backed audit queue; workers never mutate solver state."""
import json, threading, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from checkpoint_io import atomic_text, digest

class AuditQueue:
    def __init__(self, output, executable, members, workers=1, threads=None):
        if type(workers) is not int or workers<1:raise ValueError('audit workers must be positive')
        self.output=Path(output);self.executable=Path(executable);self.members=members
        from run_audit import AuditExecutionContext
        self.execution_context=AuditExecutionContext(executable)
        self.threads=threads;self.cancel=threading.Event();self.jobs={};self.futures={}
        self.pool=ThreadPoolExecutor(max_workers=workers,thread_name_prefix='alod-audit')
        self.lock=threading.Lock();self.stopped=False

    def persist(self):
        atomic_text(self.output/'audit_queue.json',json.dumps(dict(schema=1,jobs=self.jobs),indent=2)+'\n')

    def discover(self):
        for checkpoint in sorted((self.output/'checkpoints').glob('state-*-phase-0-ell-*.bin')):
            name=checkpoint.name
            if name in self.futures:continue
            with self.lock:
                self.jobs[name]=dict(status='pending',checkpoint=name,queued_at=time.time());self.persist()
            self.futures[name]=self.pool.submit(self.work,checkpoint)

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
                audit(checkpoint,found,self.executable,batch_size=8,threads=self.threads,cancel=self.cancel,execution_context=self.execution_context)
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
        self.discover();self.pool.shutdown(wait=True);self.stopped=True
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
            self.cancel.set();self.pool.shutdown(wait=True,cancel_futures=True);self.stopped=True
