import os, sys, unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from execution import runtime_environment
from run_fixed import DEFAULT,validate

class RuntimePolicy(unittest.TestCase):
    def test_no_forced_blas_cap_or_affinity(self):
        with patch.dict(os.environ,{'OPENBLAS_NUM_THREADS':'12','OMP_NUM_THREADS':'96','OMP_PLACES':'cores','GOMP_CPU_AFFINITY':'0-3'},clear=True):
            env=runtime_environment(0)
            self.assertEqual(env['OPENBLAS_NUM_THREADS'],'12')
            self.assertEqual(env['OMP_NUM_THREADS'],'96')
            self.assertNotIn('OMP_PLACES',env)
            self.assertNotIn('GOMP_CPU_AFFINITY',env)
            self.assertNotIn('MKL_NUM_THREADS',env)

    def test_python_campaign_workers_do_not_inherit_numeric_parent_identity(self):
        with patch.dict(os.environ,{'ALOD_PARENT_PID':'123'}):
            self.assertNotIn('ALOD_PARENT_PID',runtime_environment(supervised=False))
            self.assertEqual(runtime_environment()['ALOD_PARENT_PID'],str(os.getpid()))

    def test_explicit_thread_count_is_user_controlled(self):
        validate(dict(DEFAULT,threads=128))
        self.assertEqual(runtime_environment(128)['OMP_NUM_THREADS'],'128')
        validate(dict(DEFAULT,threads=0))
        with patch.dict(os.environ,{'OMP_THREAD_LIMIT':'16','OMP_DYNAMIC':'TRUE'}):
            env=runtime_environment(32)
            self.assertNotIn('OMP_THREAD_LIMIT',env)
            self.assertEqual(env['OMP_DYNAMIC'],'FALSE')
            self.assertEqual(env['OMP_NUM_THREADS'],'32')

if __name__=='__main__':unittest.main()


class AuditContextTests(unittest.TestCase):
    def test_run_scoped_provenance_and_binary_replacement(self):
        import tempfile
        from pathlib import Path
        from unittest.mock import patch
        from run_audit import AuditExecutionContext
        with tempfile.TemporaryDirectory() as directory:
            exe = Path(directory)/'alod_run'
            exe.write_bytes(b'first')
            with patch('runtime_provenance.provenance', return_value={'build': 'test'}) as probe:
                context = AuditExecutionContext(exe)
                context.verify(exe)
                context.verify(exe)
                self.assertEqual(probe.call_count, 1)
                exe.write_bytes(b'replacement')
                with self.assertRaises(ValueError):
                    context.verify(exe)


class AuditSchedulingTests(unittest.TestCase):
    def test_bounded_active_then_parallel_drain(self):
        import tempfile, threading, time
        from async_audit import AuditQueue
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'checkpoints').mkdir()
            for i in range(4):(root/'checkpoints'/f'state-{i:06d}-phase-0-ell-2.bin').write_bytes(b'x')
            gate=threading.Event();lock=threading.Lock();running=0;peak=0;started=threading.Event()
            def work(queue, checkpoint):
                nonlocal running,peak
                with lock:running+=1;peak=max(peak,running);started.set()
                gate.wait(5)
                with lock:running-=1
                return []
            with patch('run_audit.AuditExecutionContext'),patch.object(AuditQueue,'work',work):
                queue=AuditQueue(root,root/'exe',[0],workers=1,threads=32,drain_workers=2)
                try:
                    queue.discover();self.assertTrue(started.wait(2));queue.discover()
                    self.assertEqual(len(queue.futures),1);self.assertEqual(len(queue.pending),3)
                    queue.begin_drain()
                    deadline=time.monotonic()+2
                    while peak<2 and time.monotonic()<deadline:time.sleep(.01)
                    self.assertEqual(peak,2);self.assertEqual(len(queue.futures),2)
                    gate.set();accepted=[{'state_id':i} for i in range(4)]
                    self.assertEqual(queue.finish(accepted),[])
                    self.assertTrue(all(state['audit']==[] for state in accepted))
                finally:gate.set();queue.close()

    def test_close_does_not_start_pending_snapshots(self):
        import tempfile
        from async_audit import AuditQueue
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'checkpoints').mkdir()
            for i in range(3):(root/'checkpoints'/f'state-{i:06d}-phase-0-ell-2.bin').write_bytes(b'x')
            def work(queue, checkpoint):queue.cancel.wait(5)
            with patch('run_audit.AuditExecutionContext'),patch.object(AuditQueue,'work',work):
                queue=AuditQueue(root,root/'exe',[0]);queue.discover();queue.close()
                self.assertEqual(len(queue.futures),1);self.assertEqual(len(queue.pending),2)
