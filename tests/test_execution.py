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
