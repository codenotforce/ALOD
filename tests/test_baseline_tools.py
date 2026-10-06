import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from run_baseline import validate,command,main
from compare_csv import compare
from tempfile import TemporaryDirectory
from unittest.mock import patch
from types import SimpleNamespace
import json
import os


class BaselineTools(unittest.TestCase):
    def config(self):return dict(schema_version=1,problem='E1',method='AFEM',initial_level=2,states=3,theta=.15,target=0,maximum_nodes=20000,threads=1,emit_solution=False)
    def test_invalid_presets(self):
        for key,value in [('problem','E3'),('method','ALOD'),('states',True),('target',float('nan')),('theta',0),('threads',-1),('extra',1),('emit_solution',1)]:
            c=self.config();c[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):validate(c)
    def test_explicit_arguments(self):
        args=command(self.config());self.assertEqual(args[:2],['E1','AFEM']);self.assertIn('--theta=0.15',args)
        self.assertIn('--emit-solution=0',args)
    def test_csv_rejects_missing_rows_and_bad_training_values(self):
        with TemporaryDirectory() as d:
            a,b=Path(d)/'a.csv',Path(d)/'b.csv'
            a.write_text('split,x\ntrain,1\n');b.write_text('split,x\ntrain,nan\n')
            with self.assertRaises(ValueError):compare(a,b,['split','x'],audit_only_missing=['x'])
            b.write_text('split,x\n')
            with self.assertRaises(ValueError):compare(a,b,['split','x'])

    def test_run_records_child_thread_environment_without_inventing_blas_count(self):
        for settings in ({}, {'OPENBLAS_NUM_THREADS':'12','MKL_NUM_THREADS':'4',
                              'OMP_NUM_THREADS':'96','OMP_THREAD_LIMIT':'8'}):
            with self.subTest(settings=settings),TemporaryDirectory() as temp:
                root=Path(temp);config=root/'config.json';exe=root/'solver';output=root/'run'
                c=self.config();c['threads']=32
                config.write_text(json.dumps(c),encoding='utf-8');exe.write_bytes(b'fixture executable')
                observed={}
                def child(*args,**kwargs):
                    observed.update(kwargs['env'])
                    kwargs['stdout'].write(json.dumps({'stop_reason':'state_limit'})+'\n')
                    return SimpleNamespace(returncode=0)
                argv=['run_baseline.py','--config',str(config),'--executable',str(exe),'--output',str(output)]
                with patch.dict(os.environ,settings,clear=True),patch.object(sys,'argv',argv),\
                     patch('run_baseline.subprocess.run',side_effect=child):
                    main()
                record=json.loads((output/'run.json').read_text())
                self.assertIsNone(record['BLAS_threads'])
                self.assertEqual(record['status'],'completed')
                self.assertEqual(record['thread_environment']['OMP_NUM_THREADS'],'32')
                self.assertIsNone(record['thread_environment']['OMP_THREAD_LIMIT'])
                for key,value in record['thread_environment'].items():
                    self.assertEqual(value,observed.get(key))
                for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
                    self.assertEqual(record['thread_environment'][key],settings.get(key))


if __name__=='__main__':unittest.main()
