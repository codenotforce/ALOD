import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from run_baseline import validate,command
from compare_csv import compare
from tempfile import TemporaryDirectory


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


if __name__=='__main__':unittest.main()
