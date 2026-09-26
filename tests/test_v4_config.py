"""Mathematical policy changes must never become implicit checkpoint resumes."""
import json, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from checkpoint_io import validate_resume_config
from run_adaptive import DEFAULT, validate
from run_baseline import validate as validate_baseline

class Version4Config(unittest.TestCase):
    def test_paper_presets(self):
        for name in ['e1_main','e2_main']+[f'e3_k{k}' for k in (8,16,32,64,128)]:
            c=json.loads((ROOT/f'configs/adaptive/{name}.json').read_text())
            validate(c)
            self.assertEqual((c['m_ref'],c['ell_threshold'],c['enrichment_tests']),(2,.3,'kernel_lift'))

    def test_policy_changes_rejected(self):
        for key,value in [('enrichment_tests','adjoint'),('ell_ratio_mode','solution_scaled'),('ell_threshold',.3)]:
            with self.assertRaises(ValueError):validate_resume_config(DEFAULT,dict(DEFAULT,**{key:value}))

    def test_legacy_policy_is_explicit_on_resume(self):
        old={k:v for k,v in DEFAULT.items() if k not in ('ell_ratio_mode','ell_threshold','enrichment_tests')}
        validate_resume_config(old,dict(DEFAULT,enrichment_tests='adjoint'))
        with self.assertRaises(ValueError):validate_resume_config(old,DEFAULT)

    def test_baseline_gap_validation(self):
        c=json.loads((ROOT/'configs/smoke/e1_slod.json').read_text())
        self.assertEqual(validate_baseline(c)['reference_gap'],4)
        self.assertEqual(validate_baseline(dict(c,reference_gap=7))['reference_gap'],7)
        for gap in (0,9,2.5,True):
            with self.assertRaises(ValueError):validate_baseline(dict(c,reference_gap=gap))

if __name__=='__main__':unittest.main()
