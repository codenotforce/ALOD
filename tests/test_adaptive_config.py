import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from run_adaptive import DEFAULT, validate


class AdaptiveConfigTests(unittest.TestCase):
    def test_all_presets(self):
        for path in (ROOT / "configs/adaptive").glob("*.json"):
            validate(json.loads(path.read_text()))

    def test_reject_bad_contracts(self):
        for changes in [dict(m_ref=3), dict(reference_theta=[.3, 0]), dict(ell_mode="unknown"),
                        dict(method="AFEM"), dict(training_ids=[16]), dict(extra_checks=[True]),
                        dict(rank_cap=25), dict(radius=float("nan")), dict(force_promotions=[0]),
                        dict(maximum_dense_entries=0), dict(member_ids=[1])]:
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                validate(dict(DEFAULT, **changes))

    def test_control_horizons_and_samples(self):
        for method, cycles, m in [("alod", 25, 2), ("afem", 96, 1)]:
            for family in (False, True):
                name=f"e1_{'family' if family else 'nominal'}_{method}.json"
                c=json.loads((ROOT/"configs/adaptive"/name).read_text())
                self.assertEqual((c["cycles"], c["m_ref"], c["ell"], c["ell_mode"]), (cycles,m,3,"fixed"))
                self.assertEqual(c["member_ids"], list(range(48)))
                self.assertEqual(c["training_ids"], list(range(16)) if family else [0])
                self.assertTrue(c["audit"])
