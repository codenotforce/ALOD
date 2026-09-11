import copy
import unittest
from tools.run_fixed import DEFAULT, validate


class FixedConfiguration(unittest.TestCase):
    def test_presets(self):
        for problem in ("E1", "E2"):
            rows = validate(dict(DEFAULT, problem=problem))
            self.assertEqual([r["sample"] for r in rows], list(range(16)))

    def test_explicit_selection(self):
        config = dict(DEFAULT, member_ids=[0, 17], training_ids=[0])
        self.assertEqual(len(validate(config)), 2)
        for ids in ([17], [0, 0], [], [100], [True]):
            with self.assertRaises(ValueError):
                validate(dict(config, training_ids=ids))

    def test_strict_schema(self):
        for changes in (dict(unknown=True), dict(theta=float("nan")), dict(threads=True),
                        dict(graded=1), dict(interpolation="implicit"), dict(riesz_patches="automatic")):
            with self.assertRaises(ValueError):
                validate(dict(copy.deepcopy(DEFAULT), **changes))
