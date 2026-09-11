import csv
import io
import unittest
from tools.check_repository import CSV_FIELD_LIMIT, path_violations, read_audit
from tools.compare_oracle import compare


class MigrationToolsTest(unittest.TestCase):
    def test_portable_paths_and_urls(self):
        self.assertEqual(path_violations('docs/plan.md <ALOD_ROOT>/build https://github.com/example/repo'), [])

    def test_machine_paths_are_rejected(self):
        drive = chr(67) + ':' + chr(92) + 'private' + chr(92) + 'file.txt'
        unix = '/' + 'home' + '/someone/project'
        unc = chr(92)*2 + 'server' + chr(92) + 'share'
        for path in [drive, unix, unc]:
            self.assertEqual(path_violations(path), [1])

    def test_nonfinite_and_changed_solution_fail(self):
        with self.assertRaises(ValueError):
            compare({'solution': [1.0]}, {'solution': [float('nan')]})
        with self.assertRaises(ValueError):
            compare({'solution': [1.0]}, {'solution': [1.01]})

    def test_missing_and_duplicate_audits_fail(self):
        header = 'state,sample,split,exact_norm,e,E,PG_residual\n'
        row = '0,0,train,1,0.1,0.1,0\n'
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            read_audit(header+row+row,1)
        with self.assertRaisesRegex(ValueError, 'missing'):
            read_audit(header+row,1)

    def test_legacy_long_csv_field_is_not_truncated(self):
        csv.field_size_limit(CSV_FIELD_LIMIT)
        value = '12;'*50000
        parsed = list(csv.DictReader(io.StringIO('element_ids\n'+value+'\n')))
        self.assertEqual(parsed[0]['element_ids'], value)


if __name__ == '__main__':
    unittest.main()
