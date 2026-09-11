"""Long legacy fields and canonical marking-set integrity."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from export_results import CSV_LIMIT, export, legacy_rows

class ResultIO(unittest.TestCase):
    def test_long_field_and_explicit_limit(self):
        with tempfile.TemporaryDirectory() as temp:
            file=Path(temp)/'legacy.csv'
            field='1;'*100000
            file.write_text('element_ids\n'+field+'\n')
            self.assertEqual(list(legacy_rows(file))[0]['element_ids'],field)
            with file.open('w') as stream:
                stream.write('element_ids\n')
                for _ in range(CSV_LIMIT//1048576+1):
                    stream.write('1'*1048576)
                stream.write('\n')
            with self.assertRaisesRegex(ValueError,'64 MiB'):
                list(legacy_rows(file))

    def test_external_indices(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            marks=list(range(40000))
            event=dict(kind='accepted',state_id=0,coarse_free=4,reference_free=20,ell=2,rank=1,
                       coarse_marks=marks,reference_marks=[3],worst_member_id=0,reference_worst_member_id=1,training_ids=[0,1],frozen_scales=[1.,1.])
            (root/'events.jsonl').write_text(json.dumps(event)+'\n')
            export(root)
            for row in legacy_rows(root/'marking.csv'):
                file=root/'indices'/(row['set_id']+'.json')
                self.assertEqual(hashlib.sha256(file.read_bytes()).hexdigest(),row['sha256'])
                self.assertEqual(len(json.loads(file.read_text())['element_ids']),int(row['count']))
                self.assertEqual(int(row['worst_member_id']),0 if row['mesh']=='coarse' else 1)
            self.assertLess((root/'marking.csv').stat().st_size,4096)
            file.write_text('{}')
            with self.assertRaisesRegex(ValueError,'checksum'):
                export(root)
