import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from paper_delivery import adaptive,clean

class PaperDeliveryTests(unittest.TestCase):
    def fixture(self):
        return dict(config=dict(member_ids=[0,1]),states=[dict(state_id=0,coarse_free=8,rank=2,eta=[.3])],
            audits=[dict(state_id=0,sample=i,E=.1,exact_norm=2,role='train') for i in (0,1)])

    def test_missing_and_duplicate_audits(self):
        x=self.fixture();x['audits'].pop()
        with self.assertRaises(ValueError):adaptive(x)
        x=self.fixture();x['audits'].append(x['audits'][0])
        with self.assertRaises(ValueError):adaptive(x)

    def test_no_training_audit_confusion(self):
        r=adaptive(self.fixture());self.assertNotIn('test_count',r['rows'][0])
        self.assertEqual(r['rows'][0]['N_on'],10)

    def test_portable_config_whitelist(self):
        r=adaptive(self.fixture());r['config']['host_workspace']='private-host-location'
        result=clean({'e1-ALOD':r})
        self.assertNotIn('host_workspace',result['e1-ALOD']['config'])
        r['rows'][0]['E']=float('nan')
        with self.assertRaises(ValueError):clean({'e1-ALOD':r})

if __name__=='__main__':unittest.main()
