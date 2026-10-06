import json
import os
import struct
import sys
import tempfile
import unittest
from concurrent.futures import Future
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from checkpoint_io import digest, validate_resume_config
from async_audit import AuditQueue
from run_adaptive import DEFAULT, arguments, validate
from state_retention import StateRetention, accepted_hash


class RetentionTests(unittest.TestCase):
    def test_config_and_resume_contract(self):
        for mode in ('all', 'delete_after_audit', 'periodic'):
            config = dict(DEFAULT, audit=True, checkpoint_retention=mode, checkpoint_interval_cycles=5)
            validate(config)
            validate_resume_config(DEFAULT, config)
            argv = arguments(config, Path('members.txt'))
            self.assertFalse(any('checkpoint-retention' in arg or 'keep-final' in arg for arg in argv))
        for changes in (dict(checkpoint_retention='unknown'), dict(keep_final=1),
                        dict(audit=True, checkpoint_retention='periodic'),
                        dict(checkpoint_retention='delete_after_audit')):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                validate(dict(DEFAULT, **changes))

    def fixture(self, output, sid=1, phase=0, geometry='mesh-0000000000000001.bin'):
        row = dict(kind='accepted', state_id=sid, coarse_cycle=1, cycle_complete=True,
                   wall_seconds=1.0, target_reached=False)
        (output/'solver.jsonl').write_text(json.dumps(row)+'\n')
        name = f'state-{sid:06d}-phase-{phase}-ell-1.bin'
        metadata = json.dumps(dict(journal=json.dumps(row)+'\n', state_id=sid)).encode()
        payload = struct.pack('<QQ', 0x3354504B43444F4C, len(metadata)) + metadata
        payload += struct.pack('<Q', len(geometry)) + geometry.encode()
        for folder in ('checkpoints', 'audit_snapshots'):
            (output/folder/'meshes').mkdir(parents=True, exist_ok=True)
            (output/folder/'meshes'/geometry).write_bytes(b'geometry')
        path = output/'audit_snapshots'/name
        path.write_bytes(payload)
        os.link(path, output/'checkpoints'/name)
        (output/'checkpoints/latest').write_text(name+'\n')
        job = dict(status='complete', accepted_hash=accepted_hash(row))
        future = Future(); future.set_result([])
        return row, name, job, future

    def test_latest_protected_until_supervisor_finishes_then_all_links_removed(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            row, name, job, future = self.fixture(output)
            retention = StateRetention(output, dict(DEFAULT, cycles=1, checkpoint_retention='delete_after_audit'))
            retention.reclaim({name:job}, {name:future})
            self.assertTrue((output/'checkpoints'/name).exists())
            retention.reclaim({name:job}, {name:future}, producer_done=True)
            self.assertFalse((output/'checkpoints/latest').exists())
            self.assertFalse((output/'checkpoints'/name).exists())
            self.assertFalse((output/'audit_snapshots'/name).exists())
            self.assertFalse(list(output.rglob('mesh-*.bin')))

    def test_uncommitted_failed_and_wrong_generation_states_are_never_removed(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            row, name, job, future = self.fixture(output)
            retention = StateRetention(output, dict(DEFAULT, checkpoint_retention='delete_after_audit'))
            (output/'solver.jsonl').write_text(json.dumps(row))  # Partial row, not committed.
            retention.reclaim({name:job}, {name:future}, producer_done=True)
            self.assertTrue((output/'checkpoints'/name).exists())
            (output/'solver.jsonl').write_text(json.dumps(row)+'\n')
            for changed in (dict(job, status='failed'), dict(job, accepted_hash='stale')):
                retention.reclaim({name:changed}, {name:future}, producer_done=True)
                self.assertTrue((output/'checkpoints'/name).exists())

    def test_periodic_and_terminal_selection(self):
        retention = StateRetention(Path('.'), dict(DEFAULT, cycles=6, checkpoint_interval_cycles=5,
                                                  checkpoint_retention='periodic'))
        state = dict(state_id=10, coarse_cycle=5, cycle_complete=True, target_reached=False)
        self.assertTrue(retention.keep(state))
        self.assertFalse(retention.keep(dict(state, state_id=9, cycle_complete=False)))
        self.assertTrue(retention.keep(dict(state, state_id=12, coarse_cycle=6)))
        self.assertTrue(retention.keep(dict(state, state_id=2, coarse_cycle=1, target_reached=True)))
        retention.config['state_limit'] = 4
        self.assertTrue(retention.keep(dict(state, state_id=3, coarse_cycle=2, cycle_complete=False)))
        retention.mode = 'delete_after_audit'
        self.assertFalse(retention.keep(dict(state, state_id=3)))
        retention.config['keep_final'] = True
        self.assertTrue(retention.keep(dict(state, state_id=3)))

    def test_shared_and_quarantined_geometry_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            row, name, job, future = self.fixture(output)
            recovery = output/'checkpoints/recovery'
            recovery.mkdir()
            (recovery/'state-000000-phase-0-ell-1.bin').write_bytes((output/'checkpoints'/name).read_bytes())
            retention = StateRetention(output, dict(DEFAULT, checkpoint_retention='delete_after_audit'))
            retention.reclaim({name:job}, {name:future}, producer_done=True)
            self.assertEqual(len(list(output.rglob('mesh-*.bin'))), 2)
            self.assertTrue((recovery/'state-000000-phase-0-ell-1.bin').exists())

    def test_reclaimed_receipts_validate_results_and_reject_replayed_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            row, name, job, future = self.fixture(output)
            audit_dir = output/'audits/result'
            audit_dir.mkdir(parents=True)
            records = [dict(kind='sample', state_id=row['state_id'], sample=0, PG_residual=0.),
                       dict(kind='audit_complete', samples=1)]
            sample_file = audit_dir/'samples.jsonl'
            sample_file.write_text(''.join(json.dumps(r)+'\n' for r in records))
            checkpoint_hash = digest(output/'checkpoints'/name)
            manifest = dict(audit_complete=True, validation_passed=True, audit_mode='exact',
                            member_ids=[0], checkpoint_sha256=checkpoint_hash,
                            output_sha256=digest(sample_file))
            (audit_dir/'run.json').write_text(json.dumps(manifest))
            job.update(checkpoint=name, output='audits/result', checkpoint_sha256=checkpoint_hash)
            (output/'audit_queue.json').write_text(json.dumps(dict(jobs={name:job})))
            for path in (output/'checkpoints'/name, output/'audit_snapshots'/name):
                path.unlink()
            with patch('run_audit.AuditExecutionContext'):
                def queue():
                    return AuditQueue(output, Path('unused'), [0], audit_mode='exact')
                first = queue()
                try:
                    first.discover()
                    self.assertEqual(first.futures[name].result(), records[:1])
                finally:
                    first.pool.shutdown(wait=True); first.stopped=True
                sample_file.write_text('corrupted results\n')
                second = queue()
                try:
                    with self.assertRaisesRegex(ValueError, 'missing or invalid'):
                        second.discover()
                finally:
                    second.pool.shutdown(wait=True); second.stopped=True
                # Even a valid old receipt cannot authorize a replayed state
                # whose accepted row has changed.
                (output/'solver.jsonl').write_text(json.dumps(dict(row, wall_seconds=2.))+'\n')
                third = queue()
                try:
                    third.discover()
                    self.assertNotIn(name, third.futures)
                finally:
                    third.pool.shutdown(wait=True); third.stopped=True


if __name__ == '__main__':
    unittest.main()
