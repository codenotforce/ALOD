"""Accepted basis restoration, independent rebuild and runtime timing coverage."""
import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from run_adaptive import run
from run_audit import audit
from checkpoint_io import inspect


def records(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


exe = Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    for problem in ('E1', 'E2'):
        config = json.loads((ROOT / 'configs/adaptive' / f'{problem.lower()}_main.json').read_text())
        config.update(level=4, gap=2, ell=1, maximum_ell=2, cycles=0, threads=2,
                      member_ids=[0, 1, 2], training_ids=[0, 1], extra_checks=[])
        folder = root / problem
        manifest = run(config, folder, exe, audit_threads=2)
        assert manifest['audit_complete']
        checkpoint, metadata = inspect(folder, exe)
        assert metadata['format'] == 3
        saved = next((folder / 'audits').glob('*/samples.jsonl'))
        cached = records(saved)
        assert cached[-1]['lod_basis_reused'] and cached[-1]['lod_patch_rebuilds'] == 0
        timing = records(saved.parent / 'timings.jsonl')
        assert any(r['phase'] == 'lod_basis_restore' for r in timing)
        assert any(r['phase'] == 'lod_reduced_restore' for r in timing)
        assert not any(r['phase'] == 'lod_reduced_assembly' for r in timing)
        assert not any(r['phase'] == 'lod_local_correctors' for r in timing)
        assert [r['workers'] for r in timing if r['event'] == 'runtime'] == [2]
        independent = root / (problem + '-rebuild')
        audit(checkpoint, independent, exe, threads=2, reuse_basis=False)
        rebuilt = records(independent / 'samples.jsonl')
        assert not rebuilt[-1]['lod_basis_reused'] and rebuilt[-1]['lod_patch_rebuilds'] == 1
        assert [r for r in cached if r['kind'] == 'sample'] == [r for r in rebuilt if r['kind'] == 'sample']
        phases = {r['phase'] for r in records(folder / 'timings.jsonl')}
        assert {'checkpoint_accepted', 'load', 'lod_local_correctors', 'lod_reduced_factor',
                'riesz_prepare', 'riesz_estimate', 'theta', 'strong_residual'} <= phases
print('Restored and rebuilt audit samples match exactly; runtime phase traces validated')
