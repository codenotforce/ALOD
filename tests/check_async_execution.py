"""Compare mathematical trajectories and multi-batch audits across scheduling modes."""
import json
import os
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from run_adaptive import DEFAULT, run
from run_audit import audit
from checkpoint_io import inspect
from compare_adaptive import near

exe = Path(sys.argv[1]).resolve()


def records(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def trajectory(a, b):
    assert len(a) == len(b)
    for x, y in zip(a, b):
        for key in ('kind', 'state_id', 'ell', 'action', 'rank', 'coarse_marks',
                    'reference_marks', 'coarse_fingerprint', 'reference_fingerprint'):
            assert x.get(key) == y.get(key), (key, x.get(key), y.get(key))
        for key in ('theta', 'eta', 'solution', 'targets', 'frozen_scales'):
            if x.get(key) is not None:
                near(x[key], y[key], relative=1e-9, absolute=1e-11)


with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    for problem in ('E1', 'E2'):
        config = dict(DEFAULT, problem=problem, level=4, gap=2, ell=1,
                      maximum_ell=2, ell_threshold=1e-12, ell_ratio_mode="solution_scaled",
                      force_promotions=[0] if problem == "E2" else [], cycles=1, threads=4,
                      member_ids=list(range(8)), training_ids=[0, 1], rank_cap=3,
                      radius=10., emit_solution=True)
        folders = []
        for mode in ('serial', 'async'):
            os.environ['ALOD_ASYNC_DISABLE'] = '1' if mode == 'serial' else '0'
            folder = root / (problem + '-' + mode)
            run(config, folder, exe)
            timing = records(folder / 'timings.jsonl')
            teams = [r for r in timing if r['event'] == 'team']
            assert teams and all(1 <= r['workers'] <= r['requested_workers'] <= 4 for r in teams)
            assert any(r['event'] == 'budget_return' for r in timing) == (mode == 'async')
            folders.append(folder)
        expected = records(folders[0] / 'solver.jsonl')
        trajectory(expected, records(folders[1] / 'solver.jsonl'))
        resumed = root / (problem + '-resume')
        run(config, resumed, exe, pause_state=0, pause_phase='ell')
        checkpoint, _ = inspect(resumed, exe)
        run(config, resumed, exe, resume=checkpoint)
        trajectory(expected, records(resumed / 'solver.jsonl'))
        checkpoint, _ = inspect(folders[1], exe)
        for audit_mode in ('full', 'exact'):
            samples = []
            for mode in ('serial', 'async', 'bounded'):
                os.environ['ALOD_ASYNC_BUFFER_BYTES'] = '0' if mode == 'bounded' else str(1024**3)
                os.environ['ALOD_ASYNC_DISABLE'] = '1' if mode == 'serial' else '0'
                output = root / (problem + '-' + audit_mode + '-' + mode)
                audit(checkpoint, output, exe.with_name('alod_audit'),
                      threads=4, batch_size=3, audit_mode=audit_mode)
                records(output / 'timings.jsonl')  # Concurrent writes must remain valid JSONL.
                rows = records(output / 'samples.jsonl')
                samples.append([r for r in rows if r['kind'] == 'sample'])
                complete = next(r for r in rows if r['kind'] == 'audit_complete')
                assert complete['batch_pipeline'] == (mode == 'async' and not complete['accepted_values_reused'])
            assert all(len(rows) == 8 for rows in samples)
            for a, b in [(a, b) for actual in samples[1:] for a, b in zip(samples[0], actual)]:
                assert a['sample'] == b['sample']
                for key in ('e', 'f', 'g', 'E', 'F', 'G', 'exact_norm', 'energy', 'PG_residual'):
                    if a.get(key) is not None:
                        near(a[key], b[key], relative=1e-10, absolute=1e-12)
        os.environ.pop('ALOD_ASYNC_BUFFER_BYTES', None)
        print(problem, 'serial/async, promotion, resume, full/exact audit passed', flush=True)
