"""Prepare the requested E1 -> E2 -> E3 barriers; never launch implicitly."""
import argparse
import json
from pathlib import Path
from run_adaptive import validate
from campaign_schedule import validate_dependencies

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--start-at-e2', action='store_true',
                        help='Prepare a fresh E2/E3 campaign without rerunning E1.')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/'configs').mkdir()
    plan = ([] if args.start_at_e2 else [('e1', 'e1_main', [], 160)])
    plan += [('e2', 'e2_main', [] if args.start_at_e2 else ['e1'], 160),
            ('e3_k8', 'e3_k8', ['e2'], 24),
            ('e3_k16', 'e3_k16', ['e2'], 48),
            ('e3_k32', 'e3_k32', ['e3_k8', 'e3_k16'], 128),
            ('e3_k64', 'e3_k64', ['e3_k32'], 160),
            ('e3_k128', 'e3_k128', ['e3_k64'], 200)]
    jobs = []
    for name, preset, dependencies, memory in plan:
        config = json.loads((ROOT/f'configs/adaptive/{preset}.json').read_text())
        config.update(threads=32, checkpoint_interval_cycles=5,
                      ell_absolute_threshold=3.2/config['wavenumber'])
        if name in ('e1', 'e2'):
            config.update(cycles=28, state_limit=0, exact_target=-1.)
        validate(config)
        (args.output/f'configs/{name}.json').write_text(json.dumps(config,indent=2)+'\n')
        jobs.append(dict(name=name,config=f'configs/{name}.json',kind='adaptive',
                         depends_on=dependencies,memory_gib=memory,threads=32,
                         audit_workers=1,audit_threads=32,audit_drain_workers=1))
    validate_dependencies(jobs)
    spec = dict(schema=1,jobs=jobs,maximum_active=2,threads_per_job=32,
                reserve_gib=64,emergency_available_gib=32,disk_reserve_gib=35,
                ell_diagnostic_kappa=[],
                scope='E1/E2: 28 cycles; E3 existing target/horizon; tau=3.2/k; five-cycle restart points. Dependency completion includes audits.')
    (args.output/'campaign.json').write_text(json.dumps(spec,indent=2)+'\n')
    print(json.dumps(dict(prepared=len(jobs),launched=False)))


if __name__ == '__main__':
    main()
