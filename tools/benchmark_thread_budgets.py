"""Bounded, sequential scheduling comparison with numerical trajectory checks."""
import argparse
import json
import os
from pathlib import Path
import resource
import statistics
import subprocess
import time

from compare_adaptive import near
from execution import runtime_environment
from run_adaptive import DEFAULT, arguments, run, validate
from run_fixed import write_members


def records(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def compare(expected, actual):
    assert len(expected) == len(actual), "event count changed"
    for a, b in zip(expected, actual):
        for field in ("kind", "state_id", "ell", "action", "rank", "coarse_marks",
                      "reference_marks", "coarse_fingerprint", "reference_fingerprint"):
            assert a.get(field) == b.get(field), (field, a.get(field), b.get(field))
        for field in ("theta", "eta", "solution", "targets", "frozen_scales"):
            if a.get(field) is not None:
                near(b[field], a[field], relative=1e-9, absolute=1e-11)


def execute(config, folder, executable, compute_only):
    if not compute_only:
        return run(config, folder, executable, timeout=600)
    rows = validate(config)
    folder.mkdir()
    table = (folder / "members.txt").resolve()
    write_members(table, rows)
    effective = (folder / "effective.json").resolve()
    effective.write_text(json.dumps(config, indent=2)+"\n")
    env = runtime_environment(config["threads"], supervised=False)
    env["ALOD_TIMING_FILE"] = str((folder / "timings.jsonl").resolve())
    with (folder / "solver.jsonl").open("w") as out, (folder / "stderr.log").open("w") as err:
        subprocess.run([str(executable.resolve()), *arguments(config, table), f"--config-file={effective}"],
                       env=env, stdout=out, stderr=err, timeout=600, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--level", type=int, default=5)
    parser.add_argument("--gap", type=int, default=3)
    parser.add_argument("--compute-only", action="store_true", help="omit checkpoint writes and supervisor/export overhead")
    parser.add_argument("--cases", nargs="+", choices=("E1", "E2", "E3"), default=["E1", "E2", "E3"])
    args = parser.parse_args()
    if args.repeats < 1 or args.threads < 2:
        parser.error("positive repeats and at least two threads are required")
    args.output.mkdir(parents=True, exist_ok=False)
    report = dict(threads=args.threads, level=args.level, gap=args.gap, repeats=args.repeats,
                  compute_only=args.compute_only,
                  scope="Two accepted states per run; sequential fresh runs; audits disabled; E3 nominal k=8", runs=[])
    modes = {"baseline": args.baseline, "updated": args.candidate, "serial": args.candidate}
    try:
        for case in args.cases:
            nominal = case == "E3"
            config = dict(DEFAULT, problem="E1" if nominal else case, wavenumber=8 if nominal else 16,
                          threads=args.threads, level=args.level, gap=args.gap, ell=2, maximum_ell=2,
                          ell_ratio_mode="solution_scaled", ell_threshold=0.3, ell_mode="every",
                          cycles=1, m_ref=1, reference_theta=[0.3], state_limit=2, minimum_gap=2, rank_cap=8,
                          member_ids=[0] if nominal else list(range(16)),
                          training_ids=[0] if nominal else list(range(16)),
                          maximum_nodes=1000000, maximum_patch_entries=640000000,
                          maximum_dense_entries=64000000, audit=False, emit_solution=not args.compute_only)
            expected = None
            for repetition in range(args.repeats):
                order = list(modes) if repetition % 2 == 0 else list(reversed(modes))
                for mode in order:
                    os.environ["ALOD_ASYNC_DISABLE"] = "1" if mode == "serial" else "0"
                    folder = args.output / f"{case}-{mode}-{repetition}"
                    before = resource.getrusage(resource.RUSAGE_CHILDREN)
                    start = time.monotonic()
                    execute(config, folder, modes[mode], args.compute_only)
                    wall = time.monotonic() - start
                    after = resource.getrusage(resource.RUSAGE_CHILDREN)
                    cpu = after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime
                    events = records(folder / "solver.jsonl")
                    if expected is None:
                        expected = events
                    else:
                        compare(expected, events)
                    timings = records(folder / "timings.jsonl")
                    phases = {}
                    for event in timings:
                        if event["event"] == "end":
                            phase = event["phase"]
                            phases[phase] = phases.get(phase, 0.) + event["seconds"]
                    teams = sorted({event["workers"] for event in timings if event["event"] == "team"})
                    row = dict(case=case, mode=mode, repetition=repetition, wall_seconds=wall,
                               cpu_seconds=cpu, average_cpu_percent=100*cpu/wall,
                               actual_team_sizes=teams, phases=phases, numerical_parity=True)
                    report["runs"].append(row)
                    (args.output / "report.json").write_text(json.dumps(report, indent=2)+"\n")
                    print(json.dumps(row), flush=True)
        report["medians"] = [dict(case=case, mode=mode,
                                   wall_seconds=statistics.median(r["wall_seconds"] for r in report["runs"]
                                                                  if r["case"] == case and r["mode"] == mode))
                             for case in args.cases for mode in modes]
        (args.output / "report.json").write_text(json.dumps(report, indent=2)+"\n")
        print(json.dumps(report["medians"]), flush=True)
    finally:
        os.environ.pop("ALOD_ASYNC_DISABLE", None)


if __name__ == "__main__":
    main()
