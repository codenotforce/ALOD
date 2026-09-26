"""Alternating before/after E2 workflow benchmark with semantic output checks."""
import argparse
import hashlib
import json
import statistics
from pathlib import Path

from run_adaptive import run


def records(path):
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def semantic(value):
    if isinstance(value, dict):
        return {k: semantic(v) for k, v in value.items()
                if not k.endswith("_seconds") and k not in {"checkpoint", "checkpoint_bytes"}}
    if isinstance(value, list):
        return [semantic(v) for v in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--problem", choices=("E1", "E2"), default="E2")
    parser.add_argument("--level", type=int, default=6)
    parser.add_argument("--gap", type=int, default=4)
    args = parser.parse_args()
    if args.repeats < 1 or args.threads < 1:
        parser.error("repeats and threads must be positive")
    args.output.mkdir(parents=True, exist_ok=False)
    config = json.loads(Path(f"configs/adaptive/{args.problem.lower()}_main.json").read_text(encoding="utf-8"))
    config.update(cycles=1, threads=args.threads, audit=True,
                  level=args.level, gap=args.gap,
                  member_ids=list(range(16)), training_ids=list(range(16)))
    report = {"config": config, "repeats": args.repeats,
              "binary_sha256": {name: hashlib.sha256(path.read_bytes()).hexdigest()
                                for name, path in (("before", args.before), ("after", args.after))},
              "runs": [], "accepted_semantics_identical": True}
    baseline = None
    audit_baseline = None
    for repeat in range(args.repeats + 1):
        order = ("before", "after") if repeat % 2 == 0 else ("after", "before")
        for name in order:
            destination = args.output / f"{repeat}-{name}"
            manifest = run(config, destination, getattr(args, name),
                           audit_workers=1, audit_threads=args.threads)
            accepted = [r for r in records(destination / "solver.jsonl") if r["kind"] == "accepted"]
            current = semantic(accepted)
            if baseline is None:
                baseline = current
            if baseline != current:
                raise RuntimeError("accepted-state semantics differ between benchmark runs")
            audit_records = [r for path in (destination / "audits").glob("*/samples.jsonl")
                             for r in records(path)]
            audit = [r for r in audit_records if r["kind"] == "audit_complete"]
            samples = sorted((semantic(r) for r in audit_records if r["kind"] == "sample"),
                             key=lambda r: (r["state_id"], r["sample"]))
            if audit_baseline is None:
                audit_baseline = samples
            if samples != audit_baseline:
                raise RuntimeError("audit sample semantics differ between benchmark runs")
            report["audit_samples_identical"] = True
            if not manifest["audit_complete"] or len(audit) != len(accepted):
                raise RuntimeError("benchmark audits did not complete")
            row = {"version": name, "warmup": repeat == 0,
                   "solver_seconds": manifest["solver_wall_seconds"],
                   "workflow_seconds": manifest["wall_seconds"],
                   "load_seconds": sum(r["load_seconds"] for r in accepted),
                   "strong_seconds": sum(r["strong_seconds"] for r in accepted),
                   "error_seconds": sum(r["error_seconds"] for r in audit)}
            report["runs"].append(row)
            print(json.dumps(row), flush=True)
    report["medians"] = {name: {key: statistics.median(r[key] for r in report["runs"]
                          if r["version"] == name and not r["warmup"])
                         for key in report["runs"][0] if key.endswith("_seconds")}
                         for name in ("before", "after")}
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
