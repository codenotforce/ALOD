"""Alternate identical phase probes linked against before/after core libraries (Linux)."""
import argparse
import hashlib
import json
import statistics
import shutil
import subprocess
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("before", "after", "members", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    for name, default in (("level", 8), ("gap", 4), ("ell", 2), ("threads", 8),
                          ("wavenumber", 32), ("repeats", 3)):
        parser.add_argument(f"--{name}", type=int, default=default)
    parser.add_argument("--problem", choices=("E1", "E2"), default="E2")
    args = parser.parse_args()
    time_program = shutil.which("time")
    if time_program is None:
        parser.error("the GNU time executable is required for peak RSS measurements")
    if args.repeats < 1 or args.threads < 1:
        parser.error("repeats and threads must be positive")
    args.output.mkdir(parents=True, exist_ok=False)
    options = [f"--{k}={getattr(args,k)}" for k in
               ("level", "gap", "ell", "threads", "wavenumber")]
    options += [f"--members={args.members}", "--maximum-nodes=4000000"]
    result = {"configuration": {k: getattr(args, k) for k in
              ("problem", "level", "gap", "ell", "threads", "wavenumber", "repeats")},
              "members_sha256": hashlib.sha256(args.members.read_bytes()).hexdigest(),
              "binary_sha256": {name: hashlib.sha256(getattr(args, name).read_bytes()).hexdigest()
                                for name in ("before", "after")}, "runs": []}
    expected = None
    for repeat in range(args.repeats + 1):
        for name in (("before", "after") if repeat % 2 == 0 else ("after", "before")):
            stem = args.output / f"{repeat}-{name}"
            start = time.monotonic()
            completed = subprocess.run([time_program, "-f", "%M", "-o", str(stem)+".rss",
                                        str(getattr(args, name)), args.problem, *options],
                                       capture_output=True, text=True, check=True)
            row = json.loads(completed.stdout)
            elapsed = time.monotonic() - start
            semantic = {k: v for k, v in row.items() if not k.endswith("_seconds")}
            if expected is None:
                expected = semantic
            if semantic != expected:
                raise RuntimeError(f"space/training results differ: {name}, repeat {repeat}")
            row.update(version=name, warmup=repeat == 0, elapsed_seconds=elapsed,
                       peak_rss_kib=int(Path(str(stem)+".rss").read_text().strip()))
            stem.with_suffix(".json").write_text(json.dumps(row, indent=2)+"\n", encoding="utf-8")
            result["runs"].append(row)
            print(json.dumps(row), flush=True)
    result["identical_numerical_outputs"] = True
    keys = [k for k in result["runs"][0] if k.endswith("_seconds") or k == "peak_rss_kib"]
    result["medians"] = {name: {k: statistics.median(r[k] for r in result["runs"]
                                  if r["version"] == name and not r["warmup"]) for k in keys}
                         for name in ("before", "after")}
    (args.output / "summary.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
