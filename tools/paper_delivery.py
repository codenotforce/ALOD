"""Generate the version 4 paper package from native, archive or canonical data."""

import argparse
import json
from pathlib import Path
from paper.data import read, native, import_archives
from paper.bundle import deliver


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--spec", type=Path, help="Schema 1 JSON with runs mapping names to native directories."
    )
    p.add_argument("--canonical", type=Path)
    p.add_argument(
        "--executable", type=Path, default=Path(__file__).resolve().parents[1] / "build/alod_run"
    )
    p.add_argument("--archive-e1", type=Path)
    p.add_argument("--archive-e2", type=Path)
    p.add_argument("--archive-e3", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--allow-partial", action="store_true")
    a = p.parse_args()
    if sum((a.spec is not None, a.canonical is not None, a.archive_e1 is not None)) != 1:
        p.error("select one input mode")
    if a.canonical:
        data = read(a.canonical)
        if data.get("schema") != 1:
            raise ValueError("unsupported canonical schema")
        runs = data["runs"]
    elif a.spec:
        spec = read(a.spec)
        if spec.get("schema") != 1:
            raise ValueError("unsupported input spec")
        runs = {
            name: native(a.spec.parent / path, a.executable if a.executable.exists() else None)
            for name, path in spec["runs"].items()
        }
    else:
        if not a.archive_e2 or not a.archive_e3:
            p.error("all three archive roots are required")
        runs = import_archives(a.archive_e1, a.archive_e2, a.archive_e3)
    result = deliver(runs, a.output, allow_partial=a.allow_partial)
    print(json.dumps({k: result[k] for k in ("complete", "missing_runs", "issues")}))


if __name__ == "__main__":
    main()
