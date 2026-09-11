"""Run one bounded P2 fixed state using the frozen explicit member table."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = dict(problem="E1", level=2, gap=3, ell=1, graded=False, riesz_patches="n2",
               interpolation="area", threads=1, maximum_nodes=20000,
               theta=0.15, ritz_tolerance=1e-4, ritz_iterations=750,
               dense_threshold=64, member_ids=list(range(16)), training_ids=list(range(16)))


def validate(config):
    if not isinstance(config, dict) or set(config) != set(DEFAULT):
        raise ValueError("fixed-state configuration must have exactly the documented fields")
    if config["problem"] not in ("E1", "E2") or config["interpolation"] not in ("area", "arithmetic"):
        raise ValueError("invalid problem or interpolation policy")
    if config["riesz_patches"] not in ("n2", "archive"):
        raise ValueError("invalid Riesz patch policy")
    if type(config["graded"]) is not bool:
        raise ValueError("graded must be boolean")
    for key, lo, hi in [("level", 0, 12), ("gap", 1, 8), ("ell", 1, 4), ("threads", 1, 64),
                        ("maximum_nodes", 1, 200000), ("ritz_iterations", 1, 10000), ("dense_threshold", 0, 512)]:
        if type(config[key]) is not int or not lo <= config[key] <= hi:
            raise ValueError(f"invalid {key}")
    for key in ("theta", "ritz_tolerance"):
        if type(config[key]) not in (int, float) or not math.isfinite(config[key]) or not 0 < config[key] <= 1:
            raise ValueError(f"invalid {key}")
    rows = json.loads((ROOT / "data/rhs" / f"rhs_{config['problem'].lower()}.json").read_text())
    table = {r["sample"]: r for r in rows}
    for key in ("member_ids", "training_ids"):
        ids = config[key]
        if not isinstance(ids, list) or not ids or any(type(i) is not int for i in ids) or len(ids) != len(set(ids)):
            raise ValueError(f"{key} must be a nonempty list of distinct integers")
        if any(i not in table for i in ids):
            raise ValueError("unknown member ID")
    if any(i not in config["member_ids"] or table[i]["split"] != "train" for i in config["training_ids"]):
        raise ValueError("marking may select only included training members")
    return [table[i] for i in config["member_ids"]]


def member_text(rows):
    lines = []
    for r in rows:
        fields = [r["sample"], r["split"], r["wave_x"], r["wave_y"],
                  r.get("singular_coefficient", 1), r.get("wave_amplitude", 1), r.get("wave_phase", 0)]
        lines.append(" ".join(str(x) for x in fields))
    return "\n".join(lines) + "\n"


def write_members(path, rows):
    path.write_text(member_text(rows), encoding="utf-8")


def arguments(config, table):
    return [config["problem"], f"--members={table}", "--training-ids=" + ",".join(map(str, config["training_ids"]))] + [
        f"--{key.replace('_', '-')}={int(value) if type(value) is bool else value}"
        for key, value in config.items() if key not in ("problem", "member_ids", "training_ids")]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--executable", type=Path, default=ROOT / "build/alod_fixed")
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    rows = validate(config)
    args.output.mkdir(parents=True, exist_ok=False)
    table = args.output / "members.txt"
    write_members(table, rows)
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    manifest = dict(config=config, table_sha256=sha(table), executable_sha256=sha(args.executable),
                    scope="P2 fixed state; no adaptive scheduler or production trajectory", status="running")
    start = time.monotonic()
    try:
        env = {**os.environ, "OMP_NUM_THREADS": str(config["threads"]), "OPENBLAS_NUM_THREADS": "1",
               "MKL_NUM_THREADS": "1", "BLIS_NUM_THREADS": "1", "OMP_PROC_BIND": "false"}
        env.pop("OMP_PLACES", None)
        with (args.output / "state.json").open("w") as out, (args.output / "stderr.log").open("w") as err:
            run = subprocess.run([str(args.executable.resolve()), *arguments(config, table.resolve())],
                                 stdout=out, stderr=err, env=env, timeout=3600)
        manifest["status"] = "complete" if run.returncode == 0 else "failed"
        manifest["returncode"] = run.returncode
        if run.returncode:
            raise RuntimeError("fixed-state solve failed; see stderr.log")
    except BaseException:
        manifest["status"] = "failed"
        raise
    finally:
        manifest["wall_seconds"] = time.monotonic() - start
        (args.output / "run.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
