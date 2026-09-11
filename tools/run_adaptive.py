"""Run a bounded P3/P4 trajectory with explicit sample selection and event output."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time

from run_fixed import DEFAULT as FIXED, validate as validate_fixed, write_members

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = dict(FIXED, method="ALOD", cycles=1, m_ref=2, state_limit=0,
               minimum_gap=2, reference_theta=[0.3, 0.2], ell_mode="lazy",
               maximum_ell=4, extra_checks=[], force_promotions=[],
               radius=0.6, rank_cap=24, inherit=True, audit=False,
               emit_solution=False, maximum_patch_entries=8000000,
               maximum_dense_entries=8000000)


def validate(config):
    if not isinstance(config, dict) or set(config) != set(DEFAULT):
        raise ValueError("adaptive configuration must have exactly the documented fields")
    fixed = {key: config[key] for key in FIXED}
    # The common member/schema validator has a deliberately smaller P2 limit.
    cap = config["maximum_nodes"]
    if type(cap) is not int or not 1 <= cap <= 2000000:
        raise ValueError("invalid maximum_nodes")
    fixed["maximum_nodes"] = min(cap, 200000)
    rows = validate_fixed(fixed)
    for key, lo, hi in [("cycles", 0, 1000), ("m_ref", 1, 16), ("state_limit", 0, 10001),
                        ("minimum_gap", 0, 8), ("rank_cap", 1, 24),
                        ("maximum_ell", config["ell"], 4),
                        ("maximum_patch_entries", 1, 1000000000),
                        ("maximum_dense_entries", 1, 1000000000)]:
        if type(config[key]) is not int or not lo <= config[key] <= hi:
            raise ValueError(f"invalid {key}")
    for key in ("inherit", "audit", "emit_solution"):
        if type(config[key]) is not bool:
            raise ValueError(f"invalid {key}")
    if config["method"] not in ("ALOD", "AFEM") or config["ell_mode"] not in ("lazy", "every", "fixed"):
        raise ValueError("invalid method or ell_mode")
    if config["method"] == "AFEM" and config["m_ref"] != 1:
        raise ValueError("AFEM requires m_ref=1")
    for key in ("extra_checks", "force_promotions"):
        ids = config[key]
        if not isinstance(ids, list) or any(type(i) is not int or not 0 <= i <= 10000 for i in ids) or len(ids) != len(set(ids)):
            raise ValueError(f"invalid {key}")
    if config["force_promotions"] and (config["problem"] != "E2" or config["method"] != "ALOD" or config["ell_mode"] == "fixed"):
        raise ValueError("forced promotions are only an explicit E2 diagnostic")
    if 0 not in config["member_ids"]:
        raise ValueError("nominal sample 0 is required for localization")
    v = config["radius"]
    if type(v) not in (int, float) or not math.isfinite(v) or v < 0:
        raise ValueError("invalid radius")
    v = config["reference_theta"]
    if not isinstance(v, list) or len(v) != config["m_ref"] or any(type(x) not in (int, float) or not math.isfinite(x) or not 0 < x <= 1 for x in v):
        raise ValueError("reference_theta must contain m_ref fractions in (0,1]")
    return rows


def arguments(config, table):
    def text(value):
        if isinstance(value, list):
            return ",".join(map(str, value))
        return str(int(value)) if type(value) is bool else str(value)
    return ["adaptive", config["problem"], f"--members={table}"] + [
        f"--{key.replace('_', '-')}={text(value)}" for key, value in config.items()
        if key not in ("problem", "member_ids")]


def run(config, output, executable, timeout=3600):
    rows = validate(config)
    output.mkdir(parents=True, exist_ok=False)
    table = output / "members.txt"
    write_members(table, rows)
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    manifest = dict(config=config, table_sha256=sha(table), executable_sha256=sha(executable),
                    status="running", scope="P3/P4 adaptive trajectory; P5 checkpoint/resume not implemented")
    def save():
        (output / "run.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    save()
    start = time.monotonic()
    try:
        env = {**os.environ, "OMP_NUM_THREADS": str(config["threads"]), "OMP_PROC_BIND": "false",
               "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "BLIS_NUM_THREADS": "1"}
        env.pop("OMP_PLACES", None)
        with (output / "events.jsonl").open("w") as out, (output / "stderr.log").open("w") as err:
            result = subprocess.run([str(executable.resolve()), *arguments(config, table.resolve())],
                                    stdout=out, stderr=err, env=env, timeout=timeout)
        manifest["returncode"] = result.returncode
        if result.returncode:
            raise RuntimeError("adaptive run failed; see stderr.log and completed event prefix")
        events = [json.loads(line, parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
                  for line in (output / "events.jsonl").read_text().splitlines()]
        accepted = [e for e in events if e["kind"] == "accepted"]
        expected = 1 + config["cycles"] * config["m_ref"]
        expected = min(expected, config["state_limit"] or expected)
        if [e["state_id"] for e in accepted] != list(range(expected)):
            raise RuntimeError("accepted-state sequence is incomplete")
        manifest.update(status="complete", accepted_states=len(accepted),
                        ell_checks=sum(e["kind"] == "ell_check" for e in events),
                        complete_horizon=accepted[-1]["complete_horizon"],
                        interrupted_cycle=not accepted[-1]["cycle_complete"])
    except BaseException:
        manifest["status"] = "failed"
        raise
    finally:
        # Preserve the last accepted cursor even when a later resource or
        # numerical gate stops the run. This is a journal, not a checkpoint.
        log = output / "events.jsonl"
        accepted_prefix = []
        if log.exists():
            for line in log.read_text().splitlines():
                try:
                    event = json.loads(line)
                except ValueError:
                    break
                if event.get("kind") == "accepted":
                    accepted_prefix.append(event)
        if accepted_prefix:
            last = accepted_prefix[-1]
            manifest["last_accepted_cursor"] = {key: last[key] for key in
                ("state_id", "coarse_cycle", "reference_sweep", "cycle_complete", "next_phase", "mesh_revision", "ell_last_check")}
        manifest["wall_seconds"] = time.monotonic() - start
        save()
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--executable", type=Path, default=ROOT / "build/alod_run")
    parser.add_argument("--timeout", type=float, default=3600)
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("timeout must be positive and finite")
    print(json.dumps(run(json.loads(args.config.read_text()), args.output, args.executable, args.timeout)))


if __name__ == "__main__":
    main()
