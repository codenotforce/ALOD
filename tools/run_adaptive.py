"""Run or resume a checkpointed trajectory, then schedule independent audits."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time

from checkpoint_io import atomic_text, digest, inspect, recover_journal, validate_resume_config
from run_fixed import DEFAULT as FIXED, validate as validate_fixed, write_members, member_text

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


def run(config, output, executable, timeout=3600, *, resume=None, pause_state=None, pause_phase="accepted"):
    setup_start = time.monotonic()
    output, executable = Path(output), Path(executable)
    rows = validate(config)
    metadata = None
    if resume is not None:
        resume, metadata = inspect(resume, executable)
        validate_resume_config(metadata["config"], config)
        if metadata['members_text'] != member_text(rows):
            raise ValueError('resume member table differs from the frozen checkpoint')
        if not output.is_dir():
            raise ValueError("resume requires the original run directory and its committed journal")
        recover_journal(output/"solver.jsonl", metadata)
    else:
        output.mkdir(parents=True, exist_ok=False)
    if pause_state is not None and (type(pause_state) is not int or pause_state < 0 or pause_phase not in ("accepted", "training", "ell")):
        raise ValueError("invalid checkpoint pause request")
    table = output / "members.txt"
    write_members(table, rows)
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    manifest = dict(config=config, table_sha256=sha(table), executable_sha256=sha(executable),
                    status="running", scope="P5 checkpointed adaptive trajectory", solver_completed=False, validation_passed=False, audit_complete=False, paper_complete=False)
    mathematical = {key:value for key,value in config.items() if key not in
                    ('cycles','state_limit','maximum_nodes','maximum_patch_entries','maximum_dense_entries','threads','audit','emit_solution')}
    manifest['experiment_id'] = hashlib.sha256((json.dumps(mathematical,sort_keys=True)+member_text(rows)).encode()).hexdigest()[:24]
    from runtime_provenance import provenance
    manifest['provenance'] = provenance(executable)
    manifest['setup_seconds'] = time.monotonic()-setup_start
    if metadata is not None:
        manifest["resumed_from"] = dict(checkpoint=resume.name, sha256=digest(resume), state_id=metadata["state_id"], phase=metadata["phase"])
    def save():
        atomic_text(output / "run.json", json.dumps(manifest, indent=2) + "\n")
    save()
    start = time.monotonic()
    try:
        env = {**os.environ, "OMP_NUM_THREADS": str(config["threads"]), "OMP_PROC_BIND": "false",
               "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "BLIS_NUM_THREADS": "1"}
        env.pop("OMP_PLACES", None)
        effective = output / "effective.json"
        atomic_text(effective, json.dumps(config, indent=2)+"\n")
        runtime = [f"--checkpoint-dir={(output/'checkpoints').resolve()}", f"--config-file={effective.resolve()}"]
        if resume is not None:
            runtime.append(f"--resume={resume.resolve()}")
        if pause_state is not None:
            runtime += [f"--pause-state={pause_state}", f"--pause-phase={pause_phase}"]
        solver_start = time.monotonic()
        with (output / "solver.jsonl").open("a" if resume is not None else "w") as out, (output / "stderr.log").open("a" if resume is not None else "w") as err:
            result = subprocess.run([str(executable.resolve()), *arguments(config, table.resolve()), *runtime],
                                    stdout=out, stderr=err, env=env, timeout=timeout)
        manifest['solver_wall_seconds'] = time.monotonic()-solver_start
        manifest["returncode"] = result.returncode
        if result.returncode:
            raise RuntimeError("adaptive run failed; see stderr.log and completed event prefix")
        latest, last_metadata = inspect(output/"checkpoints/latest", executable)
        recover_journal(output/"solver.jsonl", last_metadata)
        atomic_text(output/"events.jsonl", (output/"solver.jsonl").read_text())
        events = [json.loads(line, parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
                  for line in (output / "events.jsonl").read_text().splitlines()]
        accepted = [e for e in events if e["kind"] == "accepted"]
        expected = 1 + config["cycles"] * config["m_ref"]
        expected = min(expected, config["state_limit"] or expected)
        paused = pause_state is not None and len(accepted) < expected
        if [e["state_id"] for e in accepted] != list(range(len(accepted))) or (not paused and len(accepted) != expected):
            raise RuntimeError("accepted-state sequence is incomplete")
        manifest.update(status="paused" if paused else "complete", solver_completed=not paused, validation_passed=True, accepted_states=len(accepted),
                        ell_checks=sum(e["kind"] == "ell_check" for e in events),
                        complete_horizon=accepted[-1]["complete_horizon"] if accepted else False,
                        interrupted_cycle=not accepted[-1]["cycle_complete"] if accepted else False, checkpoint=latest.name)
        # Audit errors never invalidate or terminate the numerical solver.
        # Each audit is a separate process and can be retried from its snapshot.
        audit_failures = []
        audit_start = time.monotonic()
        if config["audit"]:
            from run_audit import audit
            for state in accepted:
                checkpoint = output/"checkpoints"/f"state-{state['state_id']:06d}-phase-0-ell-{state['ell'] or config['ell']}.bin"
                try:
                    audit_output = output/"audits"/(digest(checkpoint)[:16]+"-shared")
                    if audit_output.exists():
                        prior = json.loads((audit_output/'run.json').read_text()) if (audit_output/'run.json').exists() else {}
                        if not prior.get('audit_complete'):
                            original = audit_output
                            attempt = 1
                            while audit_output.exists():
                                audit_output = original.with_name(original.name+f'-retry-{attempt}')
                                attempt += 1
                    if not (audit_output/"run.json").exists():
                        audit(checkpoint, audit_output, executable.with_name("alod_audit"), batch_size=8)
                    audit_manifest = json.loads((audit_output/"run.json").read_text())
                    if not audit_manifest["audit_complete"] or audit_manifest["checkpoint_sha256"] != digest(checkpoint):
                        raise ValueError("cached audit did not pass validation")
                    if audit_manifest['output_sha256'] != digest(audit_output/'samples.jsonl') or audit_manifest['member_ids'] != config['member_ids']:
                        raise ValueError('cached audit content or coverage changed')
                    state["audit"] = [json.loads(line) for line in (audit_output/"samples.jsonl").read_text().splitlines()
                                      if json.loads(line)["kind"] == "sample"]
                    if len(state['audit']) != len(config['member_ids']) or {row['sample'] for row in state['audit']} != set(config['member_ids']) or any(row['state_id'] != state['state_id'] for row in state['audit']):
                        del state['audit']
                        raise ValueError('cached audit state/sample mismatch')
                except Exception as error:
                    audit_failures.append(dict(state_id=state["state_id"], reason=str(error)))
            manifest["audit_complete"] = not audit_failures and len(accepted)==expected
            manifest["audit_failures"] = audit_failures
            if audit_failures:
                manifest["status"] = "audit_failed"
            atomic_text(output/"events.jsonl", "".join(json.dumps(event,separators=(",",":"))+"\n" for event in events))
        manifest['audit_wall_seconds'] = time.monotonic()-audit_start if config['audit'] else 0
        manifest['checkpoint_bytes'] = sum(path.stat().st_size for path in (output/'checkpoints').glob('*.bin'))
        from export_results import export
        try:
            manifest["canonical_tables"] = export(output)
        except Exception as error:
            manifest["status"] = "postprocess_failed"
            manifest["postprocess_error"] = str(error)
        manifest["paper_complete"] = bool("canonical_tables" in manifest and manifest["solver_completed"] and manifest["validation_passed"] and manifest["audit_complete"] and manifest["complete_horizon"] and config["member_ids"] == list(range(48)) and len(accepted) >= (33 if config["problem"]=="E2" else 97 if config["method"]=="AFEM" else 51))
    except BaseException:
        manifest["status"] = "failed"
        raise
    finally:
        # Preserve the last accepted cursor even when a later resource or
        # numerical gate stops the run; checkpoints remain independently usable.
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
    parser.add_argument("--resume", type=Path)
    parser.add_argument("--pause-state", type=int)
    parser.add_argument("--pause-phase", choices=["accepted", "training", "ell"], default="accepted")
    parser.add_argument("--timeout", type=float, default=3600)
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("timeout must be positive and finite")
    print(json.dumps(run(json.loads(args.config.read_text()), args.output, args.executable, args.timeout, resume=args.resume, pause_state=args.pause_state, pause_phase=args.pause_phase)))


if __name__ == "__main__":
    main()
