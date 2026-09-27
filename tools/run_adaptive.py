"""Run or resume a checkpointed trajectory, with concurrent checkpoint-backed audits."""
import argparse
import hashlib
import json
import math
import os
import sys
import tempfile
from pathlib import Path
import subprocess
import time

from adaptive_config import with_compatibility_defaults
from checkpoint_io import atomic_text, digest, inspect, recover_journal, validate_resume_config
from run_fixed import DEFAULT as FIXED, validate as validate_fixed, write_members, member_text

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = dict(FIXED, exact_target=-1., exact_scope="nominal", audit_mode="full", method="ALOD", cycles=1, m_ref=2, state_limit=0,
               minimum_gap=2, reference_theta=[0.3, 0.2], ell_mode="lazy",
               ell_ratio_mode="raw", ell_threshold=0., ell_absolute_threshold=-1., enrichment_tests="kernel_lift",
               maximum_ell=4, extra_checks=[], force_promotions=[],
               radius=0.6, rank_cap=24, inherit=True, audit=False,
               emit_solution=False, maximum_patch_entries=8000000,
               maximum_dense_entries=8000000)


def validate(config):
    if isinstance(config, dict): config = with_compatibility_defaults(config)
    if not isinstance(config, dict) or set(config) != set(DEFAULT):
        raise ValueError("adaptive configuration must have exactly the documented fields")
    if config["ell_ratio_mode"] not in ("raw","solution_scaled") or config["enrichment_tests"] not in ("kernel_lift","adjoint"):
        raise ValueError("invalid localization or enrichment policy")
    if type(config["ell_threshold"]) not in (int,float) or not math.isfinite(config["ell_threshold"]) or config["ell_threshold"]<0:
        raise ValueError("invalid ell_threshold")
    tau=config["ell_absolute_threshold"]
    if type(tau) not in (int,float) or not math.isfinite(tau) or (tau<0 and tau!=-1):
        raise ValueError("ell_absolute_threshold must be nonnegative or -1 (disabled)")
    if type(config['exact_target']) not in (int,float) or not math.isfinite(config['exact_target']) or (config['exact_target']<=0 and config['exact_target']!=-1):raise ValueError('invalid exact target')
    if config['exact_scope'] not in ('nominal','training_max') or config['audit_mode'] not in ('full','exact'):raise ValueError('invalid exact scope or audit mode')
    fixed = {key: config[key] for key in FIXED}
    # The common member/schema validator has a deliberately smaller P2 limit.
    cap = config["maximum_nodes"]
    if type(cap) is not int or not 1 <= cap <= 4000000:
        raise ValueError("invalid maximum_nodes")
    fixed["maximum_nodes"] = min(cap, 200000)
    rows = validate_fixed(fixed)
    for key, lo, hi in [("cycles", 0, 1000), ("m_ref", 1, 16), ("state_limit", 0, 10001),
                        ("minimum_gap", 0, 8), ("rank_cap", 1, 24),
                        ("maximum_ell", config["ell"], 4),
                        ("maximum_patch_entries", 1, 64000000000),
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
        if key not in ("problem", "member_ids", "audit_mode")]


def run(config, output, executable, timeout=3600, *, resume=None, pause_state=None, pause_phase="accepted", audit_workers=1, audit_threads=None, audit_drain_workers=None):
    from execution import run_lease
    with run_lease(output):
        return _run(config,output,executable,timeout,resume=resume,pause_state=pause_state,pause_phase=pause_phase,audit_workers=audit_workers,audit_threads=audit_threads,audit_drain_workers=audit_drain_workers)


def _run(config, output, executable, timeout=3600, *, resume=None, pause_state=None, pause_phase="accepted", audit_workers=1, audit_threads=None, audit_drain_workers=None):
    config = with_compatibility_defaults(config)
    setup_start = time.monotonic()
    output, executable = Path(output), Path(executable)
    rows = validate(config)
    if type(audit_workers) is not int or audit_workers<1:raise ValueError("audit workers must be positive")
    if audit_threads is not None and (type(audit_threads) is not int or not 0<=audit_threads<=2147483647):raise ValueError("invalid audit threads")
    effective_audit_threads=config["threads"] if audit_threads is None else audit_threads
    if audit_drain_workers is None:
        audit_drain_workers=audit_workers+(min(audit_workers,config["threads"]//effective_audit_threads) if effective_audit_threads>0 and config["threads"]>0 else 0)
    if type(audit_drain_workers) is not int or audit_drain_workers<audit_workers:raise ValueError("audit drain workers must be at least active workers")
    metadata = None
    startup_recovery=False
    if str(resume)=="auto" and not list((output/'checkpoints').rglob('*.bin')):
        previous=json.loads((output/'run.json').read_text())
        validate_resume_config(previous['config'],config)
        if (output/'members.txt').read_text()!=member_text(rows):raise ValueError("startup recovery member table mismatch")
        log=output/'solver.jsonl'
        if log.exists() and log.stat().st_size:raise ValueError("nonempty numerical journal without a checkpoint")
        startup_recovery=True;resume=None
    if resume is not None:
        if str(resume)=="auto":
            from checkpoint_io import recover_latest
            resume=recover_latest(output,executable,config,member_text(rows))
        resume, metadata = inspect(resume, executable)
        validate_resume_config(metadata["config"], config)
        if metadata['members_text'] != member_text(rows):
            raise ValueError('resume member table differs from the frozen checkpoint')
        if not output.is_dir():
            raise ValueError("resume requires the original run directory and its committed journal")
        recover_journal(output/"solver.jsonl", metadata)
    elif not startup_recovery:
        output.mkdir(parents=True, exist_ok=False)
    if pause_state is not None and (type(pause_state) is not int or pause_state < 0 or pause_phase not in ("accepted", "training", "ell")):
        raise ValueError("invalid checkpoint pause request")
    table = output / "members.txt"
    write_members(table, rows)
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    manifest = dict(config=config, table_sha256=sha(table), executable_sha256=sha(executable),
                    status="running", scope="P5 checkpointed adaptive trajectory", solver_completed=False, validation_passed=False, audit_complete=False, paper_complete=False)
    mathematical = {key:value for key,value in config.items() if key not in
                    ('cycles','state_limit','maximum_nodes','maximum_patch_entries','maximum_dense_entries','threads','audit','audit_mode','emit_solution')}
    manifest['experiment_id'] = hashlib.sha256((json.dumps(mathematical,sort_keys=True)+member_text(rows)).encode()).hexdigest()[:24]
    from runtime_provenance import provenance
    manifest['provenance'] = provenance(executable)
    manifest['startup_recovery']=startup_recovery
    manifest['setup_seconds'] = time.monotonic()-setup_start
    if metadata is not None:
        manifest["resumed_from"] = dict(checkpoint=resume.name, sha256=digest(resume), state_id=metadata["state_id"], phase=metadata["phase"])
    def save():
        atomic_text(output / "run.json", json.dumps(manifest, indent=2) + "\n")
    save()
    start = time.monotonic()
    queue=None;process=None;socket_directory=None;shared_endpoint=None
    try:
        from execution import runtime_environment
        env=runtime_environment(config["threads"])
        env["ALOD_TIMING_FILE"]=str((output/"timings.jsonl").resolve())
        if config["audit"]:
            from async_audit import AuditQueue
            if (sys.platform=='linux' and config['method']=='ALOD' and not os.environ.get('ALOD_REFERENCE_EXECUTION')
                    and not os.environ.get('ALOD_AUDIT_ONESHOT') and not os.environ.get('ALOD_AUDIT_DISK_ONLY')):
                socket_directory=tempfile.TemporaryDirectory(prefix='alod-audit-')
                shared_endpoint=Path(socket_directory.name)/'worker.sock'
                env['ALOD_AUDIT_SOCKET']=str(shared_endpoint)
                env['ALOD_SHARED_AUDIT_WORKERS']=str(audit_drain_workers)
            queue=AuditQueue(output,executable,config["member_ids"],audit_workers,audit_threads,drain_workers=audit_drain_workers,shared_endpoint=shared_endpoint,audit_mode=config["audit_mode"])
            queue.discover()
        effective = output / "effective.json"
        atomic_text(effective, json.dumps(config, indent=2)+"\n")
        runtime = [f"--checkpoint-dir={(output/'checkpoints').resolve()}", f"--config-file={effective.resolve()}"]
        if resume is not None:
            runtime.append(f"--resume={resume.resolve()}")
        if pause_state is not None:
            runtime += [f"--pause-state={pause_state}", f"--pause-phase={pause_phase}"]
        solver_start = time.monotonic()
        with (output / "solver.jsonl").open("a" if resume is not None else "w") as out, (output / "stderr.log").open("a" if resume is not None else "w") as err:
            process = subprocess.Popen([str(executable.resolve()), *arguments(config, table.resolve()), *runtime],
                                       stdout=out, stderr=err, env=env)
            while process.poll() is None:
                if queue:queue.discover()
                if shared_endpoint is not None and Path(str(shared_endpoint)+'.done').exists():break
                if time.monotonic()-solver_start>timeout:raise subprocess.TimeoutExpired(process.args,timeout)
                time.sleep(.1)
            result=subprocess.CompletedProcess(process.args,0) if process.poll() is None and shared_endpoint is not None else process
            if queue:queue.discover()
        manifest['solver_wall_seconds'] = time.monotonic()-solver_start
        manifest['solver_finished_at']=time.time()
        manifest["returncode"] = result.returncode
        if result.returncode:
            raise RuntimeError("adaptive run failed; see stderr.log and completed event prefix")
        if queue:queue.begin_drain()
        latest, last_metadata = inspect(output/"checkpoints/latest", executable)
        recover_journal(output/"solver.jsonl", last_metadata)
        atomic_text(output/"events.jsonl", (output/"solver.jsonl").read_text())
        events = [json.loads(line, parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
                  for line in (output / "events.jsonl").read_text().splitlines()]
        accepted = [e for e in events if e["kind"] == "accepted"]
        expected = 1 + config["cycles"] * config["m_ref"]
        expected = min(expected, config["state_limit"] or expected)
        target_reached=bool(accepted and accepted[-1].get('target_reached'))
        if target_reached:
            ratio=accepted[-1].get('exact_target_ratio')
            if config['exact_target']<=0 or not isinstance(ratio,(int,float)) or not 0<=ratio<=config['exact_target']:raise ValueError('invalid exact target evidence')
            manifest['target_state_id']=accepted[-1]['state_id']
            manifest['target_ratio']=ratio
        manifest['target_reached']=target_reached
        manifest['termination_reason']='exact_target' if target_reached else 'horizon'
        paused = not target_reached and pause_state is not None and len(accepted) < expected
        if [e["state_id"] for e in accepted] != list(range(len(accepted))) or (not paused and not target_reached and len(accepted) != expected):
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
            if not paused:manifest["status"]="auditing"
            save()
            audit_failures = queue.finish(accepted)
            if shared_endpoint is not None and process.poll() is None:
                from native_audit import shutdown_shared
                shutdown_shared(shared_endpoint)
                process.wait(timeout=30)
                if process.returncode:raise RuntimeError('shared audit server failed during drain')
            if not paused:manifest["status"]="complete"
            manifest["audit_complete"] = not audit_failures and (len(accepted)==expected or target_reached)
            manifest["audit_failures"] = audit_failures
            if audit_failures:
                manifest["status"] = "audit_failed"
            atomic_text(output/"events.jsonl", "".join(json.dumps(event,separators=(",",":"))+"\n" for event in events))
        manifest['audit_drain_seconds'] = time.monotonic()-audit_start if config['audit'] else 0
        manifest['audit_wall_seconds'] = time.monotonic()-start if config['audit'] else 0
        manifest['audit_workers']=audit_workers;manifest['audit_threads']=audit_threads
        manifest['audit_drain_workers']=audit_drain_workers
        manifest['audit_execution']='shared_memory_and_persistent_workers' if shared_endpoint else 'persistent_checkpoint_workers'
        manifest['checkpoint_bytes'] = sum(path.stat().st_size for path in (output/'checkpoints').rglob('*.bin'))
        from export_results import export
        try:
            manifest["canonical_tables"] = export(output)
        except Exception as error:
            manifest["status"] = "postprocess_failed"
            manifest["postprocess_error"] = str(error)
        manifest["paper_complete"] = bool(config["audit_mode"]=="full" and "canonical_tables" in manifest and manifest["solver_completed"] and manifest["validation_passed"] and manifest["audit_complete"] and manifest["complete_horizon"] and config["member_ids"] == list(range(48)) and len(accepted) >= (33 if config["problem"]=="E2" and config["enrichment_tests"]=="adjoint" else 97 if config["method"]=="AFEM" else 51))
    except BaseException:
        manifest["status"] = "failed"
        raise
    finally:
        if process is not None and process.poll() is None:
            process.terminate()
            try:process.wait(timeout=5)
            except subprocess.TimeoutExpired:process.kill();process.wait()
        if queue:queue.close()
        if socket_directory is not None:socket_directory.cleanup()
        # Preserve the last accepted cursor even when a later resource or
        # numerical gate stops the run; checkpoints remain independently usable.
        log = output / "solver.jsonl"
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
    import signal
    def stop(signum,frame):raise KeyboardInterrupt("run interrupted; resume with --resume auto")
    signal.signal(signal.SIGTERM,stop)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--executable", type=Path, default=ROOT / "build/alod_run")
    parser.add_argument("--resume", type=Path, help="checkpoint path or auto for validated crash recovery")
    parser.add_argument("--audit-workers", type=int, default=1)
    parser.add_argument("--audit-drain-workers", type=int, help="workers after solver exit; default reuses its explicit thread budget; set equal to --audit-workers to keep memory demand unchanged")
    parser.add_argument("--audit-threads", type=int, help="0 inherits runtime; omitted inherits run threads")
    parser.add_argument("--pause-state", type=int)
    parser.add_argument("--pause-phase", choices=["accepted", "training", "ell"], default="accepted")
    parser.add_argument("--timeout", type=float, default=3600)
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("timeout must be positive and finite")
    print(json.dumps(run(json.loads(args.config.read_text()), args.output, args.executable, args.timeout, resume=args.resume, pause_state=args.pause_state, pause_phase=args.pause_phase,audit_workers=args.audit_workers,audit_threads=args.audit_threads,audit_drain_workers=args.audit_drain_workers)))


if __name__ == "__main__":
    main()
