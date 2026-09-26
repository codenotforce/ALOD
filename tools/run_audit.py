"""Audit immutable accepted checkpoints without replaying an adaptive prefix."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

from checkpoint_io import atomic_text, digest, inspect
from run_fixed import DEFAULT as FIXED, arguments, validate, write_members

ROOT = Path(__file__).resolve().parents[1]


class AuditExecutionContext:
    """Run-scoped executable/source provenance shared by asynchronous jobs.

    This is a launch snapshot, not a global cache. Replacing the executable
    invalidates the context; checkpoint validation remains per-job.
    """
    def __init__(self, executable):
        from runtime_provenance import provenance
        self.executable = Path(executable).resolve()
        self.identity = self._identity()
        self.executable_sha256 = digest(self.executable)
        self.provenance = provenance(self.executable.with_name('alod_run'))
        self.verify(executable)

    def _identity(self):
        s = self.executable.stat()
        return (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)

    def verify(self, executable):
        if Path(executable).resolve() != self.executable or self._identity() != self.identity:
            raise ValueError('audit executable changed after launch snapshot')


def audit(checkpoint, output, executable, *, member_ids=None, batch_size=8, fresh=False, rank_zero=True, threads=None, pure=False, ell_override=0, quadrature_boost=0, refinement_steps=0, timeout=604800, cancel=None, reuse_basis=True, execution_context=None):
    setup_start=time.monotonic()
    if type(ell_override) is not int or not 0<=ell_override<=4 or type(quadrature_boost) is not int or not 0<=quadrature_boost<=8:
        raise ValueError("invalid diagnostic override")
    if type(refinement_steps) is not int or not 0<=refinement_steps<=3: raise ValueError("invalid refinement steps")
    if type(reuse_basis) is not bool:raise ValueError("reuse_basis must be boolean")
    if timeout<=0: raise ValueError("invalid audit timeout")
    output, executable = Path(output), Path(executable)
    context = execution_context or AuditExecutionContext(executable)
    context.verify(executable)
    checkpoint, metadata = inspect(checkpoint, executable.with_name("alod_run"))
    if metadata["phase"] != 0:
        raise ValueError("only accepted checkpoints may be audited")
    if type(batch_size) is not int or not 1 <= batch_size <= 50:
        raise ValueError("audit batch_size must be in 1..50")
    config = {"wavenumber":16, **metadata["config"]}
    fixed = {key: config[key] for key in FIXED}
    fixed["maximum_nodes"] = min(fixed["maximum_nodes"], 200000)
    rows = validate(fixed)
    # The accepted checkpoint owns the coefficients; a later edit of the
    # repository's sample JSON cannot silently change a past experiment.
    rows = []
    for line in metadata['members_text'].splitlines():
        fields = line.split()
        if len(fields) != 7:
            raise ValueError('invalid frozen sample row')
        rows.append(dict(sample=int(fields[0]),split=fields[1],
                         **dict(zip(('wave_x','wave_y','singular_coefficient','wave_amplitude','wave_phase'),map(float,fields[2:])))))
    if member_ids is None:
        member_ids = list(config["member_ids"])
    if not isinstance(member_ids, list) or not member_ids or any(type(i) is not int for i in member_ids) or len(member_ids) != len(set(member_ids)):
        raise ValueError("audit member IDs must be distinct integers")
    if pure:
        if config["problem"] != "E2":
            raise ValueError("pure component audit is defined only for E2")
        nominal = dict(next(r for r in rows if r["sample"] == 0))
        rows += [dict(nominal, sample=48, split="pure", singular_coefficient=1, wave_amplitude=0),
                 dict(nominal, sample=49, split="pure", singular_coefficient=0, wave_amplitude=.5)]
        member_ids = member_ids + [i for i in (48, 49) if i not in member_ids]
    if not set(member_ids) <= {r["sample"] for r in rows}:
        raise ValueError("audit cannot invent members outside the frozen table")
    if threads is not None:
        if type(threads) is not int or not 0 <= threads <= 2147483647:
            raise ValueError("invalid audit threads")
        fixed["threads"] = threads
    fixed["maximum_nodes"] = config["maximum_nodes"]
    output.mkdir(parents=True, exist_ok=False)
    table = output/"members.txt"
    write_members(table, rows)
    manifest = dict(status="running", checkpoint=checkpoint.name, checkpoint_sha256=digest(checkpoint),
                    executable_sha256=context.executable_sha256, config=config, member_ids=member_ids,
                    ell_override=ell_override, quadrature_boost=quadrature_boost, refinement_steps=refinement_steps, batch_size=batch_size, reuse_basis=reuse_basis, fresh=fresh, rank_zero=rank_zero, threads=fixed["threads"],
                    solver_completed=False, validation_passed=False, audit_complete=False, paper_complete=False)
    manifest['provenance'] = context.provenance
    manifest['setup_seconds']=time.monotonic()-setup_start
    save = lambda: atomic_text(output/"run.json", json.dumps(manifest, indent=2)+"\n")
    save()
    start = time.monotonic()
    try:
        from execution import runtime_environment
        env = runtime_environment(fixed["threads"])
        env["ALOD_TIMING_FILE"]=str((output/"timings.jsonl").resolve())

        runtime = [f"--reuse-basis={int(reuse_basis)}", f"--enrichment-tests={config.get('enrichment_tests','adjoint')}", f"--refinement-steps={refinement_steps}", f"--ell-override={ell_override}", f"--quadrature-boost={quadrature_boost}", f"--checkpoint={checkpoint.resolve()}", f"--batch-size={batch_size}",
                   "--audit-ids="+",".join(map(str, member_ids)), f"--fresh={int(fresh)}", f"--rank-zero={int(rank_zero)}",
                   f"--method={config['method']}", f"--radius={config['radius']}", f"--rank-cap={config['rank_cap']}",
                   f"--maximum-patch-entries={config['maximum_patch_entries']}", f"--maximum-dense-entries={config['maximum_dense_entries']}"]
        with (output/"samples.jsonl").open("w") as out, (output/"stderr.log").open("w") as err:
            from execution import cancellable_run
            result = cancellable_run([str(executable.resolve()), *(["audit"] if executable.stem=="alod_run" else []), *arguments(fixed, table.resolve()), *runtime],
                                    stdout=out, stderr=err, env=env, timeout=timeout, cancel=cancel)
        context.verify(executable)
        manifest["returncode"] = result.returncode
        if result.returncode:
            raise RuntimeError("independent audit failed; see stderr.log")
        manifest["solver_completed"] = True
        records = [json.loads(line, parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
                   for line in (output/"samples.jsonl").read_text().splitlines()]
        samples = [r for r in records if r["kind"] == "sample"]
        if len(samples) != len(member_ids) or {r["sample"] for r in samples} != set(member_ids):
            raise ValueError("audit sample coverage mismatch")
        if any(r["state_id"] != metadata["state_id"] or r["PG_residual"] > 1e-8 for r in samples):
            raise ValueError("audit state or PG residual mismatch")
        if digest(checkpoint) != manifest["checkpoint_sha256"]:
            raise ValueError("audit modified its immutable checkpoint")
        if not records or records[-1].get("kind") != "audit_complete" or records[-1]['samples'] != len(samples):
            raise ValueError("audit completion record missing")
        for row in samples:
            if row['exact_norm'] > 1e-12:
                for raw, normalized in [('e','E'), ('f','F'), ('g','G')]:
                    if row[raw] is not None and abs(row[raw]/row['exact_norm']-row[normalized]) > 1e-12*max(1,abs(row[normalized])):
                        raise ValueError("audit normalization mismatch")
            if config['method']=='AFEM' and any(row[k] is not None for k in ('f','g','F','G','E_ref')):
                raise ValueError("AFEM reference fields must be not applicable")
        manifest.update(status="complete", validation_passed=True, audit_complete=True,
                        samples=len(samples), timing=records[-1], output_sha256=digest(output/"samples.jsonl"))
    except BaseException:
        manifest["status"] = "failed"
        raise
    finally:
        manifest["wall_seconds"] = time.monotonic()-start
        save()
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--checkpoint", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--executable", type=Path, default=ROOT/"build/alod_run")
    p.add_argument("--member-ids", help="comma-separated frozen member IDs")
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--threads", type=int)
    p.add_argument("--refinement-steps",type=int,default=0)
    p.add_argument("--ell-override",type=int,default=0)
    p.add_argument("--quadrature-boost",type=int,default=0)
    p.add_argument("--timeout",type=float,default=604800)
    p.add_argument("--fresh", action="store_true")
    p.add_argument("--pure", action="store_true")
    p.add_argument("--no-rank-zero", action="store_true")
    p.add_argument("--rebuild-lod", action="store_true", help="independently rebuild patches instead of restoring the accepted basis")
    a = p.parse_args()
    print(json.dumps(audit(a.checkpoint, a.output, a.executable, member_ids=None if a.member_ids is None else [int(i) for i in a.member_ids.split(',')], batch_size=a.batch_size, fresh=a.fresh, rank_zero=not a.no_rank_zero, threads=a.threads, pure=a.pure, ell_override=a.ell_override, quadrature_boost=a.quadrature_boost, refinement_steps=a.refinement_steps, timeout=a.timeout, reuse_basis=not a.rebuild_lod)))


if __name__ == "__main__":
    main()
