"""Validated checkpoint metadata, atomic files and accepted-journal recovery."""
import hashlib
import json
import os
from pathlib import Path
import subprocess


def atomic_text(path, text):
    path = Path(path)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(text)
        stream.flush()
        os.fsync(stream.fileno())
    temp.replace(path)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def resolve_checkpoint(path):
    path = Path(path)
    if path.is_dir():
        path = path/"checkpoints/latest" if (path/"checkpoints").is_dir() else path/"latest"
    if path.name == "latest":
        name = path.read_text().strip()
        if Path(name).name != name or not name.endswith(".bin"):
            raise ValueError("invalid latest checkpoint pointer")
        path = path.parent/name
    return path


def inspect(path, executable):
    path = resolve_checkpoint(path)
    result = subprocess.run([str(Path(executable).resolve()), "inspect-checkpoint", str(path.resolve())],
                            capture_output=True, text=True, timeout=120)
    if result.returncode:
        raise ValueError("checkpoint validation failed: " + result.stderr.strip())
    return path, json.loads(result.stdout)


def recover_journal(path, metadata):
    """Keep only committed rows, recovering the last accepted transaction once."""
    path = Path(path)
    limit = metadata["committed_lines"]
    if type(limit) is not int or limit < 0:
        raise ValueError("invalid checkpoint log cursor")
    temp = path.with_name(path.name + ".recover")
    count = 0
    hash_value = 14695981039346656037
    with temp.open("w", encoding="utf-8", newline="\n") as out:
        if path.exists():
            with path.open() as stream:
                for line in stream:
                    if count == limit:
                        break
                    json.loads(line)
                    for byte in line.encode('utf-8'):
                        hash_value = ((hash_value ^ byte) * 1099511628211) & ((1 << 64)-1)
                    out.write(line)
                    count += 1
        if count != limit:
            raise ValueError("checkpoint requires a missing committed journal prefix")
        if str(hash_value) != metadata['journal_hash']:
            raise ValueError("committed journal prefix does not match checkpoint")
        if metadata["phase"] == 0:
            for line in metadata["journal"].splitlines():
                json.loads(line)
                out.write(line + "\n")
        out.flush()
        os.fsync(out.fileno())
    temp.replace(path)


def validate_resume_config(old, new):
    increases = {"cycles", "state_limit", "maximum_nodes", "maximum_patch_entries", "maximum_dense_entries"}
    operational = {"audit", "emit_solution", "threads"}
    for key in old:
        if key in operational:
            continue
        if key in increases:
            if key == "state_limit" and (new[key] == 0 or (old[key] == 0 and new[key] == 0)):
                continue
            if key == "state_limit" and old[key] == 0 and new[key] != 0:
                raise ValueError("resume cannot reduce an unlimited state horizon")
            if new[key] < old[key]:
                raise ValueError(f"resume may only increase {key}")
        elif old[key] != new[key]:
            raise ValueError(f"resume cannot change mathematical field {key}; create an explicit new experiment")
