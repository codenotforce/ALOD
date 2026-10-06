"""Reclaim audited state inputs without changing their numerical results.

Only the supervisor mutates state paths. Geometry collection is deferred until
the producer stops, because its background writer can publish geometry before
the corresponding state file. Recovery evidence is never removed.
"""
import hashlib
import json
import os
import re
import struct
from pathlib import Path

from checkpoint_io import atomic_text, sync_directory

STATE_NAME = re.compile(r"state-(\d+)-phase-([012])-ell-(\d+)\.bin\Z")
GEOMETRY_NAME = re.compile(r"mesh-[0-9a-f]{16}\.bin\Z")
MAGICS = {0x3154504B43444F4C, 0x3254504B43444F4C, 0x3354504B43444F4C}


def accepted_hash(row):
    # events.jsonl adds audit results to the original accepted row.
    row = {key: value for key, value in row.items() if key != 'audit'}
    return hashlib.sha256(json.dumps(row, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def state_header(path):
    """Read bounded metadata only; audit already validates the full payload."""
    with Path(path).open('rb') as stream:
        def word():
            data = stream.read(8)
            if len(data) != 8:
                raise ValueError('truncated state header')
            return struct.unpack('<Q', data)[0]

        def text():
            size = word()
            # Match the native Reader's string bound and reject truncation
            # before allocation, without imposing a new experiment-size cap.
            if size > min(2**31 - 1, os.fstat(stream.fileno()).st_size - stream.tell()):
                raise ValueError('state metadata length outside file bounds')
            data = stream.read(size)
            if len(data) != size:
                raise ValueError('truncated state metadata')
            return data.decode('utf-8')

        magic = word()
        if magic not in MAGICS:
            raise ValueError('unsupported state format')
        metadata = json.loads(text())
        geometry = text() if magic == 0x3354504B43444F4C else ''
        if geometry and not GEOMETRY_NAME.fullmatch(geometry):
            raise ValueError('invalid state geometry reference')
        return metadata, geometry


class StateRetention:
    def __init__(self, output, config):
        self.output = Path(output).resolve()
        self.config = config or {}
        self.mode = self.config.get('checkpoint_retention', 'all')
        self.accepted = {}
        self.offset = 0
        self.removed = 0
        self.removed_bytes = 0
        self.geometry_deferred = None

    def update_journal(self):
        """Consume only complete, committed journal rows, once per invocation."""
        journal = self.output / 'solver.jsonl'
        if not journal.exists():
            return
        with journal.open('rb') as stream:
            stream.seek(self.offset)
            while True:
                line = stream.readline()
                if not line or not line.endswith(b'\n'):
                    break
                row = json.loads(line)
                self.offset = stream.tell()
                if row.get('kind') == 'accepted':
                    sid = row['state_id']
                    self.accepted[sid] = dict(
                        accepted_hash=accepted_hash(row), coarse_cycle=row['coarse_cycle'],
                        cycle_complete=row['cycle_complete'], target_reached=row.get('target_reached', False))

    def paths(self, name):
        if not STATE_NAME.fullmatch(name):
            raise ValueError('invalid retained state name')
        return [self.output / folder / name for folder in ('checkpoints', 'audit_snapshots')]

    def keep(self, state):
        horizon = 1 + self.config.get('cycles', 0) * self.config.get('m_ref', 2)
        cap = self.config.get('state_limit', 0)
        final = state['target_reached'] or state['state_id'] == (min(horizon, cap) if cap else horizon) - 1
        if self.mode == 'all':
            return True
        if final and (self.mode == 'periodic' or self.config.get('keep_final', False)):
            return True
        interval = self.config.get('checkpoint_interval_cycles', 0)
        return (self.mode == 'periodic' and state['coarse_cycle'] > 0
                and state['cycle_complete'] and state['coarse_cycle'] % interval == 0)

    def reclaim(self, jobs, futures, *, producer_done=False, protected=()):
        if self.mode == 'all':
            return
        self.update_journal()
        latest = self.output / 'checkpoints/latest'
        latest_name = latest.read_text().strip() if latest.exists() else None
        protected = set(protected)
        if latest_name and not producer_done:
            protected.add(latest_name)
        remove = []
        for name, job in jobs.items():
            match = STATE_NAME.fullmatch(name)
            if not match or job.get('status') != 'complete':
                continue
            future = futures.get(name)
            if future is None or not future.done() or future.exception() is not None:
                continue
            sid = int(match[1])
            committed = self.accepted.get(sid)
            if not committed or committed['accepted_hash'] != job.get('accepted_hash'):
                continue
            if not self.keep(dict(committed, state_id=sid)) and name not in protected:
                remove.extend(self.paths(name))
            # Training/ell restart boundaries are superseded only after the
            # accepted state of that same step has been audited successfully.
            for folder in ('checkpoints', 'audit_snapshots'):
                for path in (self.output / folder).glob(f'state-{sid:06d}-phase-*-ell-*.bin'):
                    parsed = STATE_NAME.fullmatch(path.name)
                    if parsed and parsed[2] != '0' and path.name not in protected:
                        remove.append(path)
        remove = list(dict.fromkeys(path for path in remove if path.exists()))
        for path in remove:
            if path.is_symlink() or path.parent.resolve() not in (
                    self.output / 'checkpoints', self.output / 'audit_snapshots'):
                raise ValueError('refusing to reclaim a state outside the run directory')
        if remove:
            # Observing a complete row is insufficient for power-loss recovery:
            # synchronize its journal prefix before discarding numerical inputs.
            with (self.output/'solver.jsonl').open('rb') as stream:
                os.fsync(stream.fileno())
        # Once the writer has stopped, publish the surviving restart pointer
        # before unlinking its former target. A crash never leaves a stale latest.
        if producer_done and latest_name and any(path.name == latest_name for path in remove):
            survivors = [p for p in (self.output / 'checkpoints').glob('state-*-phase-*-ell-*.bin')
                         if STATE_NAME.fullmatch(p.name) and p not in remove]
            if survivors:
                def order(path):
                    parsed = STATE_NAME.fullmatch(path.name)
                    return int(parsed[1]), int(parsed[2] == '0'), int(parsed[3]), int(parsed[2])
                atomic_text(latest, max(survivors, key=order).name + '\n')
            else:
                latest.unlink()
                sync_directory(latest.parent)
        for path in remove:
            self.removed_bytes += path.stat().st_size
            path.unlink()
            self.removed += 1
        for directory in {path.parent for path in remove}:
            sync_directory(directory)
        if producer_done:
            self.collect_geometry()

    def collect_geometry(self):
        # Include quarantined states. On unreadable metadata, fail closed and
        # leave all geometry untouched rather than guessing its dependencies.
        references = set()
        for folder in ('checkpoints', 'audit_snapshots'):
            directory = self.output / folder
            for path in directory.rglob('state-*-phase-*-ell-*.bin'):
                if not STATE_NAME.fullmatch(path.name):
                    continue
                if path.is_symlink():
                    raise ValueError('symlinked state prevents geometry collection')
                try:
                    _, geometry = state_header(path)
                except (OSError, ValueError) as error:
                    self.geometry_deferred = str(error)
                    return
                if geometry:
                    references.add(geometry)
        for folder in ('checkpoints', 'audit_snapshots'):
            geometry_dir = self.output / folder / 'meshes'
            for path in geometry_dir.glob('mesh-*.bin'):
                if GEOMETRY_NAME.fullmatch(path.name) and path.name not in references:
                    if path.is_symlink() or path.parent.resolve() != geometry_dir:
                        raise ValueError('refusing to reclaim geometry outside the run directory')
                    self.removed_bytes += path.stat().st_size
                    path.unlink()
                    self.removed += 1
            if geometry_dir.exists():
                sync_directory(geometry_dir)

    def summary(self):
        return dict(mode=self.mode, keep_final=self.config.get('keep_final', False),
                    removed_paths=self.removed, removed_logical_bytes=self.removed_bytes,
                    geometry_collection_deferred=self.geometry_deferred)
