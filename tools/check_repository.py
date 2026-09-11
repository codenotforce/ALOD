"""Check portable text, frozen source hashes and canonical audit completeness."""
import argparse
import csv
import hashlib
import io
import json
import math
import re
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CSV_FIELD_LIMIT = 64 * 1024 * 1024
PATH_PATTERNS = [
    re.compile(r'(?i)(?<![a-z0-9])[a-z]:[\\/]'),
    re.compile(r'(?<!\\)\\\\[a-zA-Z0-9_.-]+\\'),
    re.compile(r'(?<![\w>/])/(?:home|Users|mnt|media|tmp|var|opt|usr|lib|lib64|etc|workspaces)/[^\s"`<>]*'),
]


def path_violations(text):
    # URLs are not local paths; no exception is made for Markdown code blocks.
    text = re.sub(r'https?://[^\s<>"`]+', '<URL>', text)
    return [i for i, line in enumerate(text.splitlines(), 1) if any(p.search(line) for p in PATH_PATTERNS)]


def canonical_bytes(data):
    return data.replace(b'\r\n', b'\n')


def read_audit(text, states, *, e2=False, marking=None, parameters=None):
    csv.field_size_limit(CSV_FIELD_LIMIT)
    rows = csv.DictReader(io.StringIO(text, newline=''))
    counts = Counter()
    identities = set()
    state_key = 'H_step' if e2 else 'state'
    energy_keys = ['exact_norm', 'exact_error', 'floor', 'gap', 'PG_residual'] if e2 else ['exact_norm', 'e', 'E', 'PG_residual']
    for row in rows:
        if None in row or any(value is None for value in row.values()):
            raise ValueError('malformed audit row')
        state, sample = int(row[state_key]), int(row['sample'])
        key = (state, sample)
        if key in identities:
            raise ValueError(f'duplicate state/sample {key}')
        identities.add(key)
        counts[state] += 1
        if not 0 <= state < states or not 0 <= sample < (50 if e2 else 48):
            raise ValueError(f'out-of-range state/sample {key}')
        split = 'train' if sample < 16 else 'test' if sample < 40 else 'shift' if sample < 48 else 'pure'
        valid_role = row['split'] == split or (split == 'shift' and row['split'] == 'shift-audit')
        if not valid_role and not (e2 and sample >= 48 and row['split'].startswith('pure')):
            raise ValueError(f'incorrect role {key}')
        if marking is not None and int(row['marking_member']) != int(sample in marking):
            raise ValueError(f'incorrect marking identity {key}')
        if not e2 and parameters is not None:
            for field in ['wave_x', 'wave_y']:
                if abs(float(row[field])-parameters[sample][field]) > 1e-14:
                    raise ValueError(f'changed sample parameter {key}/{field}')
        for field in energy_keys:
            value = float(row[field])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f'invalid {field} at {key}')
        if float(row['PG_residual']) >= 1e-7:
            raise ValueError(f'PG residual at {key} exceeds archived tolerance')
    required_samples = 50 if e2 else 48
    if counts != Counter({state: required_samples for state in range(states)}):
        raise ValueError('missing state or sample in canonical audit')
    return sum(counts.values())


def git(*args):
    return subprocess.check_output(['git', '-c', 'core.quotepath=false', *args], cwd=ROOT)


def verify(staged=False):
    names = git('ls-files', '-z', '--cached') if staged else git('ls-files', '-z', '--cached', '--others', '--exclude-standard')
    files = sorted(set(x.decode('utf-8') for x in names.split(b'\0') if x))
    def read(name):
        return git('show', ':'+name) if staged else (ROOT/name).read_bytes()
    errors = []
    for name in files:
        try:
            content = read(name).decode('utf-8-sig')
        except UnicodeDecodeError:
            errors.append(f'{name}: unexpected binary file in initial source checkpoint')
            continue
        lines = path_violations(content)
        if lines:
            errors.append(f'{name}: absolute local path at lines {lines[:8]}')
        decoded = content
        if name.endswith('.json'):
            try:
                decoded = json.dumps(json.loads(content), ensure_ascii=False)
            except ValueError as error:
                errors.append(f'{name}: invalid JSON ({error})')
        if name.endswith(('.md', '.json')) and re.search(r'[\u4e00-\u9fff]', decoded):
            errors.append(f'{name}: documentation must be written in English')
    manifest = json.loads(read('docs/provenance/source_manifest.json'))
    for entry in manifest['files']:
        try:
            digest = hashlib.sha256(canonical_bytes(read(entry['path']))).hexdigest()
            if digest != entry['sha256']:
                errors.append(f"{entry['path']}: imported-content hash mismatch")
        except (OSError, subprocess.CalledProcessError):
            errors.append(f"{entry['path']}: missing imported file")
    sample_tables = {}
    for problem in ['e1', 'e2']:
        data = json.loads(read(f'data/rhs/rhs_{problem}.json'))
        if [r['sample'] for r in data] != list(range(48)):
            errors.append(f'{problem}: formal sample identities are not exactly 0..47')
        if Counter(r['split'] for r in data) != Counter(train=16, test=24, shift=8):
            errors.append(f'{problem}: incorrect sample roles')
        if [r['sample'] for r in data if r['nominal']] != [0]:
            errors.append(f'{problem}: nominal must be sample 0 only')
        sample_tables[problem] = data
    total = 0
    for arm, states in [('nominal-alod',51), ('family-alod',51), ('nominal-afem',97), ('family-afem',97)]:
        try:
            total += read_audit(read(f'tests/fixtures/legacy/e1-{arm}.csv').decode(), states,
                                marking={0} if arm.startswith('nominal') else set(range(16)),
                                parameters=sample_tables['e1'])
        except (ValueError, csv.Error) as error:
            errors.append(f'{arm}: {error}')
    try:
        e2_total = read_audit(read('tests/fixtures/legacy/e2-accepted.csv').decode(),33,e2=True)
    except (ValueError, csv.Error) as error:
        errors.append(f'E2: {error}')
        e2_total = 0
    if errors:
        raise ValueError('\n'.join(errors))
    print(f'Validated {len(files)} text files, {len(manifest["files"])} imported hashes, '
          f'{total} E1 and {e2_total} E2 audit rows; no local absolute paths.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--staged', action='store_true', help='Inspect exact staged Git blobs')
    args = parser.parse_args()
    try:
        verify(args.staged)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error))
