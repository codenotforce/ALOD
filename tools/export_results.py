"""Portable canonical tables with bounded CSV fields and hashed marking sets."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from checkpoint_io import atomic_text

CSV_LIMIT = 64*1024*1024


def legacy_rows(path):
    csv.field_size_limit(CSV_LIMIT)
    with Path(path).open(newline='', encoding='utf-8-sig') as stream:
        try:
            for row in csv.DictReader(stream):
                if None in row or any(v is None for v in row.values()):
                    raise ValueError("malformed legacy CSV row")
                yield row
        except csv.Error as error:
            raise ValueError("legacy CSV field exceeds the 64 MiB limit or is malformed") from error


def write_csv(path, rows):
    rows = list(rows)
    keys = list(dict.fromkeys(k for row in rows for k in row))
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, keys)
    writer.writeheader()
    writer.writerows(rows)
    atomic_text(path, stream.getvalue())


def export(run):
    run = Path(run)
    indices = run/'indices'
    indices.mkdir(exist_ok=True)
    def index(values):
        encoded = json.dumps(dict(format=1, element_ids=values), separators=(',', ':'))+'\n'
        key = hashlib.sha256(encoded.encode()).hexdigest()
        file = indices/(key+'.json')
        if not file.exists():
            atomic_text(file, encoded)
        elif hashlib.sha256(file.read_bytes()).hexdigest() != key:
            raise ValueError('existing marking set failed its content checksum')
        return key
    states, markings, samples, checks, timings = [], [], [], [], []
    with (run/'events.jsonl').open() as stream:
        for line in stream:
            row = json.loads(line)
            if row['kind'] == 'ell_check':
                row['event_id'] = f"ell-{row['ell_event_id']}"
                row['parent_id'] = f"state-{row['state_id']}"
                checks.append(row)
                continue
            if row['kind'] != 'accepted':
                continue
            state = {k:v for k,v in row.items() if not isinstance(v, (dict, list))}
            state['event_id'] = f"state-{row['state_id']}"
            for key, value in row.items():
                if key.endswith('_seconds'):
                    timings.append(dict(event_id=state['event_id']+'-'+key, parent_id=state['event_id'], phase=key[:-8], seconds=value, inclusive=key=='wall_seconds'))
            state['N_H'] = row['coarse_free']
            state['N_reference'] = row['reference_free'] if row['ell'] is not None else None
            state['N_work'] = row['coarse_free']+row.get('working_rank',row['rank'])
            state['N_online'] = row['coarse_free']+row['rank']
            for family in ('coarse', 'reference'):
                marks = row[family+'_marks']
                key = index(marks)
                state[family+'_mark_count'] = len(marks)
                state[family+'_mark_set'] = key
                markings.append(dict(state_id=row['state_id'], mesh=family, count=len(marks), set_id=key, sha256=key,
                                     mean_total=row.get(family+'_mean_total'),worst_total=row.get(family+'_worst_total'),
                                     mean_captured=row.get(family+'_mean_captured'),worst_captured=row.get(family+'_worst_captured'),
                                     mean_bulk=row.get(family+'_mean_bulk'), worst_bulk=row.get(family+'_worst_bulk'),
                                     worst_member_id=row.get(family+'_worst_member_id',row['worst_member_id']), frozen_scales=json.dumps(row['frozen_scales']),
                                     training_ids=json.dumps(row['training_ids'])))
            states.append(state)
            for sample in row.get('audit', []):
                samples.append(dict(state_id=row['state_id'], **{k:v for k,v in sample.items() if k != 'state_id'}))
    for name, rows in [('states',states), ('marking',markings), ('samples',samples), ('ell_checks',checks), ('timing',timings)]:
        write_csv(run/(name+'.csv'), rows)
    return dict(states=len(states), samples=len(samples), sets=len(list(indices.glob('*.json'))))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    print(json.dumps(export(p.parse_args().run)))
