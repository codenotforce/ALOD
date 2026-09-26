"""Sample an existing supervisor and its descendants without changing affinity."""
import argparse
import json
import os
import time
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--pid', type=int, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--interval', type=float, default=5)
args = parser.parse_args()
if args.interval <= 0:
    parser.error('interval must be positive')
proc = Path(os.sep) / 'proc'
previous = {}
start_identity = None
with args.output.open('a', encoding='utf-8') as output:
    while True:
        processes = {}
        for path in proc.iterdir():
            if not path.name.isdigit():
                continue
            try:
                fields = (path / 'stat').read_text().rsplit(')', 1)[1].split()
                processes[int(path.name)] = (path, fields)
            except (OSError, IndexError):
                pass
        if args.pid not in processes:
            break
        identity = processes[args.pid][1][19]
        if start_identity is None:
            start_identity = identity
        if identity != start_identity:
            break
        selected = {args.pid}
        while True:
            children = {pid for pid, (_, fields) in processes.items() if int(fields[1]) in selected}
            if children <= selected:
                break
            selected |= children
        rows, current = [], {}
        for pid in sorted(selected):
            try:
                path, fields = processes[pid]
                status = dict(line.split(':', 1) for line in (path / 'status').read_text().splitlines())
                active = 0
                for thread in (path / 'task').iterdir():
                    stat = (thread / 'stat').read_text().rsplit(')', 1)[1].split()
                    ticks = int(stat[11]) + int(stat[12])
                    key = (pid, thread.name, stat[19])
                    active += key in previous and ticks > previous[key]
                    current[key] = ticks
                rows.append(dict(pid=pid, parent=int(fields[1]), command=status['Name'].strip(),
                                 threads=int(status['Threads']), active_threads=active,
                                 cpu_ticks=int(fields[11])+int(fields[12]),
                                 allowed_cpus=status['Cpus_allowed_list'].strip(),
                                 rss_kib=int(status.get('VmRSS', '0').split()[0]),
                                 swap_kib=int(status.get('VmSwap', '0').split()[0])))
            except (OSError, KeyError, IndexError):
                pass
        output.write(json.dumps(dict(unix_seconds=time.time(), processes=rows))+'\n')
        output.flush()
        previous = current
        time.sleep(args.interval)
