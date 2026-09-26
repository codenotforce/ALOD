"""Validate an explicit P1 baseline preset and write a portable run artifact."""
import argparse
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
FIELDS={'reference_gap','wavenumber','schema_version','problem','method','initial_level','states','theta','target','maximum_nodes','threads','emit_solution'}


def validate(config):
    if isinstance(config,dict): config={"reference_gap":4,"wavenumber":16, **config}
    if not isinstance(config,dict) or config.keys()!=FIELDS:raise ValueError('preset must contain exactly the documented schema fields')
    if type(config['schema_version']) is not int or config['schema_version']!=1:raise ValueError('unsupported schema_version')
    if config['problem'] not in ('E1','E2') or config['method'] not in ('AFEM','UFEM','SLOD'):raise ValueError('unknown problem or baseline')
    for name,lo,hi in [('reference_gap',1,8),('initial_level',0,24),('states',1,10000),('maximum_nodes',4,2_000_000_000),('threads',0,2147483647)]:
        x=config[name]
        if type(x) is not int or not lo<=x<=hi:raise ValueError(f'invalid {name}')
    for name in ['theta','target']:
        x=config[name]
        if type(x) not in (float,int) or not math.isfinite(x):raise ValueError(f'{name} must be finite')
    if not 0<config['theta']<=1 or config['target']<0:raise ValueError('invalid theta or target')
    if type(config['emit_solution']) is not bool:raise ValueError('emit_solution must be boolean')
    if type(config["wavenumber"]) not in (float,int) or not math.isfinite(config["wavenumber"]) or config["wavenumber"] <= 0: raise ValueError("invalid wavenumber")
    return config


def command(config):
    c=validate(config)
    args=[c['problem'],c['method']]
    for key in ['reference_gap','wavenumber','initial_level','states','theta','target','maximum_nodes','threads','emit_solution']:
        value=int(c[key]) if key=='emit_solution' else c[key]
        args.append(f'--{key.replace("_","-")}={value}')
    return args


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--executable',type=Path,default=ROOT/'build/alod_run')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    c=validate(json.loads(a.config.read_text(encoding='utf-8')))
    exe=a.executable.resolve()
    if not exe.is_file():raise ValueError('baseline executable is missing; build the release preset first')
    a.output.mkdir(parents=True,exist_ok=False)
    canonical=json.dumps(c,sort_keys=True,separators=(',',':')).encode()
    record={'config':c,'config_sha256':hashlib.sha256(canonical).hexdigest(),
        'executable_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),
        'arguments':command(c),'BLAS_threads':1,'status':'running'}
    manifest=a.output/'run.json'
    def save():manifest.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    save()
    from execution import runtime_environment
    env = runtime_environment(c["threads"])
    # Stream potentially large state output. Paths and host details are not serialized.
    with (a.output/'states.jsonl').open('w',encoding='utf-8') as output, (a.output/'stderr.log').open('w',encoding='utf-8') as errors:
        result=subprocess.run([str(exe),*command(c)],stdout=output,stderr=errors,env=env)
    record['status']='completed' if result.returncode==0 else 'failed'
    record['exit_code']=result.returncode
    count=0;last=None
    with (a.output/'states.jsonl').open(encoding='utf-8') as f:
        for line in f: last=json.loads(line);count+=1
    record['states_written']=count
    record['stop_reason']=last['stop_reason'] if last else None
    save()
    if result.returncode:raise SystemExit('Baseline failed; inspect its stderr.log and run.json.')
    print(f'{c["problem"]} {c["method"]}: {count} states; {record["stop_reason"]}')


if __name__=='__main__':main()
