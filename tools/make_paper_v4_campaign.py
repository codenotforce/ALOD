"""Prepare seven v4 adaptive jobs; this command never launches experiments."""
import argparse, json
from pathlib import Path
from run_adaptive import validate
ROOT=Path(__file__).resolve().parents[1]
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    names=['e1_main','e2_main']+[f'e3_k{k}' for k in (8,16,32,64,128)]
    configs={name:json.loads((ROOT/f'configs/adaptive/{name}.json').read_text()) for name in names}
    for c in configs.values():validate(c)
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'configs').mkdir()
    jobs=[]
    for name,c in configs.items():
        (args.output/f'configs/{name}.json').write_text(json.dumps(c,indent=2)+'\n',encoding='utf-8')
        jobs.append(dict(name=name,config=f'configs/{name}.json',kind='adaptive',memory_gib=160,threads=0))
    spec=dict(schema=1,jobs=jobs,maximum_active=1,threads_per_job=0,
              reserve_gib=80,emergency_available_gib=64,disk_reserve_gib=35,
              ell_diagnostic_kappa=[],scope='V4 bounded-horizon adaptive trajectories; independent audits; no inferred target completion.')
    (args.output/'campaign.json').write_text(json.dumps(spec,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(prepared=len(jobs),launched=False)))
if __name__=='__main__':main()
