"""Execute one isolated campaign job and its same-checkpoint diagnostics."""
import argparse,json,subprocess,sys
from pathlib import Path
from run_adaptive import run
from run_audit import audit
from checkpoint_io import digest
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--campaign',type=Path,required=True);p.add_argument('--name',required=True);p.add_argument('--build',type=Path,default=ROOT/'build');a=p.parse_args()
 definition=json.loads((a.campaign/'campaign.json').read_text());j=next(j for j in definition['jobs'] if j['name']==a.name)
 config_path=a.campaign/j['config'];config=json.loads(config_path.read_text());folder=a.campaign/'runs'/a.name
 if j['kind']=='baseline':
  subprocess.run([sys.executable,str(ROOT/'tools/run_baseline.py'),'--config',str(config_path),'--output',str(folder),'--executable',str(a.build/'alod_run')],check=True);return
 result=run(config,folder,a.build/'alod_run',timeout=30*86400,resume='auto' if folder.exists() else None)
 if not result['solver_completed'] or not result['audit_complete']:raise RuntimeError('trajectory or independent audit incomplete')
 if config['problem']=='E1' and config['wavenumber'] in definition['ell_diagnostic_kappa']:
  states=[json.loads(l) for l in (folder/'solver.jsonl').read_text().splitlines() if json.loads(l)['kind']=='accepted']
  samples={}
  for f in (folder/'audits').glob('*/samples.jsonl'):
   for row in map(json.loads,f.read_text().splitlines()):
    if row['kind']=='sample' and row['sample']==0:samples[row['state_id']]=row
  near=next((s['state_id'] for s in states if samples[s['state_id']]['E']<=.02),states[-1]['state_id'])
  ids=sorted({0,states[len(states)//2]['state_id'],near})
  diagnostics=[]
  for state in ids:
   cp=next((folder/'checkpoints').glob(f'state-{state:06d}-phase-0-ell-*.bin'))
   for label,kwargs in [('ell4',dict(ell_override=4)),('quadrature',dict(quadrature_boost=4,refinement_steps=2))]:
    out=folder/'diagnostics'/f'state-{state:06d}-{label}'
    audit(cp,out,a.build/'alod_run',member_ids=[0],batch_size=1,threads=config['threads'],**kwargs)
    diagnostics.append(dict(state=state,variant=label,checkpoint_sha256=digest(cp)))
  (folder/'diagnostics/index.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
 print(json.dumps(dict(job=a.name,status='complete')))
if __name__=='__main__':main()
