"""Create the predeclared E1 kappa controls and revised E1/E2 production jobs."""
import argparse,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--build',type=Path,default=ROOT/'build');a=p.parse_args()
 a.output.mkdir(parents=True,exist_ok=False);(a.output/'configs').mkdir()
 geometry=[json.loads(x) for x in subprocess.check_output([str(a.build/'alod_geometry')],text=True).splitlines()]
 (a.output/'initial_geometry.json').write_text(json.dumps(geometry,indent=2)+'\n')
 jobs=[]
 def add(name,c,kind,budget):
  (a.output/'configs'/f'{name}.json').write_text(json.dumps(c,indent=2)+'\n')
  jobs.append(dict(name=name,config=f'configs/{name}.json',kind=kind,memory_gib=budget,threads=0))
 # These are new mathematical trajectories, never archive-reproduction labels.
 for problem,filename in [('E2','e2_main'),('E1','e1_main')]:
  c=json.loads((ROOT/'configs/adaptive'/f'{filename}.json').read_text())
  c.update(wavenumber=16,threads=0,maximum_nodes=1000000,maximum_patch_entries=16000000000,maximum_dense_entries=64000000)
  add(f'{problem.lower()}-paper-lazy',c,'adaptive',112)
 base=json.loads((ROOT/'configs/adaptive/e1_main.json').read_text())
 for g in geometry:
  k=int(g['wavenumber']);level=g['level']
  c=dict(base,wavenumber=k,level=level,gap=4,ell=3,maximum_ell=4,ell_mode='fixed',
         extra_checks=[],force_promotions=[],member_ids=[0],training_ids=[0],threads=0,
         maximum_nodes=1000000,maximum_patch_entries=16000000000,maximum_dense_entries=64000000,
         cycles=25,m_ref=2,reference_theta=[.3,.2],minimum_gap=2,audit=True)
  add(f'k{k}-alod',c,'adaptive',112)
  for method in ['AFEM','UFEM','SLOD']:
   # Uniform trajectories use a predeclared horizon up to coarse level 15
   # (SLOD reference level 19), FEM level 20. Target E<=.01 stops baselines.
   states=97 if method=='AFEM' else max(1,(15 if method=='SLOD' else 20)-level+1)
   c=dict(schema_version=1,problem='E1',method=method,wavenumber=k,initial_level=level,
          states=states,theta=.15,target=.01,maximum_nodes=2000000,threads=0,emit_solution=False)
   add(f'k{k}-{method.lower()}',c,'baseline',48 if method=='SLOD' else 32)
 # Interleave inexpensive controls with the two revised production trajectories.
 names=['e2-paper-lazy','k8-alod','k8-afem','k8-ufem','k8-slod','e1-paper-lazy']
 names += [j['name'] for j in jobs if j['name'] not in names]
 jobs=sorted(jobs,key=lambda j:names.index(j['name']))
 campaign=dict(schema=1,jobs=jobs,maximum_active=3,threads_per_job=0,
               reserve_gib=80,emergency_available_gib=64,disk_reserve_gib=35,
               kappa=[8,16,32,64],targets=[.02,.01],ell_diagnostic_kappa=[32,64],
               diagnostic_states='initial; middle accepted state; first E<=.02, otherwise last',
               mathematical_policies=dict(interpolation='area',riesz_patches='n2',e2_bootstrap='empty_worst_budget',e2_ell_mode='lazy',e2_threshold=.1),
               scope='18 new trajectories, with offline audits; no manuscript edits')
 (a.output/'campaign.json').write_text(json.dumps(campaign,indent=2)+'\n')
 print(json.dumps(dict(jobs=len(jobs),geometry=geometry)))
if __name__=='__main__':main()
