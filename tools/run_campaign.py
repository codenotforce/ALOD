"""Linux tmux supervisor: restartable jobs and memory guards, no CPU binding."""
import argparse,json,os,signal,subprocess,sys,time,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
GIB=1024**3
def memory_available():
 return int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024
def processes():
 records=[]
 for p in Path('/proc').iterdir():
  if not p.name.isdecimal():continue
  try:
   fields=(p/'stat').read_text().rsplit(')',1)[1].split()
   records.append((int(fields[2]),int(fields[21])*os.sysconf('SC_PAGE_SIZE'),int(fields[11])+int(fields[12])))
  except (OSError,ValueError,IndexError):pass
 return records
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True);p.add_argument('--build',type=Path,default=ROOT/'build');p.add_argument('--resume',action='store_true');a=p.parse_args()
 from execution import run_lease
 with run_lease(a.campaign/'campaign'):
  a.campaign=a.campaign.resolve();a.build=a.build.resolve();spec=json.loads((a.campaign/'campaign.json').read_text())
  status_path=a.campaign/'status.json'
  previous={}
  if status_path.exists():
   if not a.resume:raise RuntimeError('campaign already started; use --resume')
   previous=json.loads(status_path.read_text())
   for job in previous.get('active',[]):
    proc=Path('/proc')/str(job['pid'])/'cmdline'
    if proc.exists() and str(a.campaign).encode() in proc.read_bytes():
     raise RuntimeError('a previous campaign worker is still active; refusing duplicate launch')
  elif a.resume:raise RuntimeError('no campaign status to resume')
  slots=list(range(spec['maximum_active']))
  done=[j for j in previous.get('finished',[]) if j.get('status')=='complete']
  completed={j['name'] for j in done}
  pending=[j for j in spec['jobs'] if j['name'] not in completed];active=[];start=time.time();peak={};last={};logs=a.campaign/'logs';logs.mkdir(exist_ok=True);(a.campaign/'runs').mkdir(exist_ok=True)
  from execution import runtime_environment
  env=runtime_environment(supervised=False)
  def save(reason='running'):
   state=dict(status=reason,elapsed_seconds=time.time()-start,worker_slots=slots,available_gib=memory_available()/GIB,
    pending=[j['name'] for j in pending],active=[{k:v for k,v in x.items() if k not in ('process','handle')} for x in active],finished=done)
   temp=status_path.with_suffix('.tmp');temp.write_text(json.dumps(state,indent=2)+'\n');temp.replace(status_path)
  while pending or active:
   records=processes();now=time.time()
   for x in list(active):
    pid=x['pid'];rss=sum(r for g,r,t in records if g==pid);ticks=sum(t for g,r,t in records if g==pid)
    x['rss_gib']=rss/GIB;x['peak_rss_gib']=max(x.get('peak_rss_gib',0),rss/GIB)
    if pid in last:x['cpu_percent']=100*(ticks-last[pid][0])/os.sysconf('SC_CLK_TCK')/max(.01,now-last[pid][1])
    last[pid]=(ticks,now)
    if x['process'].poll() is None and (rss>x['memory_gib']*GIB or memory_available()<spec['emergency_available_gib']*GIB or shutil.disk_usage(a.campaign).free<spec['disk_reserve_gib']*GIB):
     x['termination_reason']='resource_guard';os.killpg(pid,signal.SIGTERM)
    rc=x['process'].poll()
    if rc is not None:
     x['handle'].close();done.append({**{k:v for k,v in x.items() if k not in ('process','handle')},'exit_code':rc,'status':'complete' if rc==0 else x.get('termination_reason','failed'),'wall_seconds':now-x['started_at']});active.remove(x)
   free_slots=[i for i in range(len(slots)) if i not in {x['slot'] for x in active}]
   committed=sum(max(0,x['memory_gib']-x.get('rss_gib',0)) for x in active)
   for slot in free_slots:
    index=next((i for i,j in enumerate(pending) if memory_available()/GIB-committed-j['memory_gib']>=spec['reserve_gib']),None)
    if index is None or shutil.disk_usage(a.campaign).free<spec['disk_reserve_gib']*GIB:break
    j=pending.pop(index);cmd=[sys.executable,str(ROOT/'tools/campaign_job.py'),'--campaign',str(a.campaign),'--name',j['name'],'--build',str(a.build)]
    handle=(logs/f'{j["name"]}.log').open('a');process=subprocess.Popen(cmd,stdout=handle,stderr=subprocess.STDOUT,env=env,start_new_session=True)
    x=dict(name=j['name'],slot=slot,cpu_binding=None,pid=process.pid,process=process,handle=handle,memory_gib=j['memory_gib'],started_at=now,command=cmd)
    active.append(x);committed+=j['memory_gib'];print(json.dumps(dict(event='started',name=j['name'],pid=process.pid,cpu_binding=None)),flush=True)
   save()
   with (a.campaign/'resources.jsonl').open('a') as f:f.write(json.dumps(dict(time=now,available_gib=memory_available()/GIB,active=[dict(name=x['name'],rss_gib=x.get('rss_gib',0),cpu_percent=x.get('cpu_percent',0)) for x in active]))+'\n')
   if active or pending:time.sleep(5)
  save('complete' if all(x['exit_code']==0 for x in done) else 'finished_with_failures')
  if any(x['exit_code']!=0 for x in done):raise SystemExit(1)
if __name__=='__main__':main()
