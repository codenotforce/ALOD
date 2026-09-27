"""Versioned, path-free v4 paper data import, coverage validation and figure delivery.

Input spec schema 1 maps named runs to native run directories or portable exports.
The archive command imports the supplied September v4 export layout explicitly.
Historical imports retain their provenance class and are never new-binary evidence.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import shutil
from pathlib import Path
import xml.etree.ElementTree as ET

METHODS=('ALOD','AFEM','UFEM','LOD')
KAPPAS=(8,16,32,64,128)
FIGURES=('e1_final_meshes','e1_error_dof','e2_final_meshes','e2_error_dof','kappa_error_dof')
TABLES=('adaptive-config','e1-terminal-results','e1-family-control','e2-terminal-results','kappa-target-summary')

def read(path):
    path=Path(path)
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes())

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def adaptive(data, *, complete=True):
    states=data['states'];audits=data['audits'];seen={}
    for row in audits:
        key=(row['state_id'],row['sample'])
        if key in seen:raise ValueError('duplicate state/member audit')
        seen[key]=row
    config=data.get('config',{})
    expected=config.get('member_ids',list(range(48)))
    curves=[]
    for state in states:
        sid=state['state_id']
        if {m for s,m in seen if s==sid}!=set(expected):raise ValueError('incomplete frozen audit coverage')
        audit=seen[sid,0]
        if float(audit.get('PG_residual',0))>1e-8:raise ValueError('audit residual gate failed')
        row=dict(state=sid,N_on=state['coarse_free']+state.get('rank',0),E=float(audit['E']),
                 ell=state.get('ell'),cycle_complete=state.get('cycle_complete',True),
                 rank=state.get('rank',0),eta_relative=state['eta'][0]/float(audit['exact_norm']))
        tests=[float(a['E']) for (s,m),a in seen.items() if s==sid and a.get('role')=='test']
        if tests:
            import numpy as np
            row.update(test_count=len(tests),test_p50=float(np.quantile(tests,.5)),
                       test_p90=float(np.quantile(tests,.9)),test_max=max(tests))
        curves.append(row)
    return dict(rows=curves,config=config,complete=complete,
                audits=[{k:a.get(k) for k in ('state_id','sample','role','E','F','G','exact_norm')} for a in audits])

def native(folder,executable=None):
    folder=Path(folder);manifest=read(folder/'run.json')
    if (folder/'events.jsonl').exists():
        events=[json.loads(line) for line in (folder/'events.jsonl').read_text().splitlines()]
        states=[r for r in events if r['kind']=='accepted']
        out=adaptive(dict(states=states,audits=[a for r in states for a in r.get('audit',[])],config=manifest['config']),
                     complete=bool(manifest.get('audit_complete') and manifest.get('solver_completed')))
    else:
        rows=[json.loads(line) for line in (folder/'states.jsonl').read_text().splitlines()]
        out=dict(rows=[dict(state=r['state'],N_on=r['free_dof'],E=r['relative_energy_error'],
                           ell=r.get('ell'),cycle_complete=True) for r in rows],config=manifest['config'],
                 complete=manifest['status']=='completed')
    out.update(origin='native_run',source_sha256=sha(folder/'run.json'),source_data_sha256=sha(folder/('events.jsonl' if (folder/'events.jsonl').exists() else 'states.jsonl')),
               executable_sha256=manifest.get('executable_sha256'))
    meshes=list(folder.glob('audits/*/mesh.json'))
    if (folder/'mesh.json').exists():meshes=[folder/'mesh.json']
    if meshes:out['mesh']=read(meshes[-1])
    elif executable is not None and (folder/'checkpoints/latest').exists():
        import subprocess,tempfile
        from checkpoint_io import resolve_checkpoint
        with tempfile.TemporaryDirectory() as temp:
            file=Path(temp)/'mesh.json'
            subprocess.run([str(Path(executable).resolve()),'checkpoint-mesh',str(resolve_checkpoint(folder/'checkpoints/latest').resolve()),str(file)],check=True)
            out['mesh']=read(file)
    return out

def mesh(path):
    path=Path(path)
    if '.vtu' not in path.name:
        x=read(path);return dict(nodes=x['nodes'],elements=x.get('elements',x.get('triangles')))
    root=ET.fromstring(gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes())
    import numpy as np
    points=root.find('.//Points/DataArray')
    if points.get('format','ascii')!='ascii':raise ValueError('VTU importer requires explicit ASCII arrays')
    nodes=np.fromstring(points.text,sep=' ').reshape(-1,3)[:,:2].tolist()
    cells={a.get('Name'):np.fromstring(a.text,sep=' ',dtype=int) for a in root.findall('.//Cells/DataArray')}
    if any(cells['types']!=5):raise ValueError('VTU importer requires triangular cells')
    return dict(nodes=nodes,elements=cells['connectivity'].reshape(-1,3).tolist())

def import_archives(e1,e2,e3):
    """Named adapter for the supplied q03/kernel-lift/k8--128 exported schemas."""
    e1,e2,e3=map(Path,(e1,e2,e3));runs={}
    for exp,source in [('e1',e1/'data'),('e2',e2/'data')]:
        file=source/f'{exp}-production.json';data=read(file)
        runs[f'{exp}-ALOD']=adaptive(data)
        runs[f'{exp}-ALOD'].update(origin='historical_v4_export',source_sha256=sha(file),
            mesh={key:mesh(source/f'{exp}-final-{suffix}.json.gz') for key,suffix in [('coarse','coarse'),('reference','reference')]})
        baselines=e1/'data'/f'{exp}-baselines.json'
        for method,record in read(baselines)['runs'].items():
            name='LOD' if method=='SLOD' else method
            runs[f'{exp}-{name}']=dict(display_state_limit=({'AFEM':77,'UFEM':13,'LOD':10}[name] if exp=='e2' else len(record['trajectory'])),rows=record['trajectory'],config={},complete=True,
                origin='historical_v4_export',source_sha256=sha(baselines))
        runs[f'{exp}-AFEM']['mesh']={'coarse':mesh(e1/'data'/f'{exp}-afem-final.vtu.gz')}
    controls=e1/'data/e1-controls.json'
    for name,data in read(controls).items():
        runs['control-'+name]=adaptive(data)
        runs['control-'+name].update(origin='historical_v4_export',source_sha256=sha(controls))
    control_csv=e1/'controls/data/afem/state_statistics.csv'
    with control_csv.open(newline='',encoding='utf-8-sig') as f:statistics=list(csv.DictReader(f))
    for arm in ('nominal-afem','family-afem'):
        rows=[dict(state=int(r['state']),N_on=int(r['N_online']),E=float(r['nominal_E']),
                   test_count=24,test_p50=float(r['test_E_p50']),test_p90=float(r['test_E_p90']),
                   test_max=float(r['test_E_max']),cycle_complete=True) for r in statistics if r['arm']==arm]
        runs['control-'+arm]=dict(rows=rows,config={},complete=True,origin='historical_summary_export',source_sha256=sha(control_csv))
    summary=read(e1/'controls/analysis.json')['e1']['budgets']
    for arm,budgets in summary.items():
        runs['control-'+arm]['budget_statistics']=budgets['2000']
    for r in read(e3/'data.json')['runs']:
        name='LOD' if r['method']=='SLOD' else r['method']
        runs[f'k{r["k"]}-{name}']=dict(rows=[dict(state=x.get('state_id',x.get('state')),N_on=x.get('coarse_free',x.get('free_dof')),
                    E=x['relative_energy_error'],ell=x.get('ell'),cycle_complete=x.get('cycle_complete',True)) for x in r['rows']],
                    config=r['config'],complete=r['status'] in ('complete','target_reached'),resource_stopped=r.get('guard_reason') is not None,
                    origin='historical_v4_export',source_sha256=sha(e3/'data.json'))
    return runs

def clean(runs):
    """Whitelist portable scientific data; original archives may contain host paths."""
    allowed=('problem','method','wavenumber','level','gap','initial_level','ell','maximum_ell','ell_threshold',
             'ell_absolute_threshold','ritz_tolerance','enrichment_tests','theta','reference_theta','m_ref','training_ids','member_ids')
    result={}
    for name,run in runs.items():
        if not name or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_' for c in name):raise ValueError('unsafe run name')
        rows=[]
        for index,row in enumerate(run['rows']):
            r={k:row[k] for k in ('state','N_on','E','eta_relative','ell','rank','cycle_complete','test_count','test_p50','test_p90','test_max') if k in row}
            r.setdefault('state',index)
            if type(r['N_on']) not in (int,float) or not math.isfinite(r['N_on']) or r['N_on']<=0:raise ValueError('invalid online dimension')
            if type(r['E']) not in (int,float) or not math.isfinite(r['E']) or r['E']<=0:raise ValueError('invalid error')
            rows.append(r)
        if not rows or [r['state'] for r in rows]!=list(range(len(rows))):raise ValueError('trajectory has gaps or duplicates')
        result[name]={k:run[k] for k in ('origin','source_sha256','source_data_sha256','executable_sha256','complete','resource_stopped','mesh','audits','display_state_limit','budget_statistics') if k in run}
        result[name].update(rows=rows,config={k:v for k,v in run.get('config',{}).items() if k in allowed})
    return result

def deliver(runs, output, *, allow_partial=False):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.patches import Polygon
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    for name in ('figures','tables','data'):(output/name).mkdir()
    runs=clean(runs)
    required=[f'{e}-{m}' for e in ('e1','e2') for m in METHODS]+[f'k{k}-{m}' for k in KAPPAS for m in ('ALOD','AFEM','LOD')]+['control-'+x for x in ('nominal-alod','family-alod','nominal-afem','family-afem')]
    missing=[name for name in required if name not in runs]
    issues=[]
    for name in required:
        if name not in runs:continue
        r=runs[name]
        if name in ('e1-ALOD','e2-ALOD'):
            if len(r['rows'])<51:issues.append(name+': fewer than 51 accepted states')
            if {a['sample'] for a in r.get('audits',[])}!=set(range(48)):issues.append(name+': incomplete 48-member audit')
        if not r.get('complete') and not (name.startswith('k') and r.get('resource_stopped')):issues.append(name+': incomplete run')
    def savefig(fig,name):
        for suffix in ('pdf','png'):fig.savefig(output/'figures'/f'{name}.{suffix}',bbox_inches='tight',dpi=180)
        plt.close(fig)
    def displayed(exp,method,run):
        rows=run['rows']
        if run.get('display_state_limit'):return dict(run,rows=rows[:run['display_state_limit']])
        if exp=='e2' and method!='ALOD' and 'e2-ALOD' in runs:
            target=runs['e2-ALOD']['rows'][-1]['E'];end=next((i+1 for i,r in enumerate(rows) if r['E']<=target),len(rows));rows=rows[:end]
        return dict(run,rows=rows)
    def plotted(run):return [r for r in run['rows'] if r.get('cycle_complete',True) or r['state']==0]
    for exp in ('e1','e2'):
        fig,axes=plt.subplots(1,2,figsize=(9,3.6),layout='constrained')
        for method in METHODS:
            r=runs.get(f'{exp}-{method}')
            if r:
                rows=plotted(displayed(exp,method,r));axes[0].loglog([x['N_on'] for x in rows],[x['E'] for x in rows],'.-',label=method)
        r=runs.get(f'{exp}-ALOD')
        if r:
            rows=plotted(r);axes[1].loglog([x['N_on'] for x in rows],[x['E'] for x in rows],'.-',label='Exact relative error')
            eta=[x for x in rows if x.get('eta_relative') is not None]
            if eta:axes[1].loglog([x['N_on'] for x in eta],[x['eta_relative'] for x in eta],'--',label='Coarse indicator')
            for before,after in zip(rows,rows[1:]):
                if before.get('ell') and after.get('ell',before['ell'])>before['ell']:
                    axes[0].annotate(f"ell {before['ell']} to {after['ell']}",(after['N_on'],after['E']),xytext=(0,24),textcoords='offset points',arrowprops={'arrowstyle':'->'},fontsize=8)
        for ax in axes:ax.set(xlabel='Online dimension',ylabel='Relative energy error');ax.grid(True,alpha=.2);ax.legend(fontsize=8)
        savefig(fig,exp+'_error_dof')
        fig,axes=plt.subplots(1,3,figsize=(10,3.3),layout='constrained')
        for ax,(method,part) in zip(axes,[('ALOD','coarse'),('ALOD','reference'),('AFEM','coarse')]):
            geometry=runs.get(f'{exp}-{method}',{}).get('mesh',{}).get(part)
            if geometry:
                points=np.asarray(geometry['nodes']);tri=np.asarray(geometry['elements'],dtype=int)
                if tri.ndim!=2 or tri.shape[1]!=3 or tri.min()<0 or tri.max()>=len(points):raise ValueError('invalid mesh connectivity')
                ax.add_collection(LineCollection(points[tri[:,[0,1,2,0]]],linewidths=.12,color='#24536c',rasterized=True));ax.autoscale();ax.set_aspect('equal')
                if exp=='e2' and method=='ALOD':
                    angle=np.linspace(0,1.5*np.pi,160);arc=.6*np.c_[np.cos(angle),np.sin(angle)]
                    ax.add_patch(Polygon(np.vstack([[0,0],arc,[0,0]]),facecolor='#e8ba6b',alpha=.2))
                    ax.plot(arc[:,0],arc[:,1],'--',color='#956421',linewidth=1.)
            else:issues.append(f'{exp}-{method}: missing {part} mesh');ax.text(.5,.5,'Missing mesh',ha='center')
            ax.set_title(f'{method} {part}');ax.axis('off')
        savefig(fig,exp+'_final_meshes')
    fig,axes=plt.subplots(1,3,figsize=(12,3.5),layout='constrained')
    for ax,method in zip(axes,('ALOD','AFEM','LOD')):
        for k in KAPPAS:
            r=runs.get(f'k{k}-{method}')
            if r:
                rows=plotted(r);ax.loglog([x['N_on']/k**2 for x in rows],[x['E'] for x in rows],'.-',label=f'k={k}')
        ax.axhline(.01,color='grey',linestyle='--');ax.set(title=method,xlabel='Online dimension / k²',ylabel='Relative energy error');ax.legend(fontsize=8);ax.grid(True,alpha=.2)
    savefig(fig,'kappa_error_dof')
    def table(name,head,rows):
        with (output/'tables'/f'{name}.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.writer(f);w.writerow(head);w.writerows(rows)
        def tex(x):return str(x).replace('_',r'\_').replace('%',r'\%')
        lines=[r'\begin{tabular}{'+'l'*len(head)+'}', ' & '.join(map(tex,head))+r' \\ \hline']
        lines+=[' & '.join(map(tex,row))+r' \\' for row in rows];lines+=[r'\end{tabular}']
        (output/'tables'/f'{name}.tex').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    metrics={}
    for exp in ('e1','e2'):
        rows=[]
        for method in METHODS:
            run=runs.get(f'{exp}-{method}')
            if not run:continue
            run=displayed(exp,method,run);rr=run['rows'];last=rr[-1]
            window=6 if method=='LOD' else 5 if exp=='e1' and method=='ALOD' else 8 if exp=='e1' else 10
            tail=plotted(run)[-window:];p=r2=None
            if len(tail)>=3 and len({r['N_on'] for r in tail})>=3:
                x=np.log([r['N_on'] for r in tail]);y=np.log([r['E'] for r in tail]);a,b=np.polyfit(x,y,1);p=-float(a)
                r2=float(1-np.sum((y-a*x-b)**2)/np.sum((y-y.mean())**2)) if np.var(y)>0 else None
            crossings=[next((r['N_on'] for r in rr if r['E']<=t),'not reached') for t in (.02,.01)]
            rows.append([method,*crossings,last['N_on'],last['E'],p,r2])
        table(exp+'-terminal-results',['method','first N at 2%','first N at 1%','endpoint N','endpoint E','tail p','tail R2'],rows);metrics[exp]=rows
    rows=[]
    for arm in ('nominal-alod','family-alod','nominal-afem','family-afem'):
        run=runs.get('control-'+arm);r=next((r for r in run['rows'] if r['N_on']>=2000),None) if run else None
        if not r or r.get('test_count')!=24:issues.append('control-'+arm+': missing 2000-DOF crossing or 24-member held-out audit');continue
        rows.append([arm,r['state'],r['N_on'],*[100*r[k] for k in ('E','test_p50','test_p90','test_max')]])
    table('e1-family-first-crossing',['method','state','N','nominal %','test median %','test p90 %','test max %'],rows)
    budget_rows=[]
    for arm in ('nominal-alod','family-alod','nominal-afem','family-afem'):
        run=runs.get('control-'+arm)
        if not run:continue
        stats=run.get('budget_statistics')
        if not stats:
            rr=plotted(run);hi=next((i for i,r in enumerate(rr) if r['N_on']>=2000),None)
            if hi is not None and hi>0:
                lo,up=rr[hi-1],rr[hi];w=math.log(2000/lo['N_on'])/math.log(up['N_on']/lo['N_on'])
                lower={a['sample']:a for a in run.get('audits',[]) if a['state_id']==lo['state']}
                upper={a['sample']:a for a in run.get('audits',[]) if a['state_id']==up['state']}
                if lower.keys()==upper.keys() and 0 in lower:
                    errors={m:math.exp((1-w)*math.log(float(lower[m]['E']))+w*math.log(float(upper[m]['E']))) for m in lower}
                    held=[v for m,v in errors.items() if lower[m].get('role')=='test']
                    if len(held)==24:stats=dict(nominal=errors[0],p50=float(np.quantile(held,.5)),p90=float(np.quantile(held,.9)),maximum=max(held))
        if stats:budget_rows.append([arm,2000,*[100*stats[k] for k in ('nominal','p50','p90','maximum')]])
        else:issues.append('control-'+arm+': missing member-level interpolation evidence')
    table('e1-family-control',['method','interpolated N','nominal %','test median %','test p90 %','test max %'],budget_rows)
    rows=[]
    for k in KAPPAS:
        for method in ('ALOD','AFEM','LOD'):
            run=runs.get(f'k{k}-{method}')
            if not run:continue
            crossing=next((r for r in run['rows'] if r['E']<=.01),None);r=crossing or run['rows'][-1]
            status='target reached' if crossing else 'resource stopped' if run.get('resource_stopped') else 'target not reached'
            if not crossing and not run.get('resource_stopped'):issues.append(f'k{k}-{method}: target not reached without resource-stop evidence')
            rows.append([k,method,r['state'],r['N_on'],100*r['E'],status])
    table('kappa-target-summary',['k','method','state','N','error %','status'],rows)
    settings=[]
    for name in ('e1-ALOD','e2-ALOD',*[f'k{k}-ALOD' for k in KAPPAS]):
        for key,value in runs.get(name,{}).get('config',{}).items():settings.append([name,key,json.dumps(value)])
    inventory=read(Path(__file__).resolve().parents[1]/'configs/paper_v4_inventory.json')
    table('adaptive-config',['setting','E1','E2','E3'],inventory['settings'])
    table('observed-run-config',['run','setting','value'],settings)
    # Canonical input can regenerate the package without any original directories.
    (output/'data/runs.json.gz').write_bytes(gzip.compress(json.dumps(dict(schema=1,runs=runs),separators=(',',':')).encode(),compresslevel=4,mtime=0))
    report=dict(schema=1,manuscript='ALOD version4',manuscript_sha256=inventory['manuscript_sha256'],figures=list(FIGURES),tables=list(TABLES),
        complete=not missing and not issues,missing_runs=missing,issues=issues,scientific_notes=['The manuscript family-control text says first N >= 2000, but its supplied published table values are log-budget interpolation at N = 2000. Both estimands are exported separately; this package does not silently reconcile the manuscript.'],
        fit_policy='E1: ALOD five completed cycles, AFEM/UFEM eight states, LOD six states. E2: ALOD ten completed cycles, AFEM/UFEM ten displayed states, LOD six displayed states. E2 controls end at their first crossing of the final ALOD error.',
        origins={k:r.get('origin') for k,r in runs.items()},metrics=metrics)
    (output/'tools').mkdir();(output/'configs').mkdir()
    project=Path(__file__).resolve().parents[1]
    shutil.copyfile(Path(__file__),output/'tools/paper_delivery.py')
    shutil.copyfile(project/'configs/paper_v4_inventory.json',output/'configs/paper_v4_inventory.json')
    shutil.copyfile(project/'requirements-paper.txt',output/'requirements-paper.txt')
    (output/'README.md').write_text('# ALOD version 4 result delivery\n\n'+('Complete inventory.\n' if report['complete'] else 'Partial inventory: see delivery.json for missing evidence.\n')+'\nHistorical imports are not measurements of the current executable. Tables report actual crossings; resource-stopped endpoints are explicitly labelled. All five figures and five table generators are included. Rebuild from data/runs.json.gz with tools/paper_delivery.py --canonical.\n',encoding='utf-8')
    report['files']={str(p.relative_to(output)).replace('\\','/'):sha(p) for p in sorted(output.rglob('*')) if p.is_file()}
    (output/'delivery.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    if not report['complete'] and not allow_partial:raise ValueError('paper evidence incomplete; see delivery.json (use --allow-partial for an explicitly partial package)')
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spec',type=Path,help='Schema 1 JSON with runs mapping names to native directories.')
    p.add_argument('--canonical',type=Path)
    p.add_argument('--executable',type=Path,default=Path(__file__).resolve().parents[1]/'build/alod_run')
    p.add_argument('--archive-e1',type=Path);p.add_argument('--archive-e2',type=Path);p.add_argument('--archive-e3',type=Path)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--allow-partial',action='store_true')
    a=p.parse_args()
    if sum((a.spec is not None,a.canonical is not None,a.archive_e1 is not None))!=1:p.error('select one input mode')
    if a.canonical:
        data=read(a.canonical)
        if data.get('schema')!=1:raise ValueError('unsupported canonical schema')
        runs=data['runs']
    elif a.spec:
        spec=read(a.spec)
        if spec.get('schema')!=1:raise ValueError('unsupported input spec')
        runs={name:native(a.spec.parent/path,a.executable if a.executable.exists() else None) for name,path in spec['runs'].items()}
    else:
        if not a.archive_e2 or not a.archive_e3:p.error('all three archive roots are required')
        runs=import_archives(a.archive_e1,a.archive_e2,a.archive_e3)
    result=deliver(runs,a.output,allow_partial=a.allow_partial)
    print(json.dumps({k:result[k] for k in ('complete','missing_runs','issues')}))

if __name__=='__main__':main()
