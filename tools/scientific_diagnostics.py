"""Explicit offline scientific controls on immutable accepted checkpoints."""
import argparse
import json
from pathlib import Path
from checkpoint_io import inspect, atomic_text, digest
from run_audit import audit

MENU = {
    'localization': 'Compute Theta with the current eigenvalue-relative Ritz residual gate.',
    'full': 'Complete e/f/g audit on the accepted space.',
    'exact': 'Exact error without an auxiliary reference FEM solve.',
    'rank-zero': 'Remove enrichment on the same mesh and oversampling level.',
    'fresh': 'Retrain each RHS separately, compared with the frozen shared dictionary.',
    'pure': 'Evaluate pure singular and pure smooth E2 loads.',
    'ell': 'Rebuild oversampling on the same mesh pair; retrain the E2 shared family.',
    'rank': 'Retrain the same frozen family with each requested rank cap.',
    'radius': 'Retrain the same frozen family with each requested regional radius.',
    'reference': 'Deepen the reference mesh and inject the frozen trial/dictionary.',
    'region': 'Energy error on the origin-centred, centroid-selected physical region.',
    'two-level': 'Energy Frobenius norm of the ell+1 versus ell corrected-basis difference.',
    'quadrature': 'Increase quadrature orders on the unchanged accepted mesh.',
    'residual': 'Iterative refinement of E1 reference and coarse solves.',
    'mesh': 'Export the accepted coarse and reference triangulations.'
}

def run(checkpoint, output, executable, names, *, ells=(2,3), ranks=(1,2,4), radii=(.3,.6,.9),
        levels=(1,), region=.6, boost=2, threads=None, member_ids=None):
    checkpoint, meta=inspect(checkpoint, executable)
    if meta['phase']!=0:raise ValueError('diagnostics require an accepted checkpoint')
    cfg=meta['config'];problem=cfg['problem'];lod=cfg['method']=='ALOD'
    if names==['all']:names=list(MENU)
    if not names or len(set(names))!=len(names) or set(names)-MENU.keys():raise ValueError('invalid diagnostic selection')
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    report=dict(schema=1,checkpoint_sha256=digest(checkpoint),source_state=meta['state_id'],
                definitions=MENU,completed=[],skipped=[],failed=[],complete=False)
    def save():atomic_text(output/'diagnostics.json',json.dumps(report,indent=2)+'\n')
    for name in names:
        reason=None
        if name in ('fresh','pure','rank','radius') and (problem!='E2' or not lod):reason='Requires E2 ALOD.'
        if name in ('ell','rank-zero','two-level','reference','localization') and not lod:reason='Requires an ALOD checkpoint.'
        if name=='two-level' and meta['ell']>=4:reason='The supported oversampling maximum is four.'
        if name=='residual' and (problem!='E1' or not lod):reason='Requires E1.'
        if reason:report['skipped'].append(dict(name=name,reason=reason));continue
        jobs={
            'localization':[dict(localization=True)],'full':[{}], 'exact':[dict(audit_mode='exact')],
            'rank-zero':[dict(drop_enrichment=True)], 'fresh':[dict(fresh=True)],
            'pure':[dict(pure=True)], 'mesh':[dict(export_mesh=True)],
            'ell':[dict(ell_override=e,retrain_shared=problem=='E2') for e in ells],
            'rank':[dict(retrain_shared=True,rank_cap=r) for r in ranks],
            'radius':[dict(retrain_shared=True,radius=r) for r in radii],
            'reference':[dict(reference_levels=l) for l in levels],
            'region':[dict(region_radius=region)],'two-level':[dict(two_level=True)],
            'quadrature':[dict(quadrature_boost=boost)],'residual':[dict(refinement_steps=2)]
        }[name]
        for index,options in enumerate(jobs):
            label=f'{name}-{index}'
            try:
                result=audit(checkpoint,output/label,executable,threads=threads,member_ids=member_ids,**options)
                report['completed'].append(dict(name=name,output=label,options=options,
                    output_sha256=result['output_sha256'],audit_complete=result['audit_complete']))
            except Exception as error:report['failed'].append(dict(name=name,output=label,error=type(error).__name__))
            save()
    if digest(checkpoint)!=report['checkpoint_sha256']:raise ValueError('source checkpoint changed')
    report['complete']=not report['failed'];save()
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--menu',action='store_true')
    p.add_argument('--checkpoint',type=Path);p.add_argument('--output',type=Path)
    p.add_argument('--executable',type=Path,default=Path(__file__).resolve().parents[1]/'build/alod_run')
    p.add_argument('--diagnostics',default='all')
    for name,default in [('ells','2,3'),('ranks','1,2,4'),('radii','.3,.6,.9'),('levels','1')]:p.add_argument('--'+name,default=default)
    p.add_argument('--region',type=float,default=.6);p.add_argument('--boost',type=int,default=2)
    p.add_argument('--threads',type=int);p.add_argument('--member-ids')
    a=p.parse_args()
    if a.menu:print(json.dumps(MENU,indent=2));return
    if a.checkpoint is None or a.output is None:p.error('--checkpoint and --output are required')
    result=run(a.checkpoint,a.output,a.executable,a.diagnostics.split(','),
        ells=list(map(int,a.ells.split(','))),ranks=list(map(int,a.ranks.split(','))),
        radii=list(map(float,a.radii.split(','))),levels=list(map(int,a.levels.split(','))),
        region=a.region,boost=a.boost,threads=a.threads,
        member_ids=None if a.member_ids is None else list(map(int,a.member_ids.split(','))))
    print(json.dumps(result));raise SystemExit(0 if result['complete'] else 1)

if __name__=='__main__':main()
