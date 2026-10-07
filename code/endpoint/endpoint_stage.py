import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,time,itertools,resource
from pathlib import Path
import numpy as np
from scipy.optimize import root
from uniform import Uniform,eseq

out=Path('science')
old=json.loads((out/'uniform_evaluations.json').read_text())
u=Uniform(cpu_cap=560);u.count=old['count'];u.log=old['evaluations'];u.start=time.process_time()-old['cpu_seconds']
seed=np.array(old['extra']['roots'][-1]['x']);roots={};solvers=[]
def newton(seed,nt,nw,tol,maxiter=20):
    x=np.array(seed,dtype=float)
    for i in range(maxiter):
        r=u.evaluate(x,nt,nw);g=np.array(r['g']);jac=np.array(r['jac'])
        if np.linalg.norm(g,np.inf)<=tol:
            r['residual_target']=tol;r['newton_iterations']=i;return r
        dx=np.linalg.solve(jac,-g)
        # Recover only inside the source bracket.
        rate=1.
        while not np.all((x+rate*dx>np.array([.02,0,0,.5]))&(x+rate*dx<np.array([.10,3,6,3]))):
            rate*=.5
            if rate<1e-6:raise RuntimeError('Newton failed source bracket')
        x+=rate*dx
    raise RuntimeError('Newton failed residual target')

try:
    for nt,nw in itertools.product([64,128,256],[128,256,512]):
        r=newton(seed,nt,nw,1e-12);roots[f'{nt},{nw}']=r
        print('ROOT',nt,nw,r['x'],r['g'],flush=True)
        u.save(out/'uniform_evaluations.json',{'phase':'nine-grid roots','roots':roots})
    for tol in [1e-8,1e-10,1e-12]:
        r=newton(seed,256,512,tol);solvers.append(r)
    fine=solvers[-1];x=np.array(fine['x']);j=np.array(fine['jac']);ji=np.linalg.inv(j)
    def b(r):return np.maximum(1e-12*np.maximum(1,abs(np.array(r['x']))),abs(np.linalg.inv(np.array(r['jac'])))@abs(np.array(r['g'])))
    def seq(rs):
        er,ok=eseq([r['x'] for r in rs],[b(r) for r in rs])
        return {'error':er.tolist(),'pass':ok.tolist()}
    seqtheta=seq([roots[f'{nt},512'] for nt in [64,128,256]])
    seqomega=seq([roots[f'256,{nw}'] for nw in [128,256,512]])
    seqsolver=seq(solvers)
    def contrast(n,m):
        rows=[roots['256,512'],roots[f'{n},512'],roots[f'256,{m}'],roots[f'{n},{m}']]
        v=sum(s*np.array(r['x']) for s,r in zip([1,-1,-1,1],rows))
        err=sum(b(r) for r in rows)
        return v,err
    jl,bl=contrast(64,128);jm,bm=contrast(128,256)
    mixok=abs(jm)<=.75*abs(jl)+bl+bm
    c3=u.evaluate(x,256,512,tail5=False)
    tail=2*abs(ji)@abs(np.array(fine['g'])-np.array(c3['g']))
    rho=(np.array(seqtheta['error'])+seqomega['error']+np.array(seqsolver['error'])+4*abs(jm)+tail)
    corners=[]
    for signs in itertools.product([-1,1],repeat=4):
        xx=x+rho*np.array(signs);r=u.evaluate(xx,256,512)
        contraction=float(np.linalg.norm(np.eye(4)-ji@np.array(r['jac']),np.inf))
        corners.append({'signs':signs,'x':xx.tolist(),'jacobian_contraction':contraction,
                        'amplitude_curvature':r['amplitude_curvature']})
    passed=(all(seqtheta['pass']) and all(seqomega['pass']) and all(seqsolver['pass']) and all(mixok)
       and np.linalg.cond(j,np.inf)<=1e10 and all(rho<=.01*np.maximum(1,abs(x)))
       and max(c['jacobian_contraction'] for c in corners)<=.5
       and min(c['amplitude_curvature'] for c in corners)>0)
    cert={'status':'endpoint empirical error check','pass':bool(passed),'fine':fine,'rho':rho.tolist(),
       'sequence_theta':seqtheta,'sequence_omega':seqomega,'sequence_solver':seqsolver,
       'mixed_low':jl.tolist(),'mixed_medium':jm.tolist(),'mixed_bounds':[bl.tolist(),bm.tolist()],
       'mixed_contraction':mixok.tolist(),'tail_radius':tail.tolist(),'condition_inf':float(np.linalg.cond(j,np.inf)),
       'corners':corners,'roots':roots,'solver_ladder':solvers,
       'local_process_maxrss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'endpoint_certificate.json').write_text(json.dumps(cert,indent=2,allow_nan=False)+'\n')
    print('CERT',json.dumps({k:v for k,v in cert.items() if k not in ['roots','solver_ladder','corners','fine']}),flush=True)
finally:
    u.save(out/'uniform_evaluations.json',{'phase':'endpoint error certificate','roots':roots,'solvers':solvers})
