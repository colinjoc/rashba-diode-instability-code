import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,sys,time,hashlib,resource,importlib.util
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from uniform import Uniform

u=Uniform(cpu_cap=560)
out=Path('science')
spec=importlib.util.spec_from_file_location('prepared','evidence/engineering/primary_response.py')
prep=importlib.util.module_from_spec(spec);spec.loader.exec_module(prep)
checks=[]
for x in [[.5,1.6,0,0],[.07,1.5,.9,1.2]]:
    r=u.evaluate(x,64,128);rr,j,f=prep.uniform(*x,64,128)
    checks.append({'x':x,'prepared_differences':[r['g'][0]/x[1]-rr,r['g'][1]-j,r['F']-f]})
    step=np.array([1e-6,1e-5,1e-5,1e-5])
    jj=[]
    for i in range(4):
        delta=np.eye(4)[i]*step[i]
        plus=u.evaluate(np.array(x)+delta,64,128)
        minus=u.evaluate(np.array(x)-delta,64,128)
        jj.append((np.array(plus['g'])-np.array(minus['g']))/(2*step[i]))
    ja=np.array(r['jac']);jd=np.array(jj).T
    checks[-1].update(analytic_jacobian=ja.tolist(),central_check=jd.tolist(),
                     max_scaled_jacobian_discrepancy=float(np.max(abs(ja-jd)/np.maximum(1,abs(ja)))))
    print('CHECK',checks[-1],flush=True)

class Cached:
    def __init__(self,nt,nw):self.nt=nt;self.nw=nw;self.x=None
    def get(self,x):
        if self.x is None or not np.array_equal(x,self.x):
            self.r=u.evaluate(x,self.nt,self.nw);self.x=np.array(x)
        return self.r
    def fun(self,x):return np.array(self.get(x)['g'])
    def jac(self,x):return np.array(self.get(x)['jac'])

roots=[]
try:
    # Numerical recovery seeds inside the frozen source-reproduction bracket.
    for seed in [[.05,1.6,.9,1.25],[.05,1.4,1.2,1.2],[.06,1.7,.6,1.3]]:
        c=Cached(64,128)
        fit=least_squares(c.fun,seed,jac=c.jac,bounds=([.020000000001,.05,0,.5],[.099999999999,3,6,3]),
                          xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=85)
        r=c.get(fit.x);r['solver_success']=bool(fit.success);r['solver_message']=fit.message
        roots.append(r);print('ROOT',json.dumps(r),flush=True)
        u.save(out/'uniform_evaluations.json',{'checks':checks,'roots':roots})
        if np.linalg.norm(r['g'],np.inf)<1e-8:break
    Path(out/'recovery_checkpoint.json').write_text(json.dumps({'checks':checks,'roots':roots,
        'count':u.count,'cpu_seconds':time.process_time()-u.start,
        'local_process_maxrss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'status':'provisional roots only; branch and error checks not yet complete'},indent=2)+'\n')
finally:
    u.save(out/'uniform_evaluations.json',{'checks':checks,'roots':roots})
