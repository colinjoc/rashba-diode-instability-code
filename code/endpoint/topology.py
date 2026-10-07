import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time,resource
from pathlib import Path
import numpy as np
from scipy.optimize import root
from uniform import Uniform
out=Path('science');old=json.loads((out/'uniform_evaluations.json').read_text())
u=Uniform(cpu_cap=560);u.count=old['count'];u.log=old['evaluations'];u.start=time.process_time()-old['cpu_seconds']
ep=json.loads((out/'endpoint_certificate.json').read_text());xc=np.array(ep['fine']['x'])
rows=[]
class Pair:
    def __init__(self,T):self.T=T;self.x=None
    def get(self,x):
        if self.x is None or not np.array_equal(self.x,x):
            D1,q1,D2,q2,h=x
            self.rs=[u.evaluate([self.T,D1,q1,h],256,512),u.evaluate([self.T,D2,q2,h],256,512)]
            self.x=np.array(x)
        return self.rs
    def fun(self,x):
        a,b=self.get(x);return np.array([a['g'][0],a['g'][1],b['g'][0],b['g'][1],a['F']-b['F']])
    def jac(self,x):
        a,b=self.get(x);ja=np.array(a['jac']);jb=np.array(b['jac']);j=np.zeros((5,5))
        j[:2,:2]=ja[:2,1:3];j[:2,4]=ja[:2,3];j[2:4,2:4]=jb[:2,1:3];j[2:4,4]=jb[:2,3]
        j[4]=[2*a['g'][0],a['g'][1]/2,-2*b['g'][0],-b['g'][1]/2,a['F_h']-b['F_h']]
        return j
try:
    for fraction in [.8,.9,.95,.98]:
        T=xc[0]*fraction
        sep=np.sqrt(6*1.33*(xc[0]-T)/(2*ep['fine']['relaxed_quartic']))
        q1,q2=xc[2]-sep,xc[2]+sep
        D1,D2=xc[1]+ep['fine']['amplitude_slope']*(q1-xc[2]),xc[1]+ep['fine']['amplitude_slope']*(q2-xc[2])
        p=Pair(T);fit=root(p.fun,[D1,q1,D2,q2,xc[3]],jac=p.jac,tol=1e-10)
        a,b=p.get(fit.x);res=p.fun(fit.x)
        # Locate the intervening stationary branch, not just the two minima.
        h=fit.x[-1];midseed=np.array([(a['x'][1]+b['x'][1])/2,(a['x'][2]+b['x'][2])/2])
        cache={}
        def mid(v):
            key=tuple(v)
            if cache.get('key')!=key:
                cache['r']=u.evaluate([T,v[0],v[1],h],256,512);cache['key']=key
            return cache['r']
        fm=root(lambda v:np.array(mid(v)['g'])[:2],midseed,jac=lambda v:np.array(mid(v)['jac'])[:2,1:3],tol=1e-11)
        c=mid(fm.x)
        ok=(np.max(abs(res))<1e-8 and abs(a['x'][2]-b['x'][2])>1e-3
            and a['amplitude_curvature']>0 and b['amplitude_curvature']>0 and c['amplitude_curvature']>0
            and a['g'][2]>0 and b['g'][2]>0 and c['g'][2]<0
            and min(a['x'][2],b['x'][2])<c['x'][2]<max(a['x'][2],b['x'][2])
            and c['F']>max(a['F'],b['F']) and np.max(abs(np.array(c['g'])[:2]))<1e-8)
        row={'fraction_Tcep':fraction,'pass':bool(ok),'minima':[a,b],'intervening':c,'coexistence_residual':res.tolist(),
             'q_separation':abs(a['x'][2]-b['x'][2]),'barrier_per_N0':c['F']-max(a['F'],b['F'])}
        rows.append(row);print('TOPOLOGY',fraction,ok,'qgap',row['q_separation'],'barrier',row['barrier_per_N0'],flush=True)
        u.save(out/'uniform_evaluations.json',{'phase':'coexistence topology','rows':rows})
    result={'pass':bool(all(r['pass'] for r in rows) and all(rows[i+1]['q_separation']<rows[i]['q_separation'] for i in range(3))
                    and ep['fine']['relaxed_quartic']>0),'rows':rows,'count':u.count,
            'quartic':ep['fine']['relaxed_quartic'],'cpu_seconds':time.process_time()-u.start,
            'scope':'Local source-connected coexistence and branch coalescence checks. Not proof of global root uniqueness or unrestricted pairing minimum.'}
    (out/'branch_topology.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
finally:
    u.save(out/'uniform_evaluations.json',{'phase':'coexistence topology','rows':rows})
