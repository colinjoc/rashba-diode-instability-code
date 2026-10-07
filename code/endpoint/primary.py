import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time,itertools
from pathlib import Path
import numpy as np
from scipy.special import polygamma
from uniform import reduction,eseq

def encode(M):return {'real':np.asarray(M).real.tolist(),'imag':np.asarray(M).imag.tolist()}
def decode(m):return np.array(m['real'])+1j*np.array(m['imag'])
def norm(x):return float(np.linalg.norm(x,2))

class Primary:
    def __init__(self):self.count=0;self.start=time.process_time();self.records=[]
    def evaluate(self,x,K,nt=256,nw=512,nested=True,cost=1):
        if self.count+cost>207:raise RuntimeError('Frozen 207 matrix evaluation cap')
        self.count+=cost
        T,D,q,h=x;theta=2*np.pi*np.arange(nt)/nt
        co=np.cos(theta)[None,None,:];si=np.sin(theta)[None,None,:]
        lam=np.array([-1.,1.])[:,None,None];w=1+lam/4
        de=(q/2-lam*h)*co;om=((2*np.arange(nw)+1)*np.pi*T)[None,:,None]
        z=om+1j*de;O=np.sqrt(z*z+D*D);p=K[0]*co+K[1]*si;den=4*O*O+p*p
        response=[np.real(4*z*z/(O*den)),np.real(4*O/den),p*np.imag(z/(O*den))]
        ca=1.5*D*D+p*p/4;cp=.5*D*D+p*p/4
        da=15*D**4/8+5*D*D*p*p/8+p**4/16;dp=3*D**4/8+3*D*D*p*p/8+p**4/16
        t3=[de*de+ca,de*de+cp,-p*de/2]
        t5=[-(de**4+6*ca*de*de+da),-(de**4+6*cp*de*de+dp),p*(de**3+ca*de)]
        grids=list(itertools.product([64,128,256],[128,256,512])) if nested else [(nt,nw)]
        table={}
        for n,m in grids:
            if n>nt or m>nw or nt%n:continue
            sl=(slice(None),slice(0,m),slice(None,None,nt//n));ang=(slice(None),slice(None),slice(None,None,nt//n))
            C3=-polygamma(2,m+.5)/(2*(2*np.pi*T)**3)
            C5=-polygamma(4,m+.5)/(24*(2*np.pi*T)**5)
            val=[];b=[];tail3=[];tail5=[]
            for a in range(3):
                factor=2*np.pi*T/n if a<2 else -4*np.pi*T/n
                prim=(-w*response[a] if a<2 else w*response[a])[sl]*factor
                ct=(np.broadcast_to(1/om,(2,nw,nt))[sl]*factor if a<2 else np.zeros_like(prim))
                c0=2*np.log(T) if a<2 else 0.
                ta=factor*(w*t3[a])[ang]*C3;tb=factor*(w*t5[a])[ang]*C5
                v,bb=reduction(np.concatenate([(prim+ct).ravel(),ta.ravel(),tb.ravel(),[c0]]))
                _,bp=reduction(np.concatenate([prim.ravel(),ct.ravel(),ta.ravel(),tb.ravel(),[c0]]))
                val.append(v);b.append(max(bb,bp));tail3.append(float(ta.sum()));tail5.append(float(tb.sum()))
            # Analytical conjugate-row AP and PA formulae, no post-hoc Hermitian cleanup.
            M=np.array([[val[0],1j*val[2]],[-1j*val[2],val[1]]])
            tthree=np.array([[tail3[0],1j*tail3[2]],[-1j*tail3[2],tail3[1]]])
            tfive=np.array([[tail5[0],1j*tail5[2]],[-1j*tail5[2],tail5[1]]])
            bm=float(np.sqrt(b[0]**2+b[1]**2+2*b[2]**2+4*1e-24))
            table[f'{n},{m}']={'matrix':encode(M),'b_matrix':bm,'entry_bounds':b,
                 'tail3':encode(tthree),'tail5':encode(tfive),'hermiticity_defect':norm(M-M.conj().T)}
        result={'x':list(map(float,x)),'K':list(map(float,K)),'table':table,'evaluation_count':self.count}
        self.records.append(result);return result

def certificate(r):
    tab=r['table']
    def mat(n,m):return decode(tab[f'{n},{m}']['matrix'])
    def bb(n,m):return tab[f'{n},{m}']['b_matrix']
    axes=[]
    for axis in [0,1]:
        seqs=[]
        for fixed in ([128,256,512] if axis==0 else [64,128,256]):
            pairs=[(n,fixed) for n in [64,128,256]] if axis==0 else [(fixed,m) for m in [128,256,512]]
            v=[mat(*p) for p in pairs];b=[bb(*p) for p in pairs]
            d1=norm(v[1]-v[0]);d2=norm(v[2]-v[1]);noise=b[0]+2*b[1]+b[2]
            seqs.append({'fixed':fixed,'pass':d2<=.75*d1+noise,'d1':d1,'d2':d2,'noise':noise,'error':4*(d2+b[1]+b[2])+b[2]})
        axes.append({'pass':all(s['pass'] for s in seqs),'error':max(s['error'] for s in seqs),'sequences':seqs})
    def mixed(n,m):
        pairs=[(256,512),(n,512),(256,m),(n,m)]
        return sum(s*mat(*p) for s,p in zip([1,-1,-1,1],pairs)),sum(bb(*p) for p in pairs)
    jl,bl=mixed(64,128);jm,bm=mixed(128,256)
    mixok=norm(jm)<=.75*norm(jl)+bl+bm;emix=4*(norm(jm)+bm)
    fine=tab['256,512'];tf=decode(fine['tail5']);tt=decode(fine['tail3']);bf=fine['b_matrix']
    etail=2*norm(tf)+2*bf
    T,D,q,h=r['x'];k=np.linalg.norm(r['K']);first=(2*512+1)*np.pi*T
    tailok=(first>=4*max(D,abs(q/2-h),abs(q/2+h),k/2) and (norm(tf)<=.25*norm(tt) or max(norm(tf),norm(tt))<=bf))
    return {'pass':bool(all(a['pass'] for a in axes) and mixok and tailok),
      'axes':axes,'mixed':{'pass':mixok,'low':norm(jl),'medium':norm(jm),'bound_low':bl,'bound_medium':bm,'error':emix},
      'tail':{'pass':bool(tailok),'C3_norm':norm(tt),'C5_norm':norm(tf),'error':etail,'first_omitted_frequency':first},
      'EP':axes[0]['error']+axes[1]['error']+emix+etail+bf}

def scalar_table(r,v):
    return {key:float(np.vdot(v,decode(val['matrix'])@v).real) for key,val in r['table'].items()}
