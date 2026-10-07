"""Admitted v2 uniform functional and analytic derivatives; no changed physics.

F derivatives in D,q are radial-integrated source derivatives through order four.
T and h derivatives are obtained by the real chain rule before taking Re.
No complex-step differentiation of a real-valued function is used.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json, time, ctypes, math
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.special import polygamma

HERE=Path(__file__).resolve().parent
lib=ctypes.CDLL(str(HERE/'sums.so'))
lib.compensated.argtypes=[ctypes.POINTER(ctypes.c_double),ctypes.c_size_t,
                         ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double)]

def reduction(terms):
    a=np.ascontiguousarray(terms,dtype=np.float64).ravel()
    s=ctypes.c_double(); ab=ctypes.c_double()
    lib.compensated(a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),a.size,
                    ctypes.byref(s),ctypes.byref(ab))
    pair=float(np.sum(a))
    return s.value,max(1e-12*(1+ab.value),4*abs(pair-s.value))

# The branch of all half-integer powers is the positive-Matsubara square root.
d,z,delta=sp.symbols('d z delta')
inds=[(m,n) for k in range(5) for m in range(k+1) for n in [k-m]]
base=sp.sqrt(z*z+d*d)-z
expr=[sp.diff(base,d,m,z,n) for m,n in inds]
# Rewrite the zero-order expression to avoid Omega-z cancellation.
expr[0]=d*d/(sp.sqrt(z*z+d*d)+z)
O=sp.lambdify((d,z),expr,'numpy',cse=True)
p3=d*d*delta*delta/2+d**4/8
p5=-(d*d*delta**4/2+3*d**4*delta*delta/4+d**6/16)
tails3=sp.lambdify((d,delta),[sp.diff(p3,d,m,delta,n) for m,n in inds],'numpy',cse=True)
tails5=sp.lambdify((d,delta),[sp.diff(p5,d,m,delta,n) for m,n in inds],'numpy',cse=True)
idx={mn:i for i,mn in enumerate(inds)}

class Uniform:
    def __init__(self,cap=900,cpu_cap=600):
        self.count=0;self.start=time.process_time();self.cpu_cap=cpu_cap;self.cap=cap
        self.log=[]

    def evaluate(self,x,nt,nw,tail5=True,record=True):
        if self.count>=self.cap or time.process_time()-self.start>self.cpu_cap:
            raise RuntimeError('Frozen uniform/primary allocation limit reached')
        self.count+=1
        T,D,q,h=map(float,x)
        if T<=0 or D<=0: raise ValueError('nonpositive T or Delta')
        theta=2*np.pi*np.arange(nt)/nt
        co=np.cos(theta)[None,None,:]
        lam=np.array([-1.,1.])[:,None,None];w=1+lam/4
        de=(q/2-lam*h)*co
        om=((2*np.arange(nw)+1)*np.pi*T)[None,:,None]
        zz=om+1j*de
        vals=O(D,zz); t3=tails3(D,de);t5=tails5(D,de)
        C3=-polygamma(2,nw+.5)/(2*(2*np.pi*T)**3)
        C5=-polygamma(4,nw+.5)/(24*(2*np.pi*T)**5) if tail5 else 0.
        H={}; bH={}; Hlam={}; HT={}
        for i,(m,n) in enumerate(inds):
            fac=(.5j*co)**n
            response=-w*np.real(vals[i]*fac)
            # Keep counterterm and determinant primitive magnitudes separately.
            cd=([D*D,2*D,2.,0.,0.][m]/2 if n==0 else 0.)
            ct=cd/om*np.ones((2,1,nt))
            factor=2*np.pi*T/nt
            terms=(response+ct)*factor
            per=terms.sum(axis=(1,2))
            tail3=(2*np.pi*T*w*C3*t3[i]*(co/2)**n/nt)
            tail5v=(2*np.pi*T*w*C5*t5[i]*(co/2)**n/nt)
            tail3=np.broadcast_to(tail3,(2,1,nt));tail5v=np.broadcast_to(tail5v,(2,1,nt))
            per+=np.sum(tail3+tail5v,axis=(1,2))
            logterm=2*cd*np.log(T)
            total,bb=reduction(np.concatenate([terms.ravel(),tail3.ravel(),tail5v.ravel(),[logterm]]))
            # S is before counterterm/determinant cancellation, as required.
            _,bpre=reduction(np.concatenate([(factor*response).ravel(),(factor*ct).ravel(),tail3.ravel(),tail5v.ravel(),[logterm]]))
            H[m,n]=total;bH[m,n]=max(bb,bpre);Hlam[m,n]=per
            if m+n<=3:
                nxt=vals[idx[m,n+1]]
                # The 1/omega counterterm cancels in d/dT[T/omega].
                tt=-2*np.pi*w*np.real((vals[i]+om*nxt)*fac)/nt
                tv=float(np.sum(tt))+2*cd/T+float(np.sum(-2*tail3/T-4*tail5v/T))
                HT[m,n]=tv
        F_D=H[1,0];A=H[2,0];B=H[1,1];r=-B/A
        Jp=2*(H[0,2]-B*B/A)
        Jpp=2*(H[0,3]+3*r*H[1,2]+3*r*r*H[2,1]+r**3*H[3,0])
        # R is the canonical (undivided) gap residual = F_D/2.
        g=np.array([F_D/2,2*H[0,1],Jp,Jpp])
        def grad(m,n):
            return np.array([HT[m,n],H[m+1,n],H[m,n+1],
                             np.dot(np.array([2.,-2.]),Hlam[m,n+1])])
        ga=grad(2,0);gb=grad(1,1);gr=-gb/A+B*ga/(A*A)
        jac=np.stack([grad(1,0)/2,2*grad(0,1),
            2*(grad(0,2)-2*B*gb/A+B*B*ga/A**2),
            2*(grad(0,3)+3*r*grad(1,2)+3*r*r*grad(2,1)+r**3*grad(3,0)
               +3*gr*(H[1,2]+2*r*H[2,1]+r*r*H[3,0]))])
        quartic=(H[0,4]+4*r*H[1,3]+6*r*r*H[2,2]+4*r**3*H[3,1]+r**4*H[4,0]
                 -3*(H[1,2]+2*r*H[2,1]+r*r*H[3,0])**2/A)
        result={'x':list(map(float,x)),'nt':nt,'nw':nw,'tail5':tail5,
                'g':g.tolist(),'jac':jac.tolist(),'F':H[0,0], 'amplitude_curvature':A,
                'amplitude_slope':r,'relaxed_quartic':quartic,
                'F_h':float(np.dot(np.array([2.,-2.]),Hlam[0,1])),
                'F_derivatives':{f'{m},{n}':v for (m,n),v in H.items()},
                'primitive_bounds':{f'{m},{n}':v for (m,n),v in bH.items()},
                'count':self.count,'process_cpu_seconds':time.process_time()-self.start}
        if record:self.log.append(result)
        return result

    def save(self,path,extra=None):
        Path(path).write_text(json.dumps({'evaluations':self.log,'count':self.count,
              'cpu_seconds':time.process_time()-self.start,'extra':extra},indent=2,allow_nan=False)+'\n')

def eseq(y,b):
    y=np.asarray(y);b=np.asarray(b);ds=np.diff(y,axis=0)
    d1=np.abs(ds[0]);d2=np.abs(ds[1]);noise=b[0]+2*b[1]+b[2]
    ok=(d2<=.75*d1+noise)&((ds[0]*ds[1]>=0)|(np.minimum(d1,d2)<=noise))
    return 4*(d2+b[1]+b[2])+b[2],ok
