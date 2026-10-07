"""Finite-beta Fourier Taylor algebra. Inert; binary64 allowances are empirical.
Two same-z divided Riccati sectors, never subtraction of nearby physical roots.
Arrays retain every convolution mode; no FFT aliasing or omitted-product pruning.
"""
import math
from functools import lru_cache
import numpy as np

class InclusionFailure(ArithmeticError):
    def __init__(self,gate,**margins):
        super().__init__(gate);self.gate=gate;self.margins=margins

def require(condition,gate,**margins):
    def finite_margin(x):
        if isinstance(x,(float,np.floating)):return math.isfinite(float(x))
        if isinstance(x,(list,tuple)):return all(finite_margin(y) for y in x)
        if isinstance(x,dict):return all(finite_margin(y) for y in x.values())
        return True
    if not all(finite_margin(x) for x in margins.values()):
        raise InclusionFailure('nonfinite mathematical bound at '+gate,nonfinite_evidence=repr(margins))
    if not condition:raise InclusionFailure(gate,**margins)

def norm(a):return float(np.sum(np.abs(a)))
def scalar(a):return np.array([complex(a)],dtype=np.complex128)
def add(a,b):
    n=max(len(a),len(b));c=np.zeros(n,dtype=np.complex128)
    for x in (a,b):j=(n-len(x))//2;c[j:j+len(x)]+=x
    return c

def mul(a,b):return np.convolve(a,b)
def zero():return scalar(0)
def padd(a,b):return [add(a[k] if k<len(a) else zero(),b[k] if k<len(b) else zero()) for k in range(max(len(a),len(b)))]
def pscale(a,s):return [s*x for x in a]
def pmul(a,b,guard):
    c=[zero() for k in range(len(a)+len(b)-1)]
    for j,x in enumerate(a):
        guard()
        for k,y in enumerate(b):c[j+k]=add(c[j+k],mul(x,y))
    return c

def peval(a,b):
    v=zero()
    for x in a[::-1]:v=add(b*v,x)
    return v

@lru_cache(maxsize=256)
def _shift_weights(count,b,r):
    weights=np.zeros((count,count),dtype=np.float64)
    for j in range(count):
        for k in range(j+1):weights[j,k]=math.comb(j,k)*b**(j-k)
    powers=np.array([r**k for k in range(count)],dtype=np.float64)
    weights.setflags(write=False);powers.setflags(write=False)
    return weights,powers

def shift_bound(a,b,r):
    """Full finite binomial translation; all Fourier modes retained.
    Cached scalar weights and binary64 matrix multiplication remove Python
    array-operation overhead. Arithmetic allowances remain empirical.
    """
    if not a:return 0.
    width=max(map(len,a));coefficients=np.zeros((len(a),width),dtype=np.complex128)
    for j,x in enumerate(a):
        offset=(width-len(x))//2;coefficients[j,offset:offset+len(x)]=x
    weights,powers=_shift_weights(len(a),float(b),float(r))
    translated=weights.T@coefficients
    return math.fsum(float(x) for x in np.sum(np.abs(translated),axis=1)*powers)


def interval_add(a,b):return (add(a[0],b[0]),a[1]+b[1])
def interval_scale(a,s,rs=0):return (s*a[0],abs(s)*a[1]+rs*(norm(a[0])+a[1]))
def interval_mul(a,b):return (mul(a[0],b[0]),norm(a[0])*b[1]+norm(b[0])*a[1]+a[1]*b[1])
def isum(terms):
    v=(zero(),0.)
    for x in terms:v=interval_add(v,x)
    return v

def background(z,D,rz,rD,floor):
    O=complex(np.sqrt(z*z+D*D));O=-O if O.real<0 else O
    require(O.real>0,'nominal causal background',z=[z.real,z.imag])
    rw=2*abs(z)*rz+rz*rz+2*abs(D)*rD+rD*rD+floor*(1+abs(z*z+D*D))
    require(rw<abs(O)**2,'full-source connected background square-root disk',radius=rw,margin=abs(O)**2-rw)
    rO=rw/(abs(O)*(1+math.sqrt(1-rw/abs(O)**2)))
    require(O.real>rO,'full-source causal branch',real_margin=O.real-rO)
    s=abs(z+O);rs=rz+rO
    require(s>rs,'background a0 quotient',denominator=s,radius=rs)
    a=D/(z+O);ra=rD/(s-rs)+abs(D)*rs/(s*(s-rs))+floor*(1+abs(a))
    return O,a,rO,ra

def sector(eta,a,O,D,uK,degree,*,rO,ra,rD,ruK,floor,guard):
    """U_k recursion exact through degree, modes |m|<=k+1.
    L U0=eta-a² eta*. L Uk=-2 a eta* U{k-1}
      -D sum{i+j=k-1}Ui Uj-eta* sum{i+j=k-2}Ui Uj.
    Attached source/cell radii bound each full coefficient, not sampled corners.
    """
    star=np.conj(eta[::-1]);out=[];errs=[]
    for k in range(degree+1):
        guard()
        if k==0:numerator=interval_add((eta,0.),interval_scale((star,0.),-a*a,2*abs(a)*ra+ra*ra))
        else:
            v=interval_scale(interval_mul((star,0.),(out[k-1],errs[k-1])),-2*a,2*ra)
            quadratic=isum(interval_mul((out[i],errs[i]),(out[k-1-i],errs[k-1-i])) for i in range(k))
            v=interval_add(v,interval_scale(quadratic,-D,rD))
            if k>=2:
                quadratic=isum(interval_mul((out[i],errs[i]),(out[k-2-i],errs[k-2-i])) for i in range(k-1))
                v=interval_add(v,interval_scale(interval_mul((star,0.),quadratic),-1))
            numerator=v
        modes=np.arange(len(numerator[0]))-len(numerator[0])//2
        den=2*O+1j*modes*uK;rd=2*rO+np.abs(modes)*ruK;absden=np.abs(den)
        require(bool(np.all(absden>rd)),'full-source Taylor coefficient resolvent',order=k,min_margin=float(np.min(absden-rd)))
        x=numerator[0]/den;L=float(np.max(1/(absden-rd)))
        e=L*numerator[1]+float(np.sum(np.abs(numerator[0])*rd/(absden*(absden-rd))))+floor*(1+norm(x))
        require(math.isfinite(e) and np.all(np.isfinite(x)),'finite coefficient/radius',order=k)
        out.append(x);errs.append(e)
    return out,errs

def force_series(U,V,a,eta,degree,*,eU,eV,ra,floor,guard):
    """Formal divided force F=2(U-a²V-beta aUV)/(q(q+h)).
    Recursion retains all Fourier support. Source errors include BOTH q factors.
    """
    Z=(zero(),0.);uu=[(x,e) for x,e in zip(U,eU)];vv=[(x,e) for x,e in zip(V,eV)]
    def get(x,k):return x[k] if 0<=k<len(x) else Z
    def uv(k):return isum(interval_mul(get(uu,i),get(vv,k-i)) for i in range(k+1)) if k>=0 else Z
    q=1+a*a;rq=2*abs(a)*ra+ra*ra;require(abs(q)>rq,'first normalized-force factor',q_modulus=abs(q),radius=rq)
    iq=1/q;riq=rq/(abs(q)*(abs(q)-rq))+floor*(1+abs(iq));F=[];err=[];hh=[]
    for k in range(degree+1):
        guard()
        h=interval_add(interval_scale(interval_add(get(uu,k-1),get(vv,k-1)),a,ra),uv(k-2));hh.append(h)
        N=interval_add(interval_add(get(uu,k),interval_scale(get(vv,k),-a*a,rq)),interval_scale(uv(k-1),-a,ra))
        rhs=interval_scale(interval_scale(N,2*iq,2*riq),iq,riq)
        for j in range(1,k+1):rhs=interval_add(rhs,interval_scale(interval_mul(hh[j],(F[k-j],err[k-j])),-iq,riq))
        require(np.all(np.isfinite(rhs[0])) and math.isfinite(rhs[1]),'finite normalized-force coefficient/radius',order=k)
        F.append(rhs[0]);err.append(rhs[1]+floor*(1+norm(rhs[0])))
    # Identity residual q(q+h)F-2N, ALL beta orders and Fourier support.
    uvp=pmul(U,V,guard);N=padd(padd(U,pscale(V,-a*a)),[zero()]+pscale(uvp,-a))
    h=padd([zero()]+pscale(padd(U,V),a),[zero(),zero()]+uvp)
    denominator=padd([scalar(q*q)],pscale(h,q))
    identity=padd(pmul(denominator,F,guard),pscale(N,-2))
    projection=np.conj(eta[::-1]);coeff=[float(mul(projection,x)[len(mul(projection,x))//2].real) for x in F]
    return F,err,coeff,identity,h

def riccati_residual(P,eta,a,O,D,uK,guard):
    modes=[np.arange(len(x))-len(x)//2 for x in P]
    left=[(2*O+1j*m*uK)*x for m,x in zip(modes,P)]
    star=np.conj(eta[::-1]);sq=pmul(P,P,guard)
    right=padd([add(eta,-a*a*star)],pscale([zero()]+[mul(star,x) for x in P],-2*a))
    right=padd(right,pscale([zero()]+sq,-D));right=padd(right,pscale([zero(),zero()]+[mul(star,x) for x in sq],-1))
    return padd(left,pscale(right,-1))
