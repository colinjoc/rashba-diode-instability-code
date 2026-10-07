"""Finite-beta connected remainders and complete paired Matsubara tail.
Inert binary64 algebra; specified primitive allowances remain empirical.
"""
import math
import numpy as np
from nonlinear_response import (InclusionFailure,require,norm,add,mul,zero,scalar,
    pmul,padd,pscale,peval,shift_bound,interval_mul,interval_add,interval_scale,isum)

def intervals_product(a,b,guard):
    c=[(zero(),0.) for k in range(len(a)+len(b)-1)]
    for j,x in enumerate(a):
        guard()
        for k,y in enumerate(b):c[j+k]=interval_add(c[j+k],interval_mul(x,y))
    return c

def residual_intervals(P,eP,eta,a,O,D,uK,ra,rO,rD,ruK,guard):
    """Exact source Taylor recurrence kills orders <=degree symbolically.
    Return only genuinely omitted orders with full support and interval radii.
    Low-order binary64 defect is separately audited in runner.
    """
    pp=list(zip(P,eP));star=(np.conj(eta[::-1]),0.);sq=intervals_product(pp,pp,guard)
    terms=[(zero(),0.) for k in range(max(len(P)+1,len(sq)+2))]
    for k,x in enumerate(pp):terms[k+1]=interval_add(terms[k+1],interval_scale(interval_mul(star,x),2*a,2*ra))
    for k,x in enumerate(sq):
        terms[k+1]=interval_add(terms[k+1],interval_scale(x,D,rD))
        terms[k+2]=interval_add(terms[k+2],interval_mul(star,x))
    return [x if k>=len(P) else (zero(),0.) for k,x in enumerate(terms)]

def force_identity_intervals(U,V,F,eU,eV,eF,a,ra,guard):
    uu=list(zip(U,eU));vv=list(zip(V,eV));ff=list(zip(F,eF));Z=(zero(),0.)
    uv=intervals_product(uu,vv,guard);q=1+a*a;rq=2*abs(a)*ra+ra*ra
    qI=(scalar(q),rq);aI=(scalar(a),ra)
    count=max(len(U),len(V));h=[Z for k in range(max(count+1,len(uv)+2))]
    N=[Z for k in range(max(count,len(uv)+1))]
    for k in range(count):
        u=uu[k] if k<len(uu) else Z;v=vv[k] if k<len(vv) else Z
        h[k+1]=interval_mul(aI,interval_add(u,v));N[k]=interval_add(u,interval_scale(v,-a*a,rq))
    for k,x in enumerate(uv):h[k+2]=interval_add(h[k+2],x);N[k+1]=interval_add(N[k+1],interval_scale(interval_mul(aI,x),-1))
    den=[interval_mul(qI,qI)]+[interval_mul(qI,x) for x in h[1:]]
    out=intervals_product(den,ff,guard)
    for k,x in enumerate(N):out[k]=interval_add(out[k],interval_scale(x,-2))
    return [x if k>=len(F) else Z for k,x in enumerate(out)]

def residual_cell(terms,b,r):
    nominal=shift_bound([x[0] for x in terms],b,r)
    return nominal+math.fsum(x[1]*(abs(b)+r)**k for k,x in enumerate(terms))

def fourier_precondition(P,eta,a,O,D,uK,b,r,*,ra,rO,rD,ruK,eP,H,floor,guard):
    """Infinite Fourier block inverse. Finite inverse defects are included.
    B=diag(A_F^-1,diag(2Omega+i*m*uK)^-1 on H) with parameter-dependent
    exact high diagonal. Positive 2x2 norm bounds retain ALL cross modes.
    No projected matrix eigenvalue is substituted for full operator validity.
    """
    guard();pc=peval(P,b);ps=shift_bound(P,b,r)
    ep=math.fsum(e*(abs(b)+r)**k for k,e in enumerate(eP));drift=0.
    # Difference of full profile polynomial over beta cell, no sampled maximum.
    translated=[]
    for k in range(1,len(P)):
        translated.append(math.fsum(math.comb(j,k)*abs(b)**(j-k)*norm(P[j]) for j in range(k,len(P))))
    drift=math.fsum(x*r**(k+1) for k,x in enumerate(translated))
    star=np.conj(eta[::-1]);E=norm(eta);B=abs(b)+r;A=abs(a)+ra;Dp=abs(D)+rD
    c=add(2*b*a*star,2*b*mul(add(scalar(D),b*star),pc))
    dc=2*(r*abs(a)+B*ra)*E
    dc+=2*(r*abs(D)+B*rD+(2*abs(b)+r)*r*E)*norm(pc)
    dc+=2*B*(Dp+B*E)*(drift+ep)+floor*(1+norm(c))
    modes=np.arange(-H,H+1);d=2*O+1j*modes*uK
    degree=len(c)//2;AF=np.diag(d)
    for j,m in enumerate(modes):
        for k,n in enumerate(modes):
            index=m-n+degree
            if 0<=index<len(c):AF[j,k]+=c[index]
    inv=np.linalg.inv(AF);bn=float(np.linalg.norm(inv,1));defect=float(np.linalg.norm(np.eye(len(modes))-inv@AF,1))+floor*(1+bn*float(np.linalg.norm(AF,1)))
    require(np.all(np.isfinite(inv)) and math.isfinite(bn),'finite Fourier preconditioner')
    external=list(range(-H-degree,-H))+list(range(H+1,H+degree+1))
    FH=np.zeros((len(modes),len(external)),dtype=np.complex128)
    HF=np.zeros((len(external),len(modes)),dtype=np.complex128)
    for j,m in enumerate(modes):
        for k,n in enumerate(external):
            ij=m-n+degree;ik=n-m+degree
            if 0<=ij<len(c):FH[j,k]=c[ij]
            if 0<=ik<len(c):HF[k,j]=c[ik]
    causal=O.real-rO;require(causal>0,'preconditioned source causal margin',margin=causal)
    min_u=max(0.,abs(uK)-ruK)
    high_im=max(0.,(H+1)*min_u-2*(abs(O.imag)+rO))
    highL=1/math.hypot(2*causal,high_im)
    qff=defect+bn*(2*rO+H*ruK+dc)
    qfh=(float(np.linalg.norm(inv@FH,1))+floor*(1+bn*float(np.linalg.norm(FH,1))) if external else 0.)+bn*dc
    qhf=highL*((float(np.linalg.norm(HF,1)) if external else 0.)+dc)
    qhh=highL*(norm(c)+dc)
    radius=(qff+qhh+math.sqrt((qff-qhh)**2+4*qfh*qhf))/2
    require(radius<1,'infinite Fourier preconditioned inverse',spectral_norm_majorant=radius,qff=qff,qfh=qfh,qhf=qhf,qhh=qhh,H=H)
    # Strict positive weights avoid a zero tail weight at beta=0.
    w=math.sqrt((qhf+floor)/(qfh+floor));q=max(qff+qfh*w,qhh+qhf/w)
    require(q<1,'weighted Fourier block contraction',q=q,high_weight=w)
    gain=max(bn,highL/w);normfactor=1+w
    return dict(q=q,gain=gain,normfactor=normfactor,finite_inverse_l1=bn,high_inverse_l1=highL,
        high_weight=w,centre_profile_l1=ps+ep,operator_variation_l1=dc,finite_inverse_residual=defect)

def branch_remainder(P,eP,eta,a,O,D,uK,residual,b,r,*,ra,rO,rD,ruK,H,cap,floor,guard):
    pc=fourier_precondition(P,eta,a,O,D,uK,b,r,ra=ra,rO=rO,rD=rD,ruK=ruK,eP=eP,H=H,floor=floor,guard=guard)
    rho=residual_cell(residual,b,r)+floor*(1+pc['centre_profile_l1'])
    B=abs(b)+r;g=B*(abs(D)+rD+B*norm(eta));q=pc['q'];g0=pc['gain'];nf=pc['normfactor']
    # e is weighted-block radius; full Fourier radius <=nf*e.
    A=1-q;aa=g0*g*nf*nf;disc=A*A-4*aa*g0*rho
    require(disc>=0,'Taylor branch remainder radii discriminant',discriminant=disc,residual=rho,linear_margin=A,quadratic=aa)
    e=(2*g0*rho/(A+math.sqrt(disc))) if rho else 0.;full=nf*e
    # Common full l1 cap across adjacent cells makes the beta0 connection exact:
    # same P(beta) at shared boundaries, common uniqueness neighbourhood.
    capweighted=cap/nf
    require(full<=cap and g0*rho+q*capweighted+aa*capweighted**2<=capweighted and q+2*aa*capweighted<1,
        'common connected Taylor remainder cap',radius=full,cap=cap,cap_map=g0*rho+q*capweighted+aa*capweighted**2,weighted_cap=capweighted,uniqueness=q+2*aa*capweighted)
    inverse_full=pc['normfactor']*pc['gain']/(1-pc['q'])
    require(2*inverse_full*g*cap<1,'common full-l1 uniqueness across beta-cell overlaps',uniqueness=2*inverse_full*g*cap,full_inverse_upper=inverse_full,common_cap=cap)
    require(abs(a)+ra+B*(pc['centre_profile_l1']+cap)<1,'source-connected physical Riccati disk',upper=abs(a)+ra+B*(pc['centre_profile_l1']+cap))
    return dict(radius=full,residual=rho,preconditioner=pc,common_cap=cap,full_inverse_upper=inverse_full,overlap_uniqueness_upper=2*inverse_full*g*cap)

def force_remainder(U,V,F,identity,eta,a,ra,b,r,eU,eV,eF,rU,rV,floor):
    B=abs(b)+r;A=abs(a)+ra;E=norm(eta);q=abs(1+a*a);rq=2*abs(a)*ra+ra*ra;ql=q-rq
    P=shift_bound(U,b,r)+math.fsum(e*B**k for k,e in enumerate(eU))
    Q=shift_bound(V,b,r)+math.fsum(e*B**k for k,e in enumerate(eV))
    h=B*A*(P+Q)+B*B*P*Q;dh=B*A*(rU+rV)+B*B*(P*rV+Q*rU+rU*rV)
    second=ql-h;true_second=second-dh
    require(ql>0 and second>0 and true_second>0,'both uncertainty-inclusive force poles',first=ql,formal_second=second,true_second=true_second)
    N=P+A*A*Q+B*A*P*Q;dN=rU+A*A*rV+B*A*(P*rV+Q*rU+rU*rV)
    profile=E*(2*dN/(ql*true_second)+2*N*dh/(ql*true_second*second))
    identity_error=E*residual_cell(identity,b,r)/(ql*second)
    primitive=floor*(1+E*(N+dN)/(ql*true_second)+identity_error)
    return dict(error=profile+identity_error+primitive,profile_error=profile,
        omitted_force_series_error=identity_error,primitive_error=primitive,
        first_pole_lower=ql,second_pole_lower=true_second)

def tail_cf(T,rT,D,rD,r,E,N,floor):
    """Analytic COMPLEX beta disk majorant |F(beta)-F0|<=Cf/omega³.
    Direct all-order Riccati contraction, both force poles, every source/theta.
    It is used only above the frozen split, never as a low-frequency branch tube.
    """
    Tmin=T-rT;Tmax=T+rT;Dp=abs(D)+rD;require(Tmin>0,'tail temperature')
    w=(2*N+1)*math.pi*Tmin;A=Dp/(2*w);C=E*(1+A*A);R=C/w;G=Dp+r*E
    lip=(2*r*A*E+2*r*G*R)/(2*w);rhs=(E*(1+A*A)+2*r*A*E*R+r*G*R*R)/(2*w)
    require(rhs<=R and lip<1 and A+r*R<1,'complex-beta analytic tail disk',w=w,r=r,R=R,rhs=rhs,lip=lip,physical=A+r*R)
    Cu=r*(Dp*E*C+G*C*C/2)/(1-lip);H=r*Dp*C+r*r*C*C;q=1-A*A;dt=q*(q-H/w**2)
    require(dt>0,'complex-beta tail BOTH force factors',q=q,second=q-H/w**2)
    Nc=(1+A*A)*Cu+r*Dp*C*C/2;Cf=E*(2*Nc/dt+2*C*(1+A*A)**2*H/(dt*q*q))
    return dict(Cf=Cf+floor*(1+Cf),first_omitted=w,sum3=1/w**3+1/(4*math.pi*Tmin*w*w),Tmax=Tmax,radius=r)

def paired_tail(T,rT,D,rD,eps,E,N,r,floor,degree=0):
    c=tail_cf(T,rT,D,rD,r,E,N,floor);require(eps<r,'paired tail Cauchy disk',epsilon=eps,radius=r)
    first=2 if degree==0 else 2*(degree//2+1);ratio=eps/r
    # Sign sum cancels odd powers BEFORE absolute error. Path moment <=1/2;
    # sum_sigma and t integral cancel its factor2; DOS sum remains exactly2.
    bound=math.pi*c['Tmax']*2*c['Cf']*ratio**first/(1-ratio*ratio)*c['sum3']
    return dict(bound=bound+floor*(1+bound),first_omitted_even_order=first,
        first_omitted_frequency=c['first_omitted'],complex_beta_radius=r,
        all_higher_orders=True,infinite_frequency_sum3_upper=c['sum3'],force_Cf=c['Cf'])
