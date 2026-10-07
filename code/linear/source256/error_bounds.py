"""Analytic real-domain bounds with original empirical binary64 arithmetic floor.
No physical defaults or admission. Bounds are conservative, not fitted.
"""
import math


def derivative_bounds(w,d,D,k,Dmin=None):
    """Entrywise derivatives with sharper real-axis |Omega²|>=2wD (w<D), |Q|>=8wD."""
    if w<=0 or min(d,D,k)<0:raise ValueError('invalid positive-domain bounds')
    Z=w+d; W=Z+D; C=D*D+k*Z
    lowD=D if Dmin is None else Dmin
    mW=2*w*lowD if w<lowD else w*w+lowD*lowD
    O=math.sqrt(mW);Q=4*mW
    dd=max(8*(2*Z/(O*Q)+Z**3/(O**3*Q)+8*Z**3/(O*Q**2)),
           8*(Z/(O*Q)+8*W*Z/Q**2),4*(k/(O*Q)+C*Z/(O**3*Q)+8*C*Z/(O*Q**2)))
    dg=max(8*Z*Z*(D/(O**3*Q)+8*D/(O*Q**2)),
           8*(D/(O*Q)+8*W*D/Q**2),4*(2*D/(O*Q)+C*D/(O**3*Q)+8*C*D/(O*Q**2)))
    dk=max(16*k*Z*Z/(O*Q**2),16*W*k/Q**2,4*(Z/(O*Q)+2*C*k/(O*Q**2)))
    hh=max(8*Z*Z/(O*Q),8*W/Q,4*C/(O*Q))
    # w and delta derivatives have identical modulus majorants.
    return dict(delta=dd,gap=dg,k=dk,omega=dd,H=hh)


def omitted_tail_bound(T,N,d,D,k,vnorm):
    """All omitted terms, including leading term AND remainder, not C3 alone.

    Set x=1/w; r=1/[4(d+D+k)]. On |x|=r scalar linear Riccati
    kernels have |sqrt W|>=sqrt(.5), |4W+k²x²|>=1.9375;
    every entry of J=2xI-H has modulus<16r. Inversion pairing
    makes J odd and cancels x term. Cauchy gives |J(x)| <=
    16*x³/r²/[1-(x/r)²]. This finite-dimensional LINEAR kernel
    statement is never transferred to the nonlinear differential model.
    """
    S=d+D+k
    if T<=0 or N<1 or S<=0:raise ValueError('tail domain')
    r=1/(4*S);w=(2*N+1)*math.pi*T
    if w*r<=1:raise ValueError('tail analytic exclusion fails: increase frozen cutoff prospectively')
    sum3=1/w**3+1/(4*math.pi*T*w*w)
    bound=math.pi*T*2*vnorm*16/r**2*sum3/(1-1/(w*r)**2)
    return {'bound':bound,'inverse_frequency_radius':r,'first_omitted_frequency':w,
            'scope':'complete linear omitted tail analytic majorant; empirical arithmetic remains separate'}


def linear_complete_bounds(source,radii,Ntheta,Nw,witness,floor=1e-12):
    T,D,q,h=source;rt,rD,rq,rh=radii
    Tmin=T-rt;Dmax=D+rD;dmax=(abs(q)+rq)/2+abs(h)+rh;k=math.hypot(*witness['K'])
    v=[complex(a,b) for a,b in zip(witness['v']['real'],witness['v']['imag'])]
    vn=sum(abs(x) for x in v)**2
    if Tmin<=0:raise ValueError('source box intersects nonpositiveT')
    angular=source_error=absolute=0.
    source_error=2*rt/Tmin*sum(abs(x)**2 for x in v)
    for n in range(Nw):
        w=(2*n+1)*math.pi*Tmin;b=derivative_bounds(w,dmax,Dmax,k,D-rD)
        angular+=math.pi*(T+rt)*2*vn*math.pi/Ntheta*(b['delta']*dmax+b['k']*k)
        source_error+=math.pi*(T+rt)*2*vn*(b['delta']*(rq/2+rh)+b['gap']*rD
                      +(2/w**2+b['omega'])*w/Tmin*rt)
        source_error+=math.pi*rt*2*vn*(2/w+b['H'])
        absolute+=math.pi*(T+rt)*2*vn*(2/w+b['H'])
    tail=omitted_tail_bound(Tmin,Nw,dmax,Dmax,k,vn)
    arithmetic=floor*(1+absolute+abs(2*math.log(Tmin))*vn+tail['bound'])
    return dict(angular_continuum_bound=angular,source_box_correlated_bound=source_error,
                omitted_frequency_bound=tail['bound'],empirical_arithmetic_floor=arithmetic,
                total_excluding_original_Eep=angular+source_error+tail['bound']+arithmetic,
                arithmetic_rigorous=False,tail_details=tail,
                scope='Complete analytic linear domain majorants under ORIGINAL empirical arithmetic policy; no nonlinear certificate')


def nonlinear_residual_radius(L,P,rho,beta,a0,E,D):
    """Validated inputs required; returns sufficient scalar a-posteriori ball."""
    B=abs(beta);A=abs(a0);c=2*B*A*E;g=B*(abs(D)+B*E)
    discr=(1-L*(c+2*g*P))**2-4*L*L*g*rho
    if 1-L*(c+2*g*P)<=0 or discr<0:return {'pass':False,'reason':'no sufficient contraction enclosure'}
    e=2*L*rho/(1-L*(c+2*g*P)+math.sqrt(discr))
    lip=L*(c+2*g*(P+e))
    return {'pass':lip<1,'radius':e,'lipschitz':lip,'scope':'sufficient scalar enclosure only; input/rho validation and causal tube still mandatory'}


def angular_strip_bound(w,d,Dmin,Dmax,k,N,vnorm):
    """Continuous angular trapezoid bound, analytically separated branch/poles.

    On the real axis |z²+D²|>=2wD if w<D, else >=w²+D²;
    |4(z²+D²)+p²|>=8wD if w<D, else >=4(w²+D²).
    Perturb cos(theta) through |Imtheta|<=s by at most exp(s)-1;
    choose s with both polynomial perturbations <=half lower margin.
    Zero-free continuation inherits the physical causal root; no sampled poles.
    """
    mW=2*w*Dmin if w<Dmin else w*w+Dmin*Dmin
    mQ=4*mW
    def root(a,b,c):return 2*c/(b+math.sqrt(b*b+4*a*c)) if a else c/b
    yw=root(d*d,2*d*(w+d),mW/2) if d else 1.
    a=4*d*d+k*k; yq=root(a,2*a+8*w*d,mQ/2) if a else 1.
    s=math.log1p(min(yw,yq));Z=w+d*math.cosh(s);O=math.sqrt(mW/2);Q=mQ/2
    HH=max(8*Z*Z/(O*Q),8*(Z+Dmax)/Q,4*(Dmax*Dmax+k*math.cosh(s)*Z)/(O*Q))
    exponent=N*s
    bound=0. if exponent>700 else 2*(2/w+HH)*vnorm/math.expm1(exponent)
    return {'strip':s,'bound':bound,'Omega_modulus_min':O,'transport_product_min':Q}


def corrected_tail(source,N,witness,radii):
    T,D,q,h=source;rt,rD,rq,rh=radii;Tmin=T-rt;Dmax=D+rD
    k=math.hypot(*witness['K']);d=(abs(q)+rq)/2+abs(h)+rh
    v=[complex(a,b) for a,b in zip(witness['v']['real'],witness['v']['imag'])];vn=sum(abs(x) for x in v)**2
    S=d+Dmax+k;r=1/(4*S);w=(2*N+1)*math.pi*Tmin
    if w*r<=1:raise ValueError('linear tail exclusion insufficient')
    sum5=1/w**5+1/(8*math.pi*Tmin*w**4)
    remainder=math.pi*(T+rt)*2*vn*16/r**4*sum5/(1-1/(w*r)**2)
    coeff=0.
    for lam in (-1,1):
        dl=q/2-lam*h;weight=1+lam*.25
        aa=dl*dl+k*k/4+3*D*D;pp=dl*dl+k*k/4+D*D;ap=1j*dl*witness['K'][0]
        coeff+=weight*(v[0].conjugate()*(aa*v[0]+ap*v[1])+v[1].conjugate()*(-ap*v[0]+pp*v[1])).real
    stop=N+65536
    s3=math.fsum(1/((2*n+1)*math.pi*T)**3 for n in range(N,stop))
    omitted_s3=1/(4*math.pi*T*((2*stop-1)*math.pi*T)**2)
    correction=math.pi*T*coeff*s3
    leading_sum_error=math.pi*T*abs(coeff)*omitted_s3
    # Uniform source perturbation of analytic C3 and the T^-2 prefactor.
    abscoef=2*vn*(d*d+k*k/4+3*Dmax*Dmax+2*d*k)
    dcoef=2*vn*(2*d*(rq/2+rh)+6*Dmax*rD+2*k*(rq/2+rh))
    source_error=math.pi*(T+rt)*s3*(dcoef+2*abscoef*rt/Tmin)
    return {'correction':correction,'remainder':remainder+leading_sum_error,
            'source_error':source_error,'finite_leading_sum_error':leading_sum_error,
            'linear_C3_witness':coeff,'leading_terms':65536,
            'scope':'source-derived leading term plus complete analytic higher-order majorant, not nonlinear tail'}
