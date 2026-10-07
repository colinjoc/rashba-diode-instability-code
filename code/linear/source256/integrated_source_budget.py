"""Inert continuous source envelope. Real interval algebra under empirical binary64 policy.
No response-grid/target calculations; explicit source-domain input is mandatory.
"""
import cmath,math
# Rectangles represented by builtin ((real_lo,real_hi),(imag_lo,imag_hi)).
def addi(a,b):return (a[0]+b[0],a[1]+b[1])
def muli(a,b):
 x=[a[i]*b[j] for i in (0,1) for j in (0,1)];return (min(x),max(x))
def negi(a):return (-a[1],-a[0])
def ci(x):return ((complex(x).real,)*2,(complex(x).imag,)*2)
def ca(a,b):return (addi(a[0],b[0]),addi(a[1],b[1]))
def cn(a):return (negi(a[0]),negi(a[1]))
def cm(a,b):return (addi(muli(a[0],b[0]),negi(muli(a[1],b[1]))),addi(muli(a[0],b[1]),muli(a[1],b[0])))
def sqrange(x):return (0. if x[0]<=0<=x[1] else min(x[0]**2,x[1]**2),max(x[0]**2,x[1]**2))
def polyadd(a,b):
 r=dict(a)
 for e,v in b.items():r[e]=r.get(e,0j)+v
 return {e:v for e,v in r.items() if v!=0}
def polyscale(a,b):return {e:v*b for e,v in a.items() if v*b!=0}
def polymul(a,b):
 r={}
 for e,v in a.items():
  for f,u in b.items():
   key=tuple(x+y for x,y in zip(e,f));r[key]=r.get(key,0j)+v*u
 return {e:v for e,v in r.items() if v!=0}
def deriv(a,j):return {tuple(x-(i==j) for i,x in enumerate(e)):v*e[j] for e,v in a.items() if e[j]}
def pconst(v):return {(0,0,0,0):complex(v)}
def pvar(j):return {tuple(int(i==j) for i in range(4)):1+0j}
def formulas(K,v):
 # Variables w,D,a,c; delta=a*c and projection=K*c retain their common c.
 w,D,a,c=[pvar(j) for j in range(4)];z=polyadd(w,polyscale(polymul(a,c),1j));D2=polymul(D,D);p=polyscale(c,K)
 W=polyadd(polymul(z,z),D2);Q=polyadd(polyscale(W,4),polymul(p,p));cross=v[0].conjugate()*v[1]
 N=polyscale(polyadd(polyadd(polyscale(polymul(z,z),abs(v[0])**2+abs(v[1])**2),polyscale(D2,abs(v[1])**2+1j*cross.real)),polyscale(polymul(p,z),1j*cross.imag)),8)
 # R_j= N_j WQ -N*(W_j Q/2+W Q_j); denominator W^(3/2)Q^2.
 WQ=polymul(W,Q);R=[]
 for j in (0,1,2):
  R.append(polyadd(polymul(deriv(N,j),WQ),polyscale(polymul(N,polyadd(polyscale(polymul(deriv(W,j),Q),.5),polymul(W,deriv(Q,j)))),-1)))
 joint=polyadd(polymul(N,WQ),polymul(w,R[0]))
 # Floating symbolic coefficients: no rational/Jet scalar objects or codecs.
 return W,Q,N,(joint,R[1],R[2])
def evalpoly(p,box):
 powers=[]
 for j,x in enumerate(box):
  row=[(1.,1.)]
  for n in range(max((e[j] for e in p),default=0)):row.append(muli(row[-1],x))
  powers.append(row)
 result=ci(0)
 for e,v in p.items():
  term=(1.,1.)
  for j,n in enumerate(e):term=muli(term,powers[j][n])
  result=ca(result,cm(ci(v),(term,(0.,0.))))
 return result
def modulus_bounds(w,D,a,c,K):
 # Exact real denominator polynomials; no corner-only or sampled exclusion.
 w2,D2,a2,y=[sqrange(x) for x in (w,D,a,c)]
 realW=addi(addi(w2,D2),negi(muli(a2,y)));imW=muli(muli(muli((2.,2.),w),a),c)
 realQ=addi(muli((4.,4.),addi(w2,D2)),muli(addi((K*K,K*K),negi(muli((4.,4.),a2))),y));imQ=muli((4.,4.),imW)
 mW=math.sqrt(sqrange(realW)[0]+sqrange(imW)[0]);mQ=math.sqrt(sqrange(realQ)[0]+sqrange(imQ)[0])
 if min(mW,mQ)<=0:raise ValueError('angular cell denominator exclusion unresolved')
 return (realW,imW),(realQ,imQ),mW,mQ

def real_ratios(polys,box,K):
 w,D,a,c=box;W,Q,mW,mQ=modulus_bounds(w,D,a,c,K)
 wc=sum(w)/2;dc=sum(D)/2;ac=sum(a)/2;cc=sum(c)/2
 center=(complex(wc,ac*cc)**2+dc*dc);O0=cmath.sqrt(center)
 if O0.real<=0:raise ValueError('causal center root unresolved')
 # Re sqrt((w+i delta)^2+D^2)>=w>0; therefore |O+O0|>=wmin+ReO0.
 radius=math.hypot(max(abs(x-center.real) for x in W[0]),max(abs(x-center.imag) for x in W[1]))/(w[0]+O0.real)
 O=((O0.real-radius,O0.real+radius),(O0.imag-radius,O0.imag+radius))
 denom=cm(cm(O,W),cm(Q,Q));conj=(denom[0],negi(denom[1]))
 lower=mW**3*mQ**4;upper=sqrange(denom[0])[1]+sqrange(denom[1])[1]
 if not (0<lower<=upper and math.isfinite(upper)):raise ValueError('invalid denominator interval')
 out=[]
 for poly in polys:
  real=cm(evalpoly(poly,box),conj)[0];out.append(muli(real,(1/upper,1/lower)))
 return out,mW,mQ

def continuous_source_budget(source,radii,Nw,cells,witness,tail,guard,floor=1e-12):
 T,D,q,h=map(float,source);rt,rD,rq,rh=map(float,radii);Kx,Ky=witness['K']
 if Ky!=0:raise ValueError('selected aligned witness required; no generalization silently assumed')
 if T-rt<=0 or D-rD<=0 or Nw<1 or cells<8 or cells%4:raise ValueError('invalid continuous source domain/schedule')
 v=tuple(complex(a,b) for a,b in zip(witness['v']['real'],witness['v']['imag']));norm=sum(abs(x)**2 for x in v)
 W,Q,N,polys=formulas(Kx,v);sums=[(0.,0.)]*4;minW=minQ=float('inf');abswork=0.
 for n in range(Nw):
  guard();w=((2*n+1)*math.pi*(T-rt),(2*n+1)*math.pi*(T+rt))
  for lam in (-1,1):
   a0=q/2-lam*h;ra=rq/2+rh;a=(a0-ra,a0+ra);weight=(1+lam*.25)/cells
   for j in range(cells):
    if j%32==0:guard()
    # Cos is monotone on each half-cycle; cells include all stationary points.
    c1=math.cos(2*math.pi*j/cells);c2=math.cos(2*math.pi*(j+1)/cells);c=(min(c1,c2),max(c1,c2))
    ratios,mw,mq=real_ratios(polys,(w,(D-rD,D+rD),a,c),Kx)
    minW=min(minW,mw);minQ=min(minQ,mq)
    for i,x in enumerate((ratios[0],ratios[1],muli(ratios[2],(.5,.5)),muli(ratios[2],(-lam,-lam)))):
     sums[i]=addi(sums[i],muli(x,(weight,weight)));abswork+=weight*max(map(abs,x))
 # Integrated temperature derivative includes exact logarithm and jointly differentiated piT kernel.
 # piT*2/omega_n is T independent and contributes exactly zero here.
 tder=addi((2*norm/(T+rt),2*norm/(T-rt)),negi(muli(sums[0],(math.pi,math.pi))))
 dder=negi(muli(sums[1],(math.pi*(T-rt),math.pi*(T+rt))))
 qder=negi(muli(sums[2],(math.pi*(T-rt),math.pi*(T+rt))))
 hder=negi(muli(sums[3],(math.pi*(T-rt),math.pi*(T+rt))))
 # The reused finite C3 leading sum has a nominal truncation error. Also
 # enclose its source-box variation uniformly; no zero nominal coefficient shortcut.
 vn=sum(abs(x) for x in v)**2;dm=(abs(q)+rq)/2+abs(h)+rh;km=abs(Kx);stop=Nw+tail['leading_terms']
 abscoef=2*vn*(dm*dm+km*km/4+3*(D+rD)**2+2*dm*km)
 omitted_s3=1/(4*math.pi*(T-rt)*((2*stop-1)*math.pi*(T-rt))**2)
 uniform_leading_truncation_variation=2*math.pi*(T+rt)*abscoef*omitted_s3
 components={'temperature':rt*max(map(abs,tder)),'gap':rD*max(map(abs,dder)),
             'gradient':rq*max(map(abs,qder)),'field':rh*max(map(abs,hder)),
             'corrected_tail_source':tail['source_error'],'corrected_tail_remainder_variation':2*tail['remainder'],'corrected_tail_uniform_leading_truncation_variation':uniform_leading_truncation_variation}
 # Same empirical primitive policy, with actual envelope arithmetic scale recorded.
 arithmetic=floor*(1+abswork+math.fsum(components.values()))
 bound=math.fsum(components.values())+arithmetic
 return {'common_input_bound':bound,'components':components,'empirical_arithmetic_floor':arithmetic,
         'integrated_derivative_enclosures':{'T':list(tder),'D':list(dder),'q':list(qder),'h':list(hder)},
         'minimum_polynomial_modulus_margins':{'W':minW,'Q':minQ},'Nw':Nw,'angular_cells':cells,
         'cell_frequency_helicity_enclosures':2*Nw*cells,'physical_response_evaluations':0,
         'source_box_coverage':'Mean-value enclosure over entire rectangular source box; integrated derivatives enclose shared perturbations, not statistical correlation or corner coverage',
         'arithmetic_rigorous':False,'directed_rounding':False,'symbolic_coefficient_rounding':'Binary64 coefficient construction subject to unchanged empirical policy; not an exact arithmetic proof',
         'tail_source_semantics':'unchanged corrected_tail uniform C3/T source propagation; nonlinear tail excluded'}
