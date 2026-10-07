"""Two independently evaluated, analytically differentiated current conventions.
No finite differences, observed-discrepancy tolerances or curvature transfers.
"""
import cmath,math

def choose(a,k):
 x=1.
 for j in range(k):x*= (a-j)/(j+1)
 return x

def moment(k):return 0. if k%2 else math.comb(k,k//2)/4**(k//2)

def current_tail(T,D,q,h,N,order=9):
 C={s:0. for s in (3,5,7,9)}; CM={s:0. for s in C}
 for l in (-1,1):
  d=q/2-l*h;weight=1+l/4
  for m in range(1,5):
   for k in range(1,order-2*m+1,2):
    power=2*m+k
    # j/(2piT) = <cos theta Im(z/O)>.
    term=weight*(choose(-.5,m)*D**(2*m)*choose(-2*m,k)*(1j*d)**k).imag*moment(k+1)
    C[power]+=term;CM[power]+=m*term/(D*D)
 correction=0.;sum_error=0.;mc=0.;me=0.;stop=N+65536
 for power,c in C.items():
  sm=math.fsum(((2*n+1)*math.pi*T)**(-power) for n in range(N,stop))
  er=((2*stop-1)*math.pi*T)**(1-power)/(2*math.pi*T*(power-1))
  correction+=2*math.pi*T*c*sm;sum_error+=2*math.pi*T*abs(c)*er
  mc+=2*math.pi*T*CM[power]*sm;me+=2*math.pi*T*abs(CM[power])*er
 return correction,sum_error,C,mc,me

def evaluate(p,ref,guard):
 f=p['frozen'];T,D,q,h=f['x'];rt,rD,rq,rh=p['source_box_radii'];tm=T-rt;dm=(abs(q)+rq)/2+abs(h)+rh;Dhi=D+rD
 Nt,Nw=p['levels'][0][-1],p['levels'][1][-1];printed=[];functional=[];mixed1=[];mixed2=[];absolute=0.;angular=0.;mixed_angular=0.
 for n in range(Nw):
  guard();w=(2*n+1)*math.pi*T;pj=[];fj=[];mj=[];nj=[]
  for l in (-1,1):
   weight=1+l/4;d=q/2-l*h
   for j in range(Nt):
    c,s=ref.direction(j,Nt);z=complex(w,d*c);O=cmath.sqrt(z*z+D*D)
    if O.real<0:O=-O
    if O.real<=0 or abs(O+z)==0:raise ArithmeticError('normalization control causal denominator')
    # Printed Eilenberger current in N0 units, stable D2 form.
    pj.append(-weight*c*D*D*(1/(O*(O+z))).imag)
    # Independent analytic derivative of F=D2lnT+piT sum[D2/w-2 Re(O-z)].
    zq=1j*c/2;Oq=z*zq/O;fj.append(weight*(-2*(Oq-zq).real))
    # lim M_AP/(iK) from G source-functional derivative vs F_qD/D.
    mj.append(-weight*c*(z/O**3).imag)
    # Differentiated stable 1-z/O = D2/[O(O+z)] identity independently.
    L=O*(O+z);OD=D/O;LD=OD*(2*O+z)
    nj.append(-weight*c*(2/L-D*LD/(L*L)).imag)
  printed.append(math.fsum(pj)/Nt);functional.append(math.fsum(fj)/Nt)
  mixed1.append(math.fsum(mj)/Nt);mixed2.append(math.fsum(nj)/Nt);absolute+=math.fsum(abs(x) for x in pj+fj+mj+nj)/Nt
  _,s=ref.strip_error((2*n+1)*math.pi*tm,dm,D-rD,Dhi,0.,Nt)
  den=math.expm1(Nt*s) if Nt*s<700 else float('inf');ww=(2*n+1)*math.pi*tm
  Omin=math.sqrt((2*ww*(D-rD) if ww<D-rD else ww*ww+(D-rD)**2)/2);Z=ww+dm*math.cosh(s)
  angular+=2*math.pi*(T+rt)*2*2*math.cosh(s)*(1+Z/Omin)/den
  mixed_angular+=math.pi*(T+rt)*2*2*math.cosh(s)*Z/Omin**3/den
 jc=2*math.pi*T*math.fsum(printed);fd=2*math.pi*T*math.fsum(functional)
 # functional is Fq; printed j=2Fq. The common 2piT prefactor applies both arrays.
 mc=math.pi*T*math.fsum(mixed1);md=math.pi*T*math.fsum(mixed2)
 correction,leading_error,coeff,mixed_correction,mixed_leading_error=current_tail(T,D,q,h,Nw)
 r=1/(4*(dm+Dhi));wo=(2*Nw+1)*math.pi*tm
 if wo*r<=1:raise ArithmeticError('current tail Cauchy exclusion')
 sum11=wo**-11+wo**-10/(20*math.pi*tm)
 # On inverse-frequency radius r, |z/O|<4; cosine<=1 real.
 # angular/current inversion retains odd powers; remove through x9.
 rem=2*math.pi*(T+rt)*2*4/r**11*sum11/(1-1/(wo*r)**2)
 current_error=angular+leading_error+rem+p['arithmetic_floor']*(1+absolute+abs(jc)+abs(fd))
 # Differentiate the same explicit odd tail coefficients analytically.
 # Cauchy derivative in D within radius Dlo/4 supplies a separate bound;
 # use enlarged gap3Dhi/2 for the inverse-frequency analytic disk.
 rm=1/(4*(dm+1.5*Dhi));mixed_rem=2*math.pi*(T+rt)*2*32/(D-rD)**2/rm**11*sum11/(1-1/(wo*rm)**2)
 mixed_error=mixed_angular+mixed_leading_error+mixed_rem+p['arithmetic_floor']*(1+absolute+abs(mc)+abs(md))
 mc+=mixed_correction;md+=mixed_correction
 tolerance=p['normalization_tolerance'];scale=p['source_control_relative_tolerance']*max(1.,abs(jc+correction),abs(fd+correction))
 return {'printed_j_over_N0':jc+correction,'analytic_2Fq_over_N0':fd+correction,'difference':abs(jc-fd),
  'current_error_each':current_error,'current_adequacy_threshold':scale,
  'current_pass':abs(jc-fd)<=tolerance+2*current_error and current_error<=scale,
  'mixed_source_derivative':mc,'mixed_functional_FqD_over_D':md,'mixed_difference':abs(mc-md),
  'mixed_error_each':mixed_error,'mixed_pass':abs(mc-md)<=tolerance+2*mixed_error and mixed_error<=scale,
  'mixed_scope':'Independent analytic differential identity at common finite regulator; continuum tail coverage conservative and recorded; this normalization guard contributes no scalar current-to-curvature error',
  'current_tail':{'correction':correction,'remainder':rem,'leading_sum_error':leading_error,'coefficients':coeff,'order':9},
  'normalization':'printed j/N0=2 partial_q(F/N0); DOS weights sum two; piT explicit',
  'domain':'Uniform source-box analytic angular/causal/tail bounds; comparison at identical frozen nominal source. Shared parameter variation cancels in normalization identity; no use of saved256 for these guard errors.',
  'tolerance_frozen':tolerance,'arithmetic_rigorous':False,'discrepancy_in_tolerance':False}
