"""Independent source-functional Hessian evaluation; no Riccati imports.
Binary64 errors are empirical, never directed-rounding certificates.
"""
import cmath, math
from tail_sum import full_sum

def direction(j,N):
 if j%(N//4)==0:return ((1.,0.),(0.,1.),(-1.,0.),(0.,-1.))[j//(N//4)]
 t=2*math.pi*j/N;return math.cos(t),math.sin(t)

def kernel(z,D,p):
 O=cmath.sqrt(z*z+D*D)
 if O.real<0:O=-O
 Q=4*(z*z+D*D)+p*p
 if O.real<=0 or abs(Q)==0:raise ArithmeticError('reference causal/denominator exclusion')
 # Source functional derivative, distinct code/algebra from Riccati sectors.
 return ((8*z*z/(O*Q),4*z*p/(O*Q)),(-4*z*p/(O*Q),8*O/Q))

def strip_error(w,d,Dlo,Dhi,k,N):
 W=2*w*Dlo if w<Dlo else w*w+Dlo*Dlo;Q=4*W
 def r(a,b,c):return 2*c/(b+math.sqrt(b*b+4*a*c)) if a else c/b
 y1=r(d*d,2*d*(w+d),W/2) if d else 1.
 a=4*d*d+k*k;y2=r(a,2*a+8*w*d,Q/2) if a else 1.
 s=math.log1p(min(y1,y2));Z=w+d*math.cosh(s);O=math.sqrt(W/2);q=Q/2
 M=((2/w+8*Z*Z/(O*q),4*Z*k*math.cosh(s)/(O*q)),(4*Z*k*math.cosh(s)/(O*q),2/w+8*(Z+Dhi)/q))
 den=math.expm1(N*s) if N*s<700 else float('inf')
 return [[2*x/den for x in row] for row in M],s

def tail(source,radii,N,K):
 T,D,q,h=source;rt,rD,rq,rh=radii;tm=T-rt;dm=(abs(q)+rq)/2+abs(h)+rh;k=math.hypot(*K);S=dm+D+rD+k;r=1/(4*S);w=(2*N+1)*math.pi*tm
 if w*r<=1:raise ArithmeticError('independent reference tail exclusion')
 # Source C3 angular moments, no response tail call.
 C=[[0j,0j],[0j,0j]];C5=[[0j,0j],[0j,0j]]
 for l in (-1,1):
  d=q/2-l*h;weight=1+l/4
  C[0][0]+=weight*(d*d+k*k/4+3*D*D);C[1][1]+=weight*(d*d+k*k/4+D*D)
  C[0][1]+=weight*1j*d*K[0];C[1][0]-=weight*1j*d*K[0]
  angularcross=d*d*(3*K[0]**2+K[1]**2)
  C5[0][0]+=weight*(-.75*d**4-9*D*D*d*d-3*angularcross/8-3.75*D**4-.625*D*D*k*k-3*k**4/64)
  C5[1][1]+=weight*(-.75*d**4-3*D*D*d*d-3*angularcross/8-.75*D**4-.375*D*D*k*k-3*k**4/64)
  ap5=-1j*(1.5*d**3*K[0]+3*D*D*d*K[0]+.375*d*K[0]*k*k)
  C5[0][1]+=weight*ap5;C5[1][0]-=weight*ap5
 s3=full_sum(T,N,3);s5=full_sum(T,N,5)
 sum3=s3['value'];sum5=s5['value'];sum3error=s3['error'];sum5error=s5['error']
 sum7=1/w**7+1/(12*math.pi*tm*w**6)
 # On |1/w|=r each counterterm-subtracted entry <=16r;
 # inversion pair is odd, x term cancels; after independently angular-integrated C3/C5, remainder begins x7.
 rem=math.pi*(T+rt)*2*16/r**6*sum7/(1-1/(w*r)**2)
 corr=[[math.pi*T*(C[i][j]*sum3+C5[i][j]*sum5) for j in range(2)] for i in range(2)]
 err=[[rem+math.pi*T*(abs(C[i][j])*sum3error+abs(C5[i][j])*sum5error) for j in range(2)] for i in range(2)]
 return corr,err,{'radius':r,'remainder_entry_majorant':rem,'leading_sum_error':sum3error,'leading_terms':'infinite_EM_B6','C3_sum':s3,'C5_sum':s5,'orders':[3,5]}

def evaluate(source,K,Ntheta,Nw,radii,guard,floor=1e-12):
 T,D,q,h=source;rt,rD,rq,rh=radii;tm=T-rt;dm=(abs(q)+rq)/2+abs(h)+rh
 if tm<=0 or D-rD<=0:raise ArithmeticError('reference box exclusion')
 dirs=[direction(j,Ntheta) for j in range(Ntheta)];rows=[[] for _ in range(4)];absolute=[0.]*4;angular=[[0.,0.],[0.,0.]]
 for n in range(Nw):
  guard();w=(2*n+1)*math.pi*T;terms=[[] for _ in range(4)]
  for l in (-1,1):
   weight=1+l/4;d=q/2-l*h
   for c,s in dirs:
    G=kernel(complex(w,d*c),D,K[0]*c+K[1]*s)
    for i in range(2):
     for j in range(2):terms[2*i+j].append(weight*((2/w if i==j else 0)-G[i][j]))
  for ij,t in enumerate(terms):
   rows[ij].append(complex(math.fsum(x.real for x in t),math.fsum(x.imag for x in t))/Ntheta)
   absolute[ij]+=math.fsum(abs(x) for x in t)/Ntheta
  b,_=strip_error((2*n+1)*math.pi*tm,dm,D-rD,D+rD,math.hypot(*K),Ntheta)
  for i in range(2):
   for j in range(2):angular[i][j]+=math.pi*(T+rt)*2*b[i][j]
 raw=[[((2*math.log(T) if i==j else 0)+math.pi*T*complex(math.fsum(x.real for x in rows[2*i+j]),math.fsum(x.imag for x in rows[2*i+j]))) for j in range(2)] for i in range(2)]
 correction,tailerror,td=tail(source,radii,Nw,K)
 M=[[raw[i][j]+correction[i][j] for j in range(2)] for i in range(2)]
 primitive=[[floor*(1+math.pi*T*absolute[2*i+j]+abs(M[i][j])+abs(2*math.log(T))) for j in range(2)] for i in range(2)]
 errors=[[angular[i][j]+tailerror[i][j]+primitive[i][j] for j in range(2)] for i in range(2)]
 return {'raw_matrix':raw,'matrix':M,'entry_errors':errors,'angular_entry_errors':angular,'tail_entry_errors':tailerror,'primitive_entry_errors':primitive,'tail_details':td,'evaluations':2*Ntheta*Nw,'arithmetic_rigorous':False}

def witness_error(errors,v):return math.fsum(abs(v[i])*abs(v[j])*errors[i][j] for i in range(2) for j in range(2))

def pack(matrix):return [[[complex(x).real,complex(x).imag] for x in row] for row in matrix]
