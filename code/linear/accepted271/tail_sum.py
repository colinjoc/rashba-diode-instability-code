"""Euler-Maclaurin p>1 positive power sum; empirical binary64, analytic remainder."""
import math

def full_sum(T,N,p=3,floor=1e-12):
 if T<=0 or N<1 or p not in (3,5):raise ValueError('tail domain')
 a=N+.5
 rising=lambda n:math.prod(p+j for j in range(n))
 terms=[a**(1-p)/(p-1),.5*a**-p,rising(1)/12*a**(-p-1),-rising(3)/720*a**(-p-3),rising(5)/30240*a**(-p-5)]
 scale=(2*math.pi*T)**-p
 value=math.fsum(terms)*scale
 # EM through B8 has the explicit B8 boundary term and periodic B8 remainder.
 # We omit that boundary term from value, and bound BOTH by triangle inequality.
 # |periodic B8|<=1/30; integral |f^(8)|=(p)_7*a^(-p-7).
 analytic=2*rising(7)/1209600*a**(-p-7)*scale
 empirical=floor*math.fsum(abs(x) for x in terms)*scale
 return dict(value=value,analytic_error=analytic,empirical_error=empirical,error=analytic+empirical,lower=value-analytic-empirical,upper=value+analytic+empirical,arithmetic_rigorous=False)
