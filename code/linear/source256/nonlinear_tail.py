"""Uniform finite-beta asymptotic residual majorant, NO physics or target release.
Smooth finite Fourier source, real-positive omega only. No 1/omega Cauchy claim.
"""
import math


def uniform_residual_majorant(w,Dmax,dmax,Kmax,eta_l1,beta_disk,J=6):
    G=Dmax+beta_disk*eta_l1
    A=[0.,G/2]
    for m in range(1,J):
        A.append((m*Kmax+2*dmax)*A[m]/2+G/2*math.fsum(A[j]*A[m-j] for j in range(1,m)))
    profile=math.fsum(A[j]/w**j for j in range(1,J+1))
    rho=(J*Kmax+2*dmax)*A[J]/w**J+G*math.fsum(
        A[j]*A[k]/w**(j+k) for j in range(1,J+1) for k in range(1,J+1) if j+k>=J)
    L=1/(2*w);a=1-2*L*G*profile;disc=a*a-4*L*L*G*rho
    if a<=0 or disc<0:return {'pass':False,'reason':'asymptotic residual ball fails','rho':rho,'profile':profile}
    e=2*L*rho/(a+math.sqrt(disc));R=profile+e;lip=2*L*G*R
    if R>=1 or lip>=1:return {'pass':False,'reason':'causal unit disk/contraction bound fails','radius':e,'R':R}
    force_error=2*(1+R*R)/(1-R*R)**2*e
    return {'pass':True,'J':J,'coefficient_majorants':A,'residual_bound':rho,
            'profile_error':e,'profile_l1_upper':R,'lipschitz':lip,'f_error':force_error,
            'scope':'uniform beta-disk finite-asymptotic profile/error majorant; coefficients/normalization/path and cutoff integration remain separate'}


def correlated_divided_error(tail,beta_disk,epsilon_max):
    if not tail['pass'] or beta_disk<=epsilon_max:raise ValueError('valid analytic beta disk required')
    return tail['f_error']/(beta_disk-epsilon_max)
