"""Known source-compatible radial-integrated kernel. No executable target search."""
import numpy as np
from scipy.special import polygamma


def tails(T, M):
    return (-polygamma(2,M+.5)/(2*(2*np.pi*T)**3),
            -polygamma(4,M+.5)/(24*(2*np.pi*T)**5))


def geometry(T,D,q,h,K,nt,nw):
    theta=2*np.pi*np.arange(nt)/nt
    n=np.array([-1.,1.])[:,None,None]
    weights=1+n/4
    delta=(q/2-n*h)*np.cos(theta)[None,None,:]
    omega=((2*np.arange(nw)+1)*np.pi*T)[None,:,None]
    z=omega+1j*delta
    om=np.sqrt(z*z+D*D)
    p=(K[0]*np.cos(theta)+K[1]*np.sin(theta))[None,None,:]
    return weights,delta,omega,z,om,p


def uniform(T,D,q,h,nt,nw):
    weights,delta,omega,z,om,p=geometry(T,D,q,h,(0.,0.),nt,nw)
    c3,c5=tails(T,nw)
    # Weighted angular asymptotics of the canonical gap residual / Delta.
    t3=np.sum(np.mean(weights*(delta*delta+D*D/2),axis=2))
    t5=np.sum(np.mean(weights*(delta**4+3*D*D*delta*delta+3*D**4/8),axis=2))
    residual=np.log(T)+np.pi*T*np.sum(np.mean(1/omega-weights*np.real(1/om),axis=2))+np.pi*T*(t3*c3-t5*c5)
    costheta=np.cos(2*np.pi*np.arange(nt)/nt)
    current=2*np.pi*T*np.sum(np.mean(weights*costheta*np.imag(z/om),axis=2))
    current3=np.sum(np.mean(weights*costheta*D*D*delta,axis=2))
    current5=np.sum(np.mean(weights*costheta*(2*D*D*delta**3+1.5*D**4*delta),axis=2))
    current+=2*np.pi*T*(current3*c3-current5*c5)
    f3=np.sum(np.mean(weights*(D*D*delta*delta/2+D**4/8),axis=2))
    f5=np.sum(np.mean(weights*(D*D*delta**4/2+3*D**4*delta*delta/4+D**6/16),axis=2))
    free=D*D*np.log(T)+2*np.pi*T*np.sum(np.mean(D*D/(2*omega)-weights*np.real(D*D/(om+z)),axis=2))+2*np.pi*T*(f3*c3-f5*c5)
    return float(residual),float(current),float(free)


def kernel(T,D,q,h,K,nt,nw):
    weights,delta,omega,z,om,p=geometry(T,D,q,h,K,nt,nw)
    denominator=4*om*om+p*p
    c3,c5=tails(T,nw)
    ca=1.5*D*D+p*p/4
    cp=.5*D*D+p*p/4
    da=15*D**4/8+5*D*D*p*p/8+p**4/16
    dp=3*D**4/8+3*D*D*p*p/8+p**4/16
    diag=[]
    for response,c,d in [(4*z*z/(om*denominator),ca,da),(4*om/denominator,cp,dp)]:
        tail3=np.sum(np.mean(weights*(delta*delta+c),axis=2))
        tail5=np.sum(np.mean(weights*(delta**4+6*c*delta*delta+d),axis=2))
        diag.append(float(2*np.log(T)+2*np.pi*T*np.sum(np.mean(1/omega-weights*np.real(response),axis=2))+2*np.pi*T*(tail3*c3-tail5*c5)))
    cross=np.sum(np.mean(weights*p*np.imag(z/(om*denominator)),axis=2))
    cross3=np.sum(np.mean(weights*(-p*delta/2),axis=2))
    cross5=np.sum(np.mean(weights*p*(delta**3+ca*delta),axis=2))
    off=-4j*np.pi*T*(cross+cross3*c3+cross5*c5)
    return np.array([[diag[0],off],[np.conj(off),diag[1]]])
