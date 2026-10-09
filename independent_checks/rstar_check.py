import numpy as np
from scipy.optimize import minimize
# (1) pendulum: zero-width segment, V-comp alpha>=0, W-comp d>=alpha; S_qq=(d-alpha)^2/2
def Rpend(x,s):
    a,d=abs(x[0]),abs(x[1])
    lam=(1+np.exp(-(d-a)**2/4))/2
    return lam*(d*d-a*a)/(d*d+a*a+2*s*s)
# (2) quadratic stretch a=kappa-b q^2
def Rquad(x,s,k,b):
    a,d=abs(x[0]),abs(x[1])
    return (k-b*(d-a)**2/2)*(d*d-a*a)/(d*d+a*a+2*s*s)
def best(F,args,starts):
    return max(-minimize(lambda x:-F(x,*args),x0,method='Nelder-Mead',options=dict(xatol=1e-14,fatol=1e-16,maxiter=40000)).fun for x0 in starts)
for s in [1e-3,1e-2,0.1,0.3]:
    st=[[0.3*s,2*np.sqrt(s)],[0.1*s,1.5*np.sqrt(s)],[s,3*np.sqrt(s)]]
    Rs=best(Rpend,(s,),st)
    ser=s-7/8*s**2+65/96*s**3-61/128*s**4+3151/10240*s**5-8483/46080*s**6
    print(f"pendulum sigma={s}: 1-R*={1-Rs:.12f}  series(6)={ser:.12f}  lead 2sqrt(2*beta*kappa)*s={1.0*s:.3g}")
for k,b,name in [(1,1.5,'double well'),(0.5,1.5,'quartic local')]:
    for s in [1e-3,1e-2]:
        sq=np.sqrt(s)
        Rs=best(Rquad,(s,k,b),[[0.3*s,2*sq*(k/b)**.25],[0.1*s,1.5*sq*(k/b)**.25],[s,3*sq*(k/b)**.25]])
        eps=s*np.sqrt(b/k); Psi=2*eps-5/2*eps**2+7/4*eps**3-1/8*eps**4-53/64*eps**5
        beta=b/2
        print(f"{name} sigma={s}: (k-R*)/s={(k-Rs)/s:.7f} pred 2sqrt(2 beta k)={2*np.sqrt(2*beta*k):.7f}; 1-R*/k={1-Rs/k:.12f} Psi={Psi:.12f}")
