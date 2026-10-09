# my own best single Gaussian R*(sigma), from the proven Gaussian law R = lambda * phi
import numpy as np
from scipy.optimize import minimize
V=np.array([1,1])/np.sqrt(2);W=np.array([1,-1])/np.sqrt(2)
def R1(x,sig):
    mq,a,b,c=x
    L=np.array([[a,0],[b,c]]);S=L@L.T
    lam=(1-np.cos(mq)*np.exp(-S[0,0]/2))/2
    T=S+sig**2*np.eye(2)
    Tv=V@T@V;Tw=W@T@W
    return lam*(Tw-Tv)/(Tw+Tv)
rng=np.random.default_rng(1)
for sig in [0.1,0.3]:
    best=-1
    for k in range(400):
        x0=[np.pi+rng.normal(0,.3),rng.normal(0,1),rng.normal(0,1),rng.normal(0,1)*rng.choice([1,1e-2,1e-4])]
        r=minimize(lambda x:-R1(x,sig),x0,method='Nelder-Mead',options=dict(xatol=1e-12,fatol=1e-15,maxiter=20000))
        best=max(best,-r.fun)
    print('sigma',sig,'my R* =',repr(best))
