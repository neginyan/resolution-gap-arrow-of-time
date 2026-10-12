# Independent check with one momentum direction u = (0,0,cos phi, sin phi) of infinite variance (exact limit):
# F = (1/sigma) Q (L L^T + sigma I)^-1 Q^T on the complement Q; law 1 R = tr(F SymJbar)/tr F.
import numpy as np, sys
from scipy.optimize import minimize
eps=float(sys.argv[1]); sig=float(sys.argv[2]); ntry=int(sys.argv[3])
def symJ(mq,Sq):
    E=lambda a: np.cos(a@mq)*np.exp(-0.5*a@Sq@a)
    c1,c2,c12=E(np.array([1.,0])),E(np.array([0,1.])),E(np.array([1.,-1]))
    A=-np.array([[c1+eps*c12,-eps*c12],[-eps*c12,c2+eps*c12]])
    J=np.zeros((4,4)); J[:2,2:]=np.eye(2); J[2:,:2]=A; return 0.5*(J+J.T)
kap=(1+eps+np.sqrt(1+eps*eps))/2
def gap(x):
    mq=np.array([np.pi,0])+sig*x[:2]; ph=x[2]
    Q=np.zeros((4,3)); Q[0,0]=1; Q[1,1]=1; Q[2,2]=-np.sin(ph); Q[3,2]=np.cos(ph)
    L=np.zeros((3,3)); L[np.tril_indices(3)]=x[3:]; LL=L@L.T
    G=Q@np.linalg.inv(LL+sig*np.eye(3))@Q.T; Sq=sig*LL[:2,:2]
    return (kap-np.trace(G@symJ(mq,Sq))/np.trace(G))/sig
rng=np.random.default_rng(2); best=None
for k in range(ntry):
    x0=np.concatenate([[0,0],[rng.uniform(0,np.pi)],rng.normal(0,1.5,6)])
    r=minimize(gap,x0,method='Nelder-Mead',options=dict(maxiter=60000,xatol=1e-10,fatol=1e-13))
    r=minimize(gap,r.x,method='BFGS',options=dict(gtol=1e-11))
    if best is None or r.fun<best.fun: best=r
print(f"eps {eps} sigma {sig}: min (kappa-R*)/sigma = {best.fun:.6f}  (phi {best.x[2]%np.pi:.4f})")
