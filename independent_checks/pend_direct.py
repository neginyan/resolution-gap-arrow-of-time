# Independent direct-sum check for pendulum mixtures (no posterior formulas, no Stein/Tweedie):
# nodes from each Gaussian component, pushed by the exact pendulum flow (RK4), blurred with the exact Gaussian kernel,
# entropy by quadrature on a grid, Sigma by central time difference; R = -Sigma/(sigma^2 trF).
import json,sys,numpy as np
d=json.load(open(sys.argv[1]))[int(sys.argv[2])]; ws,mus,Ss,sig=d['ws'],[np.array(m) for m in d['mus']],[np.array(S) for S in d['Ss']],d['sigma']
hfac=float(sys.argv[3]); dxfac=float(sys.argv[4]); dt=float(sys.argv[5])
gh_x,gh_w=np.polynomial.hermite_e.hermegauss(7); gh_w=gh_w/gh_w.sum()
X=[];Wt=[]
for w,m,S in zip(ws,mus,Ss):
    l,U=np.linalg.eigh(S); a=np.sqrt(l[1]); b=np.sqrt(max(l[0],0))
    h=sig/hfac; n=max(int(np.ceil(16*a/h)),3); u=np.linspace(-8*a,8*a,n); wu=np.exp(-.5*(u/a)**2); wu/=wu.sum()
    for xv,wv in zip(gh_x,gh_w):
        X.append(m[None,:]+np.outer(u,U[:,1])+xv*b*U[:,0][None,:]); Wt.append(w*wv*wu)
X=np.vstack(X); Wt=np.concatenate(Wt)
f=lambda P: np.stack([P[:,1],-np.sin(P[:,0])],1)
def push(P,t,steps=4):
    hh=t/steps
    for _ in range(steps):
        k1=f(P);k2=f(P+hh/2*k1);k3=f(P+hh/2*k2);k4=f(P+hh*k3);P=P+hh/6*(k1+2*k2+2*k3+k4)
    return P
dx=sig/dxfac; pad=8*sig
lo=X.min(0)-pad-3*dx; hi=X.max(0)+pad+3*dx
nq=int(np.ceil((hi[0]-lo[0])/dx)); npp=int(np.ceil((hi[1]-lo[1])/dx))
R_=int(np.ceil(7*sig/dx)); off=np.arange(-R_,R_+1)
def density(P,grad=False):
    p=np.zeros(nq*npp); g1=np.zeros_like(p) if grad else None; g2=np.zeros_like(p) if grad else None
    iq0=np.round((P[:,0]-lo[0])/dx).astype(int); ip0=np.round((P[:,1]-lo[1])/dx).astype(int)
    for s in range(0,len(P),3000):
        sl=slice(s,s+3000)
        IQ=iq0[sl,None]+off[None,:]; IP=ip0[sl,None]+off[None,:]
        dq=(lo[0]+IQ*dx-P[sl,0,None])/sig; dp=(lo[1]+IP*dx-P[sl,1,None])/sig
        eq=np.exp(-.5*dq**2); ep=np.exp(-.5*dp**2)
        Kv=Wt[sl,None,None]*eq[:,:,None]*ep[:,None,:]/(2*np.pi*sig**2)
        idx=(IQ[:,:,None]*npp+IP[:,None,:]).ravel()
        p+=np.bincount(idx,Kv.ravel(),minlength=len(p))
        if grad:
            g1+=np.bincount(idx,(Kv*(-dq[:,:,None]/sig)).ravel(),minlength=len(p))
            g2+=np.bincount(idx,(Kv*(-dp[:,None,:]/sig)).ravel(),minlength=len(p))
    return p,g1,g2
def Hent(p):
    m=p>0; return -np.sum(p[m]*np.log(p[m]))*dx*dx
p0,g1,g2=density(X,True); m=p0>1e-300
trF=np.sum((g1[m]**2+g2[m]**2)/p0[m])*dx*dx
Hp=Hent(density(push(X,dt))[0]); Hm=Hent(density(push(X,-dt))[0])
R=-(Hp-Hm)/(2*dt)/(sig**2*trF)
print(f"K {d['K']} sigma {sig} nodes {len(X)} grid {nq}x{npp} h=sig/{hfac} dx=sig/{dxfac} dt {dt}: mass {p0.sum()*dx*dx:.12f}  my (1-R)/sigma {(1-R)/sig:.6f}   lab {d['gap']:.6f}",flush=True)
