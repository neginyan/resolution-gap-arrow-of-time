# Independent direct-sum check: nodes from each Gaussian component, pushed by the exact twist flow (RK4),
# blurred by an exact Gaussian kernel on local patches, entropy by quadrature, Sigma by central time difference.
import json,sys,numpy as np
obs=sys.argv[1]; sys.path.insert(0,obs)
from i3_flows import unpack, scales
sig=float(sys.argv[2]); fn=sys.argv[3]
d=json.load(open(obs+'/results/'+fn))
best=max(d,key=lambda r:r['gstar']); K=best['K']; x=np.array(best['x'])
sv0,sw0=scales('twist',sig); ws,mus,Ss=unpack(x,K,'twist',sv0,sw0)
print('K',K,'gstar(lab, sigma=1e-3)',best['gstar'],flush=True)
gh_x,gh_w=np.polynomial.hermite_e.hermegauss(7); gh_w=gh_w/gh_w.sum()
X=[];Wt=[]
for w,m,S in zip(ws,mus,Ss):
    l,U=np.linalg.eigh(np.array(S)); a=np.sqrt(l[1]); b=np.sqrt(max(l[0],0))
    h=sig/4; n=int(np.ceil(16*a/h)); u=np.linspace(-8*a,8*a,n); wu=np.exp(-.5*(u/a)**2); wu/=wu.sum()
    for xv,wv in zip(gh_x,gh_w):
        pts=np.array(m)[None,:]+np.outer(u,U[:,1])+xv*b*U[:,0][None,:]
        X.append(pts); Wt.append(w*wv*wu)
X=np.vstack(X); Wt=np.concatenate(Wt); print('nodes',len(X),'mass',Wt.sum(),flush=True)
def f(P):
    den=1+P[:,0]**2+P[:,1]**2; return np.stack([P[:,1]/den,-P[:,0]/den],1)
def push(P,t,steps=4):
    hh=t/steps
    for _ in range(steps):
        k1=f(P);k2=f(P+hh/2*k1);k3=f(P+hh/2*k2);k4=f(P+hh*k3);P=P+hh/6*(k1+2*k2+2*k3+k4)
    return P
dx=sig/float(sys.argv[4]); pad=9*sig
lo=X.min(0)-pad-0.01;hi=X.max(0)+pad+0.01
gq=np.arange(lo[0],hi[0],dx);gp=np.arange(lo[1],hi[1],dx); print('grid',len(gq),len(gp),flush=True)
R_=int(np.ceil(7*sig/dx)); off=np.arange(-R_,R_+1)
def density(P,grad=False):
    p=np.zeros((len(gq),len(gp))); g1=np.zeros_like(p) if grad else None; g2=np.zeros_like(p) if grad else None
    iq0=np.round((P[:,0]-lo[0])/dx).astype(int); ip0=np.round((P[:,1]-lo[1])/dx).astype(int)
    for s in range(0,len(P),4000):
        sl=slice(s,s+4000)
        IQ=(iq0[sl,None]+off[None,:]); IP=(ip0[sl,None]+off[None,:])
        dq=(lo[0]+IQ*dx-P[sl,0,None])/sig; dp=(lo[1]+IP*dx-P[sl,1,None])/sig
        eq=np.exp(-.5*dq**2); ep=np.exp(-.5*dp**2)
        Kv=Wt[sl,None,None]*eq[:,:,None]*ep[:,None,:]/(2*np.pi*sig**2)
        II=np.broadcast_to(IQ[:,:,None],Kv.shape); JJ=np.broadcast_to(IP[:,None,:],Kv.shape)
        np.add.at(p,(II.ravel(),JJ.ravel()),Kv.ravel())
        if grad:
            np.add.at(g1,(II.ravel(),JJ.ravel()),(Kv*(-dq[:,:,None]/sig)).ravel())
            np.add.at(g2,(II.ravel(),JJ.ravel()),(Kv*(-dp[:,None,:]/sig)).ravel())
    return p,g1,g2
def Hent(p):
    m=p>0; return -np.sum(p[m]*np.log(p[m]))*dx*dx
dt=float(sys.argv[5])
p0,g1,g2=density(X,True); m=p0>1e-300
trF=np.sum((g1[m]**2+g2[m]**2)/p0[m])*dx*dx; print('mass on grid',p0.sum()*dx*dx,flush=True)
Hp=Hent(density(push(X,dt))[0]); Hm=Hent(density(push(X,-dt))[0])
Sig=(Hp-Hm)/(2*dt); R=-Sig/(sig**2*trF)
print('my R',repr(R),' (kappa-R)/sigma',(0.25-R)/sig,flush=True)
