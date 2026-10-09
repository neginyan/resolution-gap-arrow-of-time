"""Independent check (written from scratch, no Stein/Tweedie/posterior formulas).
R = -Sigma/(sigma^2 trF), Sigma = d h(Y)/dt with window fixed.
Method: rho on a grid, transported exactly by Liouville rho_t(x)=rho0(phi_{-t}x) (pendulum, RK4),
blurred by FFT Gaussian, entropy by quadrature, time derivative by central difference."""
import numpy as np, sys
V=np.array([1,1])/np.sqrt(2); W=np.array([1,-1])/np.sqrt(2)
def shape(n,l):return n*n*np.outer(V,V)+l*l*np.outer(W,W)
def rho0(Q,P,ws,mus,Ss):
    out=0
    for w,m,S in zip(ws,mus,Ss):
        Si=np.linalg.inv(S);dq,dp=Q-m[0],P-m[1]
        out=out+w*np.exp(-.5*(Si[0,0]*dq*dq+2*Si[0,1]*dq*dp+Si[1,1]*dp*dp))/(2*np.pi*np.sqrt(np.linalg.det(S)))
    return out
def back(Q,P,t,steps=4):
    h=-t/steps
    f=lambda q,p:(p,-np.sin(q))
    for _ in range(steps):
        k1=f(Q,P);k2=f(Q+h/2*k1[0],P+h/2*k1[1]);k3=f(Q+h/2*k2[0],P+h/2*k2[1]);k4=f(Q+h*k3[0],P+h*k3[1])
        Q,P=Q+h/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),P+h/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    return Q,P
def run(ws,mus,Ss,sig,L=6.5,N=2048,dt=2e-3):
    x=np.linspace(-L,L,N,endpoint=False);dx=x[1]-x[0]
    Q,P=np.meshgrid(x,x,indexing='ij')
    k=2*np.pi*np.fft.fftfreq(N,dx);KQ,KP=np.meshgrid(k,k,indexing='ij')
    G=np.exp(-.5*sig**2*(KQ**2+KP**2))
    def pY(t):
        r=rho0(*back(Q,P,t),ws,mus,Ss)
        return np.clip(np.real(np.fft.ifft2(np.fft.fft2(r)*G)),1e-300,None)
    def H(p):return -np.sum(p*np.log(p))*dx*dx
    Sig=(H(pY(dt))-H(pY(-dt)))/(2*dt)
    p=pY(0.0)
    gq=np.real(np.fft.ifft2(1j*KQ*np.fft.fft2(p)));gp=np.real(np.fft.ifft2(1j*KP*np.fft.fft2(p)))
    m=p>1e-200
    trF=np.sum((gq[m]**2+gp[m]**2)/p[m])*dx*dx
    return -Sig/(sig**2*trF), np.sum(p)*dx*dx
if __name__=='__main__':
    n1=0.09
    for ratio,delta in [(4.4,0.0),(1.0,0.0),(8.0,0.0),(4.4,1.0)]:
        u=W
        ws=[.5,.5];mus=[-delta/2*u,delta/2*u];Ss=[shape(n1,1.9),shape(ratio*n1,0.63*1.9)]
        R,mass=run(ws,mus,Ss,0.1)
        print(f"ratio={ratio} delta={delta}  R_indep={R:.5f}  mass={mass:.6f}",flush=True)
