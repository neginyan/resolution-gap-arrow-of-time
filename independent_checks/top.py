import json,sys,numpy as np
sys.path.insert(0,'.')
import indep
d=json.load(open(sys.argv[1]))
best=max(d,key=lambda r:r['R'])
ws=best['ws'];mus=[np.array(m) for m in best['mus']];Ss=[np.array(S) for S in best['Ss']]
sig=best['sigma'];print('lab R',best['R'],'K',best['K'],'sigma',sig,flush=True)
# centre grid on the state (indep grid is centred at 0): shift means
c=np.mean(mus,0)
mus=[m-c for m in mus]
# pendulum is not translation invariant in q -> shift the force instead
def back(Q,P,t,steps=4):
    h=-t/steps;f=lambda q,p:(p,-np.sin(q+c[0]))
    for _ in range(steps):
        k1=f(Q,P);k2=f(Q+h/2*k1[0],P+h/2*k1[1]);k3=f(Q+h/2*k2[0],P+h/2*k2[1]);k4=f(Q+h*k3[0],P+h*k3[1])
        Q,P=Q+h/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),P+h/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    return Q,P
indep.back=back
for L,N in [(4.0,3072),(4.0,4096)]:
    R,m=indep.run(ws,mus,Ss,sig,L=L,N=N,dt=5e-4);print(L,N,'indep R',round(R,6),'mass',round(m,6),flush=True)
