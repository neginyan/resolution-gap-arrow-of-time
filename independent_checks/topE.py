import json,sys,numpy as np
sys.path.insert(0,'.')
import indep
d=json.load(open(sys.argv[1]))
def walk(x):
    if isinstance(x,dict):
        if 'ws' in x and 'R' in x: yield x
        for v in x.values(): yield from walk(v)
    elif isinstance(x,list):
        for v in x: yield from walk(v)
b=max(walk(d),key=lambda r:r['R'])
ws=b['ws'];mus=[np.array(m) for m in b['mus']];Ss=[np.array(S) for S in b['Ss']];sig=b['sigma']
print('lab R',b['R'],flush=True)
lo=np.min([m-3.8*np.sqrt(np.diag(S)) for m,S in zip(mus,Ss)],0);hi=np.max([m+3.8*np.sqrt(np.diag(S)) for m,S in zip(mus,Ss)],0)
c=(lo+hi)/2;L=float(np.max(hi-lo)/2)+0.6;print('center',c,'L',L,flush=True)
mus=[m-c for m in mus]
def back(Q,P,t,steps=4):
    h=-t/steps;f=lambda q,p:(p,-np.sin(q+c[0]))
    for _ in range(steps):
        k1=f(Q,P);k2=f(Q+h/2*k1[0],P+h/2*k1[1]);k3=f(Q+h/2*k2[0],P+h/2*k2[1]);k4=f(Q+h*k3[0],P+h*k3[1])
        Q,P=Q+h/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),P+h/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    return Q,P
# p-shift: force does not depend on p, and flow q'=p uses absolute p -> shift p back
def back2(Q,P,t,steps=4):
    Qa,Pa=back(Q,P+c[1],t,steps);return Qa,Pa-c[1]
indep.back=back2
for N in [int(sys.argv[2])]:
    R,m=indep.run(ws,mus,Ss,sig,L=L,N=N,dt=2e-4);print(N,'spacing',2*L/N,'indep R',round(R,6),'mass',round(m,6),flush=True)
