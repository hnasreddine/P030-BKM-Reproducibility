import sys,math,time,os
import mpmath as mp
import numpy as np
sys.path.insert(0,os.path.dirname(__file__))
import p030lib as L
mp.mp.dps=55
lam=mp.mpf(sys.argv[1]);c=mp.mpf(sys.argv[2]);N=int(sys.argv[3]);out=sys.argv[4]
gf=mp.mpf(sys.argv[5]) if len(sys.argv)>5 else mp.mpf(3)
b=L.bath_poly(2,gs=1,gf=float(gf))
p=[mp.mpf(str(float(v))) for v in b.pcoef]; eta=mp.mpf('1.8')
R=mp.mpf('2'); T=c*mp.log(1/lam)/(lam*lam);u=R*T
K=int(float(u)+32*float(u)**(1/3)+150)
start=K+100; jp=mp.mpf(0);jc=mp.mpf('1e-400');vals=[None]*(K+1);norm=mp.mpf(0)
for k in range(start,0,-1):
 jm=2*k/u*jc-jp;jp,jc=jc,jm
 if k-1<=K:vals[k-1]=jc
 if (k-1)%2==0:norm+=(1 if k==1 else 2)*jc

def bump(t):
 return mp.exp(1-1/(1-t*t)) if abs(t)<1 else mp.mpf(0)

dx=mp.mpf(2)/(N-1)
t0=time.time()
with open(out,'w') as f:
 f.write('%d %d %s %s %s\n'%(N,K,mp.nstr(R,52),mp.nstr(T,52),mp.nstr(lam,52)))
 for i in range(N):
  x=-1+i*dx; ph=bump(x);te=bump((abs(x)-mp.mpf('.5'))/mp.mpf('.4'))
  theta=bump((x-mp.mpf('.5'))/mp.mpf('.4'))
  pv=p[0]+p[1]*x+p[2]*x*x
  j00=ph/mp.pi+2*te; j11=gf*ph/mp.pi+2*te; j01=eta*theta*pv
  sd=mp.sqrt(max(j00*j11-j01*j01,mp.mpf(0)))
  den=mp.sqrt(j00+j11+2*sd) if j00+j11+2*sd>0 else mp.mpf(1)
  a=mp.sqrt(dx)*(j00+sd)/den
  z=mp.sqrt(dx)*j01/den
  d=mp.sqrt(dx)*(j11+sd)/den
  f.write(' '.join(mp.nstr(t,42) for t in [x,a,z,d])+'\n')
 for v in vals:f.write(mp.nstr(v/norm,43)+'\n')
print('HP PREP',str(lam),str(c),'N',N,'K',K,'sec',time.time()-t0,'Jtail',mp.nstr(vals[-1]/norm,8),flush=True)
