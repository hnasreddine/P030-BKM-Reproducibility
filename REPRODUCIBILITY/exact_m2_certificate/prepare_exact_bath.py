"""Prepare near-exact order-two polynomial bath for the stored cheb_dd runner.
The polynomial is defined analytically by orthogonality, numerical coefficients
are from check_orthogonal.py; do NOT interpret this as validated coefficient enclosure.
"""
from pathlib import Path
import mpmath as m, time
m.mp.dps=68
root=Path(str(Path(__file__).resolve().parent))
old=(Path(__file__).resolve().parent / 'input_005_5_8500.txt')
new=root/'input_005_5_8500_ideal_m2.txt'
p=list(map(m.mpf,(root/'orthogonal_coefs.txt').read_text().splitlines()))

def bump(t):
 return m.exp(1-1/(1-t*t)) if abs(t)<1 else m.mpf(0)
N=8500
with old.open() as f:
 header=f.readline().strip()
 for i in range(N):f.readline()
 bessel=f.read()
N0,K,R,T,lam=header.split();assert int(N0)==N
eta=m.mpf('1.8');gf=m.mpf('3');dx=m.mpf(2)/(N-1);tm=time.time()
with new.open('w') as f:
 f.write(header+'\n')
 for j in range(N):
  x=-1+j*dx
  ph=bump(x)
  te=bump((abs(x)-m.mpf('.5'))/m.mpf('.4'))
  th=bump((x-m.mpf('.5'))/m.mpf('.4'))
  v=p[0]+p[1]*x+p[2]*x*x
  j00=ph/m.pi+2*te
  j11=gf*ph/m.pi+2*te
  j01=eta*th*v
  det=j00*j11-j01*j01
  if det<0 and abs(det)>m.mpf('1e-50'):raise Exception((j,det))
  sd=m.sqrt(max(m.mpf(0),det))
  den=m.sqrt(j00+j11+2*sd) if j00+j11+2*sd>0 else m.mpf(1)
  a=m.sqrt(dx)*(j00+sd)/den
  b=m.sqrt(dx)*j01/den
  d=m.sqrt(dx)*(j11+sd)/den
  f.write(' '.join(m.nstr(z,42) for z in [x,a,b,d])+'\n')
  if j%2000==0:print('grid',j,'elapsed',time.time()-tm,flush=True)
 f.write(bessel)
print('saved',new,'bytes',new.stat().st_size,'elapsed',time.time()-tm,flush=True)
