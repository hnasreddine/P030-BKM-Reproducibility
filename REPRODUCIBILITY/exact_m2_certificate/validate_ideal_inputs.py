from pathlib import Path
"""Validate decimal input tokens against analytic bath samples using mpmath.iv.
A separate DD parse error budget is required. For historical fixed SVD pcoef:
 p(x)=0.8406338664123079-4.503497377672212 x+5.200634288262571 x².
"""
import mpmath as mp, pathlib, time
mp.iv.dps=55
I=mp.iv.mpf
N=8500;dx=I(2)/(N-1);f=pathlib.Path(str(Path(__file__).resolve().parent / 'input_005_5_8500_ideal_m2.txt')).open()
head=f.readline().split();assert int(head[0])==N
lam=I('.05');eta=I('1.8');gf=I(3)
p=[I(s) for s in Path(__file__).resolve().parent.joinpath('orthogonal_coefs.txt').read_text().splitlines()]
def bump(t):
 if t.a >= I(1).a or t.b <= I(-1).b:return I(0)
 return mp.iv.exp(1-1/(1-t*t))
maxerr=I(0);sumx=I(0);sumg=I(0);mx={'x':I(0),'a':I(0),'b':I(0),'d':I(0)}; where={}
t0=time.time()
for j in range(N):
 x=I(-1)+j*dx
 if j==N-1:x=I(1)
 if j==0:x=I(-1)
 xv=float(-1+2*j/(N-1))
 xx=I(str(xv)) # only branch via rational/comparison: branch not near boundaries
 ph=bump(x)
 te=bump((abs(x)-I('.5'))/I('.4'))
 theta=bump((x-I('.5'))/I('.4')) if j>(N-1)*.55 and j<(N-1)*.95 else I(0)
 pv=p[0]+p[1]*x+p[2]*x*x
 j00=ph/mp.iv.pi+2*te
 j11=gf*ph/mp.iv.pi+2*te
 j01=eta*theta*pv
 rad=j00*j11-j01*j01
 if rad.a < I(0).a:
  if rad.b >=I(0).a and rad.a > I('-1e-100'):
   rad=I([0,rad.b])
  else:raise RuntimeError(f'negative bath determinant at {j}: {rad}')
 sd=mp.iv.sqrt(rad)
 den=mp.iv.sqrt(j00+j11+2*sd) if j00.b > I(0).a or j11.b>I(0).a else I(1)
 fac=mp.iv.sqrt(dx)
 predicted=[x,fac*(j00+sd)/den,fac*j01/den,fac*(j11+sd)/den]
 if any('inf' in str(z) for z in predicted):
  print('INF AT',j,'x',x,'ph',ph,'te',te,'theta',theta,'j00',j00,'j11',j11,'j01',j01,'rad',rad,'sd',sd,'den',den,flush=True);raise Exception('inf')
 nums=[I(v) for v in f.readline().split()]
 for tag,num,exact in zip(['x','a','b','d'],nums,predicted):
  err=abs(exact-num)
  if err.b>mx[tag].b:mx[tag]=err;where[tag]=(j,str(num),str(exact),str(pv))
 sumx+=2*abs(predicted[0]-nums[0])**2
 sumg+=2*lam**2*(abs(predicted[1]-nums[1])**2+2*abs(predicted[2]-nums[2])**2+abs(predicted[3]-nums[3])**2)
 if j%2000==0:print('step',j,'elapsed',round(time.time()-t0,1),flush=True)
cert=mp.iv.sqrt(sumx+sumg)
print('DD token vs exact bath H_N, Frobenius upper',cert.b)
print('max_absolute_errors', {k:str(v.b) for k,v in mx.items()})
T_exact=2000*mp.iv.log(I(20));T_record=I(head[3]);print('T_decimal_vs_formula',abs(T_exact-T_record).b)
print('elapsed_seconds',round(time.time()-t0,2))

print('where',where)
