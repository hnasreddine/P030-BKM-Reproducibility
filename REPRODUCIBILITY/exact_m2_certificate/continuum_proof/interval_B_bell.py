from pathlib import Path
HERE=Path(__file__).resolve().parent
import sys, time, mpmath as mp,math
mp.iv.dps=40;mp.mp.dps=45

def bder_interval(lo,hi,k):
 t=mp.iv.mpf([lo,hi]); f=[]
 for r in range(1,k+1):
  v=-(math.factorial(r)/2)*(1/(1-t)**(r+1)+(-1)**r/(1+t)**(r+1))
  f.append(v)
 D=[mp.iv.mpf(1)]
 for n in range(1,k+1):
  D.append(sum(math.comb(n-1,j-1)*f[j-1]*D[n-j] for j in range(1,n+1)))
 return float(abs(D[k]*mp.iv.exp(1-1/(1-t*t))).b)

def tail(k):
 c=mp.mpf('.24');d=mp.mpf('.002')
 return 2*mp.e*mp.factorial(k)*4**k*c**(1-k)*mp.gammainc(k-1,c/d,mp.inf)
N=int(sys.argv[1]) if len(sys.argv)>1 else 4000
start=time.time();a=mp.mpf('-.998');h=mp.mpf('1.996')/N
sums={12:mp.mpf(0),13:mp.mpf(0),14:mp.mpf(0)}
for i in range(N):
 lo=a+h*i;hi=lo+h
 for k in sums:sums[k]+=mp.mpf(bder_interval(mp.nstr(lo,25),mp.nstr(hi,25),k))*h
 if (i+1)%1000==0:print('progress',i+1,'sec',time.time()-start,flush=True)
for k,v in sums.items():
 v+=tail(k)
 print('B_UPPER',k,mp.nstr(v,20),flush=True)
with open(str(HERE / 'interval_B_bell_')+str(N)+'.txt','w') as f:
 for k,v in sums.items():f.write(f'{k} {mp.nstr(v+tail(k),35)}\n')
