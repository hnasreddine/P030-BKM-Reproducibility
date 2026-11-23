from pathlib import Path
HERE=Path(__file__).resolve().parent
import sys,time, mpmath as mp
mp.iv.dps=35
mp.mp.dps=40

def add(a,b):
 o=[0]*max(len(a),len(b))
 for i,x in enumerate(a):o[i]+=x
 for i,x in enumerate(b):o[i]+=x
 return o

def mul(a,b):
 o=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):o[i+j]+=x*y
 return o

def der(a):return [i*v for i,v in enumerate(a)][1:] or [0]
polys=[[1]]
for k in range(14):
 p=polys[-1]
 polys.append(add(mul(der(p),[1,0,-2,0,1]),mul(p,[0,4*k-2,0,-4*k])))

def bder_interval(lo,hi,k):
 t=mp.iv.mpf([lo,hi]); u=1-t*t
 v=mp.iv.mpf(0)
 for co in reversed(polys[k]):v=v*t+co
 v=v/u**(2*k)*mp.iv.exp(1-1/u)
 return float(abs(v).b)

def ktail(k,delta):
 c=mp.mpf('.24');d=mp.mpf(delta)
 return 2*mp.e*mp.factorial(k)*4**k*c**(1-k)*mp.gammainc(k-1,c/d,mp.inf)

N=int(sys.argv[1]) if len(sys.argv)>1 else 2000
start=time.time(); a=mp.mpf('-.998');h=mp.mpf('1.996')/N
sums={12:mp.mpf(0),13:mp.mpf(0),14:mp.mpf(0)}
for i in range(N):
 lo=a+h*i;hi=lo+h
 for k in sums:
  sums[k]+=mp.mpf(bder_interval(mp.nstr(lo,25),mp.nstr(hi,25),k))*h
 if (i+1)%1000==0:print('progress',i+1,'sec',time.time()-start,flush=True)
for k,v in sums.items():
 v+=ktail(k,mp.mpf('.002'))
 print('B_CERTIFIED_UPPER',k,mp.nstr(v,20),'time',time.time()-start,flush=True)
with open(str(HERE / 'interval_B_bounds_')+str(N)+'.txt','w') as f:
 for k,v in sums.items():f.write(f'k={k}: {mp.nstr(v+ktail(k,mp.mpf(".002")),30)}\n')
