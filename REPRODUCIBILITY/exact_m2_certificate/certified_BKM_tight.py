"""Validated interval Riemann upper/lower bounds for the BKM defect of the
EXACT DECIMAL MATRIX in the corrected DD output (not yet V_N exact).

Positive scalar Schur kernel:
 Delta(W) = int_0^inf g( (alpha*u + kappa*detC)/((u^2+s*u+detC)*(a+u)) ) du,
 g(t)=sum_{k>=2}t^k/k, alpha = ||W b*||^2, C=0.2*WW*, s=tr(C),
 detC=0.04*abs(det(W))^2, kappa=b C0^-1 b*=0.05, a=1-s.
Monotone f(u)=g(x(u)), because x decreases for positive u.

Interval arithmetic is used for every quantity/operation; tails use analytic bounds.
"""
import mpmath as mp
from pathlib import Path
import time,sys
mp.iv.dps=55
I=mp.iv.mpf
f=Path(sys.argv[1] if len(sys.argv)>1 else str(Path(__file__).resolve().parent / 'out_ideal_m2_lean.txt'))
rows=[[I(z) for z in ln.split()] for ln in f.read_text().splitlines() if ln.strip()]
assert len(rows)==4 and all(len(r)==4 for r in rows)
v=lambda j: (rows[j][0]+rows[j][1],rows[j][2]+rows[j][3])
V=[[v(0),v(2)],[v(1),v(3)]]
abs2=lambda z:z[0]**2+z[1]**2
mul=lambda z,w:(z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0])
sub=lambda z,w:(z[0]-w[0],z[1]-w[1])
det=sub(mul(V[0][0],V[1][1]),mul(V[0][1],V[1][0]))
s=I('.2')*sum((abs2(V[i][j]) for i in range(2) for j in range(2)),I(0))
alpha=I('.01')*(abs2(V[0][1])+abs2(V[1][1]))
d=I('.04')*abs2(det)
a=I(1)-s
kappa=I('.05')
x=lambda u: (alpha*u+kappa*d)/((u*u+s*u+d)*(a+u))
def g(t):
  # positive series with certified remainder, using max t<0.051
  M=16; power=t*t; val=power/2
  for k in range(3,M+1):
    power=power*t
    val+=power/k
  val += I([0, (t**(M+1)/((M+1)*(1-t))).b])
  return val
lo=lambda z: z.a
hi=lambda z: z.b
u0=I('1e-55')
upper=I('.01')
r=I('1.001')
# Monotonic lower and upper Darboux sums using geometric interval radii.
t0=time.time(); count=0
fleft=g(x(u0))
u=u0; l=I(0);h=I(0)
while u.b < upper.a:
  nxt=u*r
  if nxt.b > upper.a: nxt=upper
  fright=g(x(nxt));du=nxt-u
  l += du*fright
  h += du*fleft
  fleft=fright;u=nxt;count+=1
  if count%20000==0:print('intervals',count,'seconds',round(time.time()-t0,2),flush=True)
  if count>500000: raise RuntimeError('grid too long')
# left tail <=u0*g(x(0)); right tail <=alpha^2/(6*a^2*upper^3)/(1-x(upper))
xmax=kappa/a
lefttail=u0*g(xmax)
# x(u) <= alpha/(a*u*u) + kappa*d/(a*u^3), u >= upper;
# and its ratio to alpha/u²: at u=upper denominator u², numerator <= alpha*u+kappa*d.
# simple x(u) <= (alpha+kappa*d/upper)/(a*u*u) with u >= upper
coef=(alpha+kappa*d/upper)/a
righttail=coef*coef/(I(6)*upper**3*(1-x(upper)))
lo_bound=l.a
hi_bound=(h+lefttail+righttail).b
print('N_intervals',count,'elapsed_seconds',round(time.time()-t0,3))
print('s',s,'alpha',alpha,'detC',d,'a',a)
print('low',lo_bound)
print('high',hi_bound)
print('ratio_high_low',hi_bound/lo_bound)
print('lefttail_up',lefttail.b,'righttail_up',righttail.b)
print('mid_est', (lo_bound+hi_bound)/2)
assert lo_bound>0
assert hi_bound/lo_bound < I('1.0011').a
