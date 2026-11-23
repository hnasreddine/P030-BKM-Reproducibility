from pathlib import Path
HERE=Path(__file__).resolve().parent
"""Rigorous total-variation bounds on L1 norms of derivatives of standard Cinf bump.
SymPy Sturm isolates every real root of exact integer derivative numerator P_k;
mpmath.iv propagates outward interval enclosures through b^(k-1) evaluated at those roots.
"""
import sympy as sy,mpmath as mp,time
mp.iv.dps=55
x=sy.symbols('x')
polys=[sy.Poly(1,x,domain=sy.ZZ)]
for k in range(14):
 p=polys[-1]
 polys.append(p.diff()*sy.Poly((1-x*x)**2,x)+p*sy.Poly(4*k*x*(1-x*x)-2*x,x))

def ivrat(s):return mp.iv.mpf(str(s.p))/mp.iv.mpf(str(s.q))
def eval_bder(k,lo,hi):
 u=mp.iv.mpf([ivrat(lo).a,ivrat(hi).b]); den=1-u*u
 v=mp.iv.mpf(0)
 for co in polys[k].all_coeffs():v=v*u+int(co)
 return v/(den**(2*k))*mp.iv.exp(1-1/den)

for k in (12,13,14):
 intervals=polys[k].intervals(eps=sy.Rational(1,10**20))
 inner=[(lo,hi) for (lo,hi),mult in intervals if hi>=-1 and lo<=1]
 vals=[mp.iv.mpf(0)]
 for lo,hi in inner:vals.append(eval_bder(k-1,lo,hi))
 vals.append(mp.iv.mpf(0))
 total=mp.iv.mpf(0)
 for a,b in zip(vals[:-1],vals[1:]): total+=abs(b-a)
 print('B',k,'roots',len(inner),'interval',total,flush=True)
 with open(str(HERE / 'B')+str(k)+'_cert.txt','w') as f:f.write(str(total)+'\n')
