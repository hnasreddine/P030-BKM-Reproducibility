from pathlib import Path
HERE=Path(__file__).resolve().parent
import sympy as s,time,sys
x=s.symbols('x');P=s.Poly(1,x,domain=s.ZZ)
for k in range(14):
 P=P.diff()*s.Poly((1-x*x)**2,x)+P*s.Poly(4*k*x*(1-x*x)-2*x,x)
 if k+1 in (12,13,14):
  print('k',k+1,'degree',P.degree(),'height',max(len(str(abs(int(v)))) for v in P.all_coeffs()),flush=True)
  start=time.time()
  ints=P.intervals(eps=s.Rational(1,10**13))
  ints=[(l,r,m) for ((l,r),m) in ints if r>=-1 and l<=1]
  print('roots',len(ints),'elapsed',time.time()-start,'sample',ints[:2],flush=True)
  with open(str(HERE / 'roots')+str(k+1)+'.txt','w') as f:
   for l,r,m in ints:f.write(str(l)+' '+str(r)+' '+str(m)+'\n')
