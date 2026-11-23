import mpmath as mp, time
from pathlib import Path
mp.mp.dps=63
z=mp.mpf(0); l=mp.mpf('.1');r=mp.mpf('.9')
def theta(x):
 t=(x-mp.mpf('.5'))/mp.mpf('.4');return mp.exp(1-1/(1-t*t)) if abs(t)<1 else z
br=[l,mp.mpf('.16'),mp.mpf('.26'),mp.mpf('.38'),mp.mpf('.5'),mp.mpf('.62'),mp.mpf('.74'),mp.mpf('.84'),r]
t=time.time()
I={}
for k in [-2,-1,0,1]:
 I[k]=mp.quad(lambda x: theta(x)*x**k,br,method='gauss-legendre',maxdegree=10)
 print('moment',k,mp.nstr(I[k],63),'elapsed',round(time.time()-t,1),flush=True)
# p2(x)=x²+a x+b, ∫θ p/x²=0 and ∫θ p/x=0
A=mp.matrix([[I[-2],I[-1]],[I[-1],I[0]]]); B=mp.matrix([-I[0],-I[1]])
b,a=mp.lu_solve(A,B)
p2=lambda x:x*x+a*x+b
maxx=max([mp.mpf('.1'),mp.mpf('.9'),-a/2],key=lambda x:abs(p2(x)))
M=abs(p2(maxx)); p=[b/M,a/M,1/M]
print('monic',*(mp.nstr(y,55) for y in [b,a,mp.mpf(1)]));print('norm',mp.nstr(M,55), 'maxx',mp.nstr(maxx,35))
print('pcoefs',*(mp.nstr(y,55) for y in p));print('residuals',*(mp.nstr(mp.fsum([p[0]*I[-2+j],p[1]*I[-1+j],p[2]*I[j]]),30) for j in range(2)))
with open(str(Path(__file__).resolve().parent / 'orthogonal_coefs.txt'),'w') as f:
 for y in p:f.write(mp.nstr(y,60)+'\n')
print('TOTAL SEC',round(time.time()-t,1))
