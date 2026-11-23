from pathlib import Path
HERE=Path(__file__).resolve().parent
import mpmath as mp, time
mp.mp.dps=45
z=mp.mpf(0); a=mp.mpf('0.1'); b=mp.mpf('0.9'); x0=mp.mpf('.5');w=mp.mpf('.4')
def theta(x):
 t=(x-x0)/w
 return mp.exp(1-1/(1-t*t)) if abs(t)<1 else z
bounds=[a,mp.mpf('.16'),mp.mpf('.26'),mp.mpf('.38'),mp.mpf('.5'),mp.mpf('.62'),mp.mpf('.74'),mp.mpf('.84'),b]
def I(k):return mp.quad(lambda x:theta(x)*x**k,bounds,method='gauss-legendre',maxdegree=10)
t=time.time(); Iminus2=I(-2);Iminus1=I(-1);I0=I(0);I1=I(1)
A=mp.matrix([[I0,Iminus1],[Iminus1,Iminus2]]);B=mp.matrix([-I1,-I0]);sol=mp.lu_solve(A,B)
coeff=[sol[1],sol[0],mp.mpf(1)] # a in polynomial = sol[0] ; b=sol[1]
for tag,k in [('I-2',-2),('I-1',-1),('I0',0),('I1',1)]:print(tag,mp.nstr({-2:Iminus2,-1:Iminus1,0:I0,1:I1}[k],40),flush=True)
print('monic p low-to-high:',*[mp.nstr(v,40) for v in coeff],flush=True)
for k in range(3):
 val=sum(coeff[j]*I(j-k-1) for j in range(3))
 print('moment',k,mp.nstr(val,30),flush=True)
print('elapsed',time.time()-t,flush=True)
with open(str(HERE / 'high_precision_monic_coefficients.txt'),'w') as f:
 f.write('\n'.join(mp.nstr(x,44) for x in coeff)+'\n')
