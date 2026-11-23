from pathlib import Path
HERE=Path(__file__).resolve().parent
import mpmath as mp
mp.mp.dps=65
# Based on b(t)=exp(1-1/(1-t^2)), |b^(r)(t)| <= e r! (4/d)^r exp(-0.24/d) , d=1-|t|
# rigorous because Re 1/(1-z²) >= 0.24/d on |z-t|=d/4 for t real |t|<1

def K(r):
 if r==0:return mp.mpf(2)
 if r==1:return mp.mpf(2)  # total variation of bump, monotone each side
 c=mp.mpf('0.24')
 return 2*mp.e*mp.factorial(r)*4**r*c**(1-r)*mp.gamma(r-1)

coef=[mp.mpf(s.strip()) for s in open(str(HERE / 'high_precision_monic_coefficients.txt'))]
# normalise to sup exactly at right support endpoint x=0.9
p_at_09=sum(a*mp.mpf('.9')**i for i,a in enumerate(coef))
coef=[a/p_at_09 for a in coef]
def pmaxder(j):
 lo,hi=mp.mpf('.1'),mp.mpf('.9')
 if j==0:return mp.mpf(1)  # certified from extrema endpoints+vertex
 if j==1:return max(abs(coef[1]+2*coef[2]*lo),abs(coef[1]+2*coef[2]*hi))
 if j==2:return 2*abs(coef[2])
 return mp.mpf(0)

w=mp.mpf('.4');eta=mp.mpf('1.8');lam=mp.mpf('.05');N=8500;c=mp.mpf(5)
h=mp.mpf(2)/(N-1);P=2*mp.pi/h;T=c*mp.log(1/lam)/(lam**2)
print('normalised p coefficients',*[mp.nstr(a,25) for a in coef]);print('p deriv max',[mp.nstr(pmaxder(j),18) for j in range(3)])
print('P',mp.nstr(P,20),'T',mp.nstr(T,20),'margin',mp.nstr(P-T,20))
best=(mp.inf,None,None)
for r in range(2,41):
 # spectral matrix op norm <= |j00|+|j11|+2|j01|; derivative integrated
 # diag base: (1+3)/pi K(r), off-res: each diag 2*theta(abs x) -> sum diag 8*w^(1-r)*K(r)
 C=(4/mp.pi)*K(r)+8*w**(1-r)*K(r)
 off=sum(mp.binomial(r,j)*pmaxder(j)*w**(1-(r-j))*K(r-j) for j in range(3) if j<=r)
 C+=2*eta*off
 E=lam**2*T**2*mp.zeta(r)*C/(P-T)**r
 if E<best[0]:best=(E,r,C)
 if r in [2,4,6,8,10,12,14,16,18,20,24,28,32,36,40]:print('r',r,'C_analytic',mp.nstr(C,5),'E_bound',mp.nstr(E,7))
print('BEST',mp.nstr(best[0],12),'r',best[1],'C',mp.nstr(best[2],10))
with open(str(HERE / 'analytic_bounds.txt'),'w') as f:f.write('Best rigorous analytic continuum error bound: E <= '+mp.nstr(best[0],30)+' using r='+str(best[1])+'\n')
