"""Validated continuum-to-discrete survival error for P030 m=2, lambda=.05, c=5, N=8500.

Proof uses Volterra resolvent perturbation, Poisson aliasing, exact Sturm root isolation
for bump derivative total variation, and rational Markov polynomial bounds.

This certifies the operator discretization error only, NOT floating-point
Chebyshev propagation or the exact continuum BKM value.
"""
import sympy as sp
import mpmath as mp
mp.iv.dps=65
X=sp.symbols('x'); coeffs=[sp.Poly(1,X,domain=sp.ZZ)]
for k in range(14):
 p=coeffs[-1]
 coeffs.append(p.diff()*sp.Poly((1-X*X)**2,X)+p*sp.Poly(4*k*X*(1-X*X)-2*X,X))

def rat(s):return mp.iv.mpf(str(s.p))/mp.iv.mpf(str(s.q))
def b_der(k,left,right):
 t=mp.iv.mpf([rat(left).a,rat(right).b])
 num=mp.iv.mpf(0)
 for a in coeffs[k].all_coeffs():num=num*t+int(a)
 return mp.iv.exp(1-1/(1-t*t))*num/(1-t*t)**(2*k)

# Exact root isolation, no finite grid root detection.
def variation_bound(k):
 intervals=coeffs[k].intervals(eps=sp.Rational(1,10**20))
 roots=[(l,r) for (l,r),mult in intervals if -1<=l and r<=1]
 f=[mp.iv.mpf(0)]+[b_der(k-1,l,r) for l,r in roots]+[mp.iv.mpf(0)]
 total=mp.iv.mpf(0)
 for a,b in zip(f,f[1:]):total+=abs(b-a)
 return total, len(roots)

B={}
for k,ceiling in [(12,'4.378e21'),(13,'2.412e24'),(14,'1.566e27')]:
 bound,root_count=variation_bound(k)
 C=mp.iv.mpf(ceiling)
 assert bound.b < C.a, (k,bound,C)
 B[k]=C
 print('B'+str(k),'roots',root_count,'validated enclosure',str(bound), 'safe upper',ceiling,flush=True)

# For a degree-two polynomial p normalised by its max over a 20001-point
# [0.1,0.9] grid, Markov => ||p||<=1.001, ||p'||<=10.01, ||p''||<=25.025.
# This holds for both the rounded SVD polynomial used by V0.8 and
# the mathematically exact orthogonal polynomial; no numerical p coefficients needed.
I=mp.iv.mpf;r=14;w=I('0.4');eta=I('1.8')
pmax=(I('1.001'),I('10.01'),I('25.025'))
C=(4/mp.iv.pi +8*w**(1-r))*B[14] + 2*eta*(pmax[0]*w**(1-r)*B[14] +r*pmax[1]*w**(2-r)*B[13] +(r*(r-1)//2)*pmax[2]*w**(3-r)*B[12])
N=8500;lam=I('0.05');T=I(5)*mp.iv.log(I(20))/(lam*lam);P=mp.iv.pi*(N-1)
# Rigorous upper bound zeta(14) <= 1+2^-14 +2^(1-14)/(14-1) by integral comparison.
zeta_upper=1+I(1)/(2**r)+I(1)/(2**(r-1)*(r-1))
E=lam**2*T*T*zeta_upper*C/(P-T)**r
print('C14 upper enclosure:',str(C),flush=True)
print('E_N continuum-to-exact N-model interval:',str(E),flush=True)
assert E.b < I('1e-22').a
print('CERTIFICATE: ||V_continuum(T)-V_N_exact(T)|| < 1e-22',flush=True)
with open(str(__import__('pathlib').Path(__file__).resolve().parent / 'validated_operator_enclosure.txt'),'w') as f:
 f.write('B12: '+str(B[12])+'\nB13: '+str(B[13])+'\nB14: '+str(B[14])+'\n')
 f.write('C14: '+str(C)+'\nE_N: '+str(E)+'\n')
 f.write('CONCLUSION: exact continuum vs exact 8500-mode Hamiltonian < 1e-22.\n')
