"""Rigorous (conditional on imported certified B12..B14 bounds) trapezoid+interval
certificate of exact order-2 orthogonal polynomial coefficients.
With composite trapezoid over [.1,.9] for compactly supported flat bump,
Poisson remainder <= 2*zeta(14)*||f^(14)||_1*(h/(2*pi))^14.
All numerical quadrature and algebra is outward mpmath.iv at 72 digits.
"""
from pathlib import Path
import mpmath as m,sys,time,math
m.iv.dps=72
I=m.iv.mpf
N=int(sys.argv[1]) if len(sys.argv)>1 else 150000
z=I(0); w=I('.4');h=I('.8')/N
# Imported validated upper bounds of B_k=||b^(k)||_1 (documented in CERTIFICATE...)
B={k:I(2)**(12-k)*I('4.378e21') for k in range(13)}
B[13]=I('2.412e24');B[14]=I('1.566e27')
r=14

def deriv_sup(j,k):
 if j==-2:return I(math.factorial(k+1)*10**(k+2))
 if j==-1:return I(math.factorial(k)*10**(k+1))
 if j==0:return I(1 if k==0 else 0)
 if j==1:return I('.9') if k==0 else (I(1) if k==1 else I(0))
 raise ValueError(j)

def C14(j):
 return sum((I(math.comb(14,l))*w**(1-l)*B[l]*deriv_sup(j,14-l) for l in range(15)),I(0))
# zeta(14)<1+2^-14+2^-13/13 <1.001; pi > 3.14.
errors={j:I('2.002')*C14(j)*(h/I('6.28'))**14 for j in [-2,-1,0,1]}
print('N',N,'h',h,'Conservative quadrature bounds (upper)',{j:str(errors[j].b) for j in errors},flush=True)
t0=time.time();S={j:I(0) for j in errors}
for k in range(1,N):
 t=I(-1)+I(2*k)/N
 y=I(1)-t*t
 v=m.iv.exp(I(1)-I(1)/y)
 x=I('.1')+I('.8')*k/N
 inv=I(1)/x
 S[-2] += v*inv*inv
 S[-1] += v*inv
 S[0] += v
 S[1] += v*x
 if k%25000==0:print('integrated',k,'seconds',round(time.time()-t0,1),flush=True)
intv={j:S[j]*h+I([-errors[j].b,errors[j].b]) for j in errors}
for j in [-2,-1,0,1]:print('I',j,intv[j],flush=True)
# p(x)=x²+a*x+b, with b I(-2)+a I(-1)+I0=0, b I(-1)+a I0+I1=0
I2,I1,I0,Ip=intv[-2],intv[-1],intv[0],intv[1]
det=I2*I0-I1*I1
assert det.a>z.a
b=(I1*Ip-I0*I0)/det
a=(I1*I0-I2*Ip)/det
P=lambda x:x*x+a*x+b
vp=P(I('.1'));vq=P(I('.9'));vt=b-a*a/4
assert vq.a>abs(vp).b and vq.a>abs(vt).b
M=vq
coef=[b/M,a/M,I(1)/M]
written=[I(q) for q in Path(str(Path(__file__).resolve().parent / 'orthogonal_coefs.txt')).read_text().splitlines()]
diff=[abs(v-u) for v,u in zip(coef,written)]
diffup=sum((q.b for q in diff),z)
print('det_lower',det.a,flush=True)
print('b, a intervals',b,a,flush=True)
print('max proof: p(.9) > max abs p(.1), abs p(vertex)',vq.a,abs(vp).b,abs(vt).b,flush=True)
for j,q in enumerate(coef):print('normalized_coef',j,q,'delta upper',diff[j].b,flush=True)
print('TOTAL_SUM_COEFFICIENT_ERROR_UPPER',diffup,flush=True)
assert diffup < I('1e-27').a, 'certified orthogonal polynomial difference exceeds threshold'
print('CERTIFIED: for all |x| <=1, |p_exact(x)-p_written(x)|<1e-27',flush=True)
print('elapsed_seconds',round(time.time()-t0,2))
