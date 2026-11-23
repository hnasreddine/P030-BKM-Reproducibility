"""Analytic check that the scalar Schur integrand g(x(t)) is decreasing
for the exact-decimal stored matrix. No sampling or derivative approximations.
"""
from pathlib import Path
import mpmath as mp
mp.iv.dps=65;I=mp.iv.mpf
rows=[[I(z) for z in line.split()] for line in (Path(__file__).resolve().parent / 'out_ideal_m2_lean.txt').read_text().splitlines()]
v=lambda i:(rows[i][0]+rows[i][1],rows[i][2]+rows[i][3])
V=[[v(0),v(2)],[v(1),v(3)]]
norm2=lambda z:z[0]*z[0]+z[1]*z[1]
prod=lambda a,b:(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
det=lambda a,b,c,d: (prod(a,d)[0]-prod(b,c)[0],prod(a,d)[1]-prod(b,c)[1])
s=I('.2')*sum((norm2(z) for row in V for z in row),I(0))
a=1-s;d=I('.04')*norm2(det(V[0][0],V[0][1],V[1][0],V[1][1]))
alpha=I('.01')*(norm2(V[0][1])+norm2(V[1][1]));kappa=I('.05')
# numerator of x'(t), with x(t)=(alpha*t+kappa*d)/((t^2+s*t+d)*(a+t))
# polynomial coefficients for t^3, t^2, t, 1:
coefs=[-2*alpha,-alpha*(s+a)-3*kappa*d,-2*kappa*d*(s+a),d*(alpha*a-kappa*(d+s*a))]
for i,z in enumerate(coefs):
 print(f'coefficient_t^{3-i}',z)
 assert z.b<0
x0=kappa/a
assert x0.b<I('.051').a
print('x(0) <0.051 PASS')
print("g(x(t)) nonincreasing t>=0: PASS")
