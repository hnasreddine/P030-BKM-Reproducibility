from pathlib import Path
"""Validated Bessel coefficient enclosure via first-order ODE residual.
For F(t)=e^{iu sin t}=sum_k J_k(u)e^{ikt}, define approximate
Fhat(t)=sum_{k=-K}^K jhat_k e^{ikt}, jhat_-k=(-1)^k jhat_k.
Because F'=iu cos(t) F, a variation-of-constants estimate yields
 sup_{|t|<=pi}|Fhat(t)-F(t)| <= |Fhat(0)-1|+2*pi*sum_k |r_k|,
 r_k=k*jhat_k-u/2*(jhat_{k-1}+jhat_{k+1}).
Every Fourier coefficient has this error. Combined Chebyshev series error
<= (2K+1)*error, provided ||A||<=1. Note that u token in input is
used; time-string truncation must be bounded separately.
"""
import mpmath as mp, pathlib
I=mp.iv.mpf
mp.iv.dps=65
path=pathlib.Path(str(Path(__file__).resolve().parent / 'input_005_5_8500_ideal_m2.txt'))
with path.open() as f:
 header=f.readline().split();N,K=int(header[0]),int(header[1]); R=I(header[2]); T=I(header[3]);u=R*T
 for _ in range(N):f.readline()
 J=[I(f.readline().strip()) for k in range(K+1)]
 j0=J[0]
 normal=j0 + 2*sum((J[k] for k in range(2,K+1,2)),I(0))
 res=I(0)
 rmax=I(0)
 for k in range(1,K+1):
  term=k*J[k] - u*(J[k-1]+(J[k+1] if k<K else I(0)))/2
  res+=abs(term)
  if abs(term).b>rmax.b:rmax=abs(term)
 res+=abs(u*J[K]/2) # k=K+1 endpoint
 bound=abs(normal-1)+4*mp.iv.pi*res  # sum +/- k, t length 2pi
 # safer: phase t from 0 to +-pi; 2*pi overkill
 coeffbound=(2*K+1)*bound
 print('N',N,'K',K)
 print('u',u)
 print('normalization_interval',normal)
 print('sum_positive_k_r',res)
 print('max_r',rmax)
 print('sup_F_error_upper',bound.b)
 print('Chebyshev_Bessel_total_error_bound',coeffbound.b)
