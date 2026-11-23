import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np,math
from p030lib_original import bump_derivs
from pathlib import Path
p=np.array([float(x) for x in (Path(__file__).resolve().parent / 'high_precision_monic_coefficients.txt').read_text().splitlines()],dtype=np.longdouble)
p/=np.polyval(p[::-1],np.longdouble('.9'))
x=np.linspace(-1,1,8001,dtype=np.longdouble); w=np.longdouble('.4');eta=np.longdouble('1.8')
# scaled symmetric bumps
D0=bump_derivs(x,20);DP=bump_derivs((x-np.longdouble('.5'))/w,20);DM=bump_derivs((-x-np.longdouble('.5'))/w,20)
px=p[0]+p[1]*x+p[2]*x*x
for r in [2,4,6,8,10,12,14,16,18,20]:
 diag=(1/np.pi)*D0[r]+2*(DP[r]+(-1)**r*DM[r])/w**r
 diag2=(3/np.pi)*D0[r]+2*(DP[r]+(-1)**r*DM[r])/w**r
 off=eta*(DP[r]/w**r*px+r*DP[r-1]/w**(r-1)*(p[1]+2*p[2]*x)+r*(r-1)*DP[r-2]/w**(r-2)*p[2])
 eig=np.sqrt(((diag-diag2)/2)**2+off**2)
 norm=np.maximum(abs((diag+diag2)/2+eig),abs((diag+diag2)/2-eig))
 C=np.trapezoid(norm,x)
 print(r,'numerical_Cr',float(C),flush=True)
