import json,mpmath as mp
mp.mp.dps=90
raw=json.load(open(__import__('os').path.join(__import__('os').path.dirname(__file__), 'direct_m2_results_reference.json')))
z=next(row for row in raw if row['lam']==0.05 and row['c']==5 and row['N']==8500)
def parse(s):
 s=s.strip('()')
 if ' + ' in s:a,b=s.split(' + ');return mp.mpc(a,b[:-1])
 if ' - ' in s:a,b=s.split(' - ');return mp.mpc(a,'-'+b[:-1])
 raise ValueError(s)
V=mp.matrix([[parse(z['V'][i][j]) for j in range(2)] for i in range(2)])
# singular values of full 2x2 survival matrix, stable via exact determinant
G=V.conjugate().transpose()*V
s1=mp.sqrt(((G[0,0]+G[1,1]).real+mp.sqrt((G[0,0]+G[1,1]).real**2-4*abs(mp.det(V))**2))/2)
s2=abs(mp.det(V))/s1
sC=mp.mpf('.2')*mp.fsum([abs(v)**2 for row in z['V'] for v in []]) if False else mp.mpf('.2')*mp.fsum(abs(V[i,j])**2 for i in range(2) for j in range(2))
aW=1-sC
rho0=mp.mpf('.01')/mp.mpf('.2')
kappa=rho0/aW
E=mp.mpf('9.35499814562769272405542008118e-23')
delta=E/s2
am=(1-delta)**2;ap=(1+delta)**2
amin=1-ap*sC;amax=1-am*sC
cup=max(mp.mpf(1),aW/amin)
clow=min(mp.mpf(1),aW/amax)
Phi=lambda c:c*c*(1-kappa)/(1-c*kappa)
lo=am*Phi(clow);hi=ap*Phi(cup)
print('sigma_min_estimate',mp.nstr(s2,35))
print('delta',mp.nstr(delta,35),'aW',mp.nstr(aW,35),'kappa',mp.nstr(kappa,35))
print('multiplicative_bounds_for_continuum_defect_relative_to_exact_N',mp.nstr(lo,35),mp.nstr(hi,35))
print('relative_lower_percent',mp.nstr((1-lo)*100,20),'upper_percent',mp.nstr((hi-1)*100,20))
print('IMPORTANT conditional: V_N entries from previous non-interval computation, not certified ball')
