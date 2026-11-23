"""Static BKM/Schur identities and inequalities: randomized high-precision tests + adversarial cases."""
import mpmath as mp, random, json
mp.mp.dps = 50
random.seed(12345)
def rand_state(d, small=None):
    n = d+1
    A = mp.matrix(n, n)
    for i in range(n):
        for j in range(n): A[i, j] = mp.mpc(random.gauss(0, 1), random.gauss(0, 1))
    R = A*A.H
    if small is not None:           # make the C block tiny / ill conditioned
        Sg = mp.diag([1]+[mp.mpf(10)**(-random.uniform(0, small)) for _ in range(d)])
        R = Sg*R*Sg
    tr = mp.re(sum(R[i, i] for i in range(n)))
    return R/tr
def logm_h(M):
    w, U = mp.eighe(M); return U*mp.diag([mp.log(x) for x in w])*U.H
def tr(M): return mp.re(sum(M[i, i] for i in range(M.rows)))
def Dplus(X, Y): return tr(X*(logm_h(X)-logm_h(Y)) - X + Y)
def g(x):
    # -log(1-x)-x, evaluated without cancellation for tiny x
    if abs(x) < mp.mpf(10)**-12: return x*x/2 + x**3/3 + x**4/4 + x**5/5
    return -mp.log1p(-x)-x
def pquad(f, pts):
    # separate adaptive quadrature on each subinterval (a single mp.quad call over many breakpoints can stop too early)
    return mp.fsum(mp.quad(f, [p, q]) for p, q in zip(pts[:-1], pts[1:]))
worst = dict(det_id=0, def_id=0, dplus_id=0, viol=0, scal=0, inv=0)
ntests = 0
for trial in range(60):
    d = random.choice([1, 2, 3, 4]); small = random.choice([None, 3, 8, 15])
    rho = rand_state(d, small)
    a = mp.re(rho[0, 0]); b = rho[0, 1:]; C = rho[1:, 1:]
    Pi = mp.matrix(d+1, d+1); Pi[0, 0] = a; Pi[1:, 1:] = C
    Dval = tr(rho*(logm_h(rho)-logm_h(Pi)))
    cs, U = mp.eighe(C); beta = b*U
    Bval = mp.fsum(abs(beta[j])**2*(mp.log(a)-mp.log(cs[j]))/(a-cs[j]) for j in range(d))
    x = lambda u: mp.fsum(abs(beta[j])**2/(cs[j]+u) for j in range(d))/(a+u)
    pts = [0]+[mp.mpf(10)**k for k in range(-40, 14, 2)]+[mp.inf]
    I1 = pquad(lambda u: -mp.log(1-x(u)), pts); I2 = pquad(lambda u: g(x(u)), pts)
    worst['det_id'] = max(worst['det_id'], abs(I1/Dval-1)); worst['def_id'] = max(worst['def_id'], abs(I2/(Dval-Bval)-1))
    K = b.H*b/a; kap = max(mp.re(e) for e in mp.eighe(mp.diag([c**-0.5 for c in cs])*U.H*K*U*mp.diag([c**-0.5 for c in cs]), eigvals_only=True))
    assert 0 <= kap < 1
    Dp = Dplus(C-K, C); trK = tr(K); eps_ = max(cs)
    # scalar-kernel version of D_+(C-K||C) = int Tr g(Y_u): rank one => g(y(u)), y(u)=sum |beta_j|^2/(a (c_j+u))
    y = lambda u: mp.fsum(abs(beta[j])**2/(cs[j]+u) for j in range(d))/a
    I3 = pquad(lambda u: g(y(u)), pts); worst['dplus_id'] = max(worst['dplus_id'], abs(I3/Dp-1))
    defect = Dval-Bval
    v1 = defect >= 0; v2 = defect <= Dp*(1+mp.mpf(10)**-30); v3 = Dp <= trK*(1+mp.mpf(10)**-30)
    bound28 = kap/(1-kap)*eps_/a*mp.log(1+a/eps_)*trK
    v4 = (Dp-defect) <= bound28*(1+mp.mpf(10)**-30) and Dp-defect >= -mp.mpf(10)**-40
    ok5 = True
    if eps_ < min(a, mp.e**-1):
        Rr = Dval - (-tr(K*logm_h(C)) + Dp)
        bound29 = trK*((d*eps_+eps_*mp.log(1/eps_))/(a-eps_) + kap/(1-kap)*eps_/a*mp.log(1+a/eps_))
        ok5 = abs(Rr) <= bound29*(1+mp.mpf(10)**-30)
    if not (v1 and v2 and v3 and v4 and ok5): worst['viol'] += 1; print('VIOLATION', trial, v1, v2, v3, v4, ok5)
    # invariants used in the threshold proof: rho_1+rho_2 = b C^{-1} b^dag /1 ; (sum k_i/mu_i)
    inv = mp.fsum(abs(beta[j])**2/cs[j] for j in range(d)); inv2 = (b*mp.inverse(C)*b.H)[0]
    worst['inv'] = max(worst['inv'], abs(inv/mp.re(inv2)-1))
    ntests += 1
# scalar integral bound: int_0^inf u c /((a+u)(c+u)^2) du <= (c/a) log(1+a/c)
for _ in range(400):
    a = mp.mpf(10)**random.uniform(-6, 0); c = mp.mpf(10)**random.uniform(-30, 2)
    I = mp.quad(lambda u: u*c/((a+u)*(c+u)**2), [0, c, a, mp.inf]) if c < a else mp.quad(lambda u: u*c/((a+u)*(c+u)**2), [0, a, c, mp.inf])
    worst['scal'] = max(worst['scal'], I/((c/a)*mp.log(1+a/c)))
out = {k: (float(v) if k != 'viol' else v) for k, v in worst.items()}; out['ntests'] = ntests
print(out)
# Counterexample: an absolute propagator estimate does NOT give relative entropy/defect accuracy without a singular-value floor.
mp.mp.dps = 60
import p030lib as L, numpy as np
C0 = mp.matrix([[0.2, 0], [0, 0.2]]); b = mp.matrix([[0, 0.1]])
ce = []
for delta in ['1e-6', '1e-9', '1e-12']:
    dl = mp.mpf(delta)
    V1 = mp.matrix([[mp.mpf('1e-3'), 0], [0, dl]]); V2 = mp.matrix([[mp.mpf('1e-3'), 0], [0, 2*dl]])
    D1, B1, F1 = L.entropy_and_defect(V1, C0, b, dps=60); D2, B2, F2 = L.entropy_and_defect(V2, C0, b, dps=60)
    ce.append(dict(abs_diff=float(dl), D_ratio=float(D2/D1), defect_ratio=float(F2/F1)))
print('counterexample abs-vs-rel:', ce)
json.dump(dict(random_tests=out, abs_vs_rel_counterexample=ce), open('../data/testS_static.json', 'w'), indent=1)
