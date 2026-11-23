"""Test B: order m=2.
 Part 1: reproduce the V0.3 'failed' m=2 test with its own bath (three narrow bumps, h=0.075, eta=0.025),
         extend to smaller coupling, and identify the discrepancy via the first-Born memory tail.
 Part 2: polynomially modulated off-resonant band (bath_poly): thresholds c_D=2m/dg, c_BKM=4m/dg."""
import numpy as np, mpmath as mp, json, sys, time, p030lib as L
from numpy.polynomial.legendre import leggauss
mp.mp.dps = 80
part = sys.argv[1]

def v03_bath():
    x, wx = leggauss(2600); centers = np.array([.22, .50, .78]); h = .075
    B = np.vstack([L.bump((x-c)/h) for c in centers]).T
    mom = np.vstack([wx@(B/x[:, None]**(k+1)) for k in range(4)])
    _, _, v = np.linalg.svd(mom[:2, :]); coeff = v[-1, :]; q = B@coeff; coeff = coeff/np.max(abs(q))
    eta = .025
    mu = mom@coeff
    def qf(xx):
        return sum(xx.dtype.type(coeff[i])*L.bump((xx-xx.dtype.type(centers[i]))/xx.dtype.type(h)) for i in range(3))
    def Jf(xx):
        ph = L.bump(xx); J = np.zeros(xx.shape+(2, 2), dtype=ph.dtype)
        J[..., 0, 0] = ph/np.pi; J[..., 1, 1] = 2*ph/np.pi
        J[..., 0, 1] = J[..., 1, 0] = eta*qf(xx)
        return J
    b = L.Bath(Jf, 2, 'v03_m2'); b.mu = mu; b.eta = eta; b.qf = qf
    b.kappa = {k: eta*(-1)**k/(1j)**(k+1)*mu[k] for k in range(4)}
    return b

def born_tail(bath_q, eta, eps, T, Srange=6000.0, ds=0.25, nx=40001):
    """-eps * int_T^inf (s-T) F_sf(s) ds, F_sf(s)=eta*int q(x) e^{-ixs} dx (trapezoid in x, Simpson in s)."""
    x = np.linspace(-1, 1, nx); dx = x[1]-x[0]; q = bath_q(x)
    s = T + np.arange(0, Srange+ds/2, ds)
    if len(s) % 2 == 0: s = s[:-1]
    F = np.empty(len(s), complex)
    for i in range(0, len(s), 500):
        F[i:i+500] = eta*dx*(np.exp(-1j*np.outer(s[i:i+500], x))@q)
    w = np.ones(len(s)); w[1:-1:2] = 4; w[2:-1:2] = 2; w *= ds/3
    return -eps*np.sum(w*(s-T)*F), abs(F[0])

if part == '1':
    b = v03_bath(); D = np.zeros((2, 2)); gs, gf, dg, m = 1.0, 2.0, 1.0, 2
    kap = b.kappa[2]
    print('v0.3 bath: mu', b.mu, 'kappa_2', kap, 'kappa_3', b.kappa[3], 'min eig', b.min_eig(20001), flush=True)
    rows = []
    for lam in [0.2, 0.15, 0.12, 0.10, 0.08, 0.065, 0.05, 0.04, 0.032, 0.026]:
        eps = lam*lam; tau = 2*np.log(1/lam); T = tau/eps; t0 = time.time()
        N = int((T + 45/eps)/np.pi) + 300
        fr = L.Friedrichs(b, D, lam, N, dtype=np.longdouble)
        fr.moments(int(fr.R*T + 12*(fr.R*T)**(1/3) + 150))
        V = fr.V(T, dps=40)
        Vsf = complex(V[0, 1])
        pole = eps**m*kap*((-gf)**m*np.exp(-gf*tau) - (-gs)**m*np.exp(-gs*tau))/dg
        v03pred = lam**4*(-kap)*(np.exp(-tau)-np.exp(-2*tau))
        bt, Fabs = born_tail(lambda x: b.qf(x), b.eta, eps, T)
        rec = dict(lam=lam, tau=tau, T=T, N=N, K=fr.K, Vsf=[Vsf.real, Vsf.imag], Vsf_abs=abs(Vsf),
                   pole_abs=abs(pole), ratio_v03=abs(Vsf)/abs(v03pred), ratio_pole=abs(Vsf)/abs(pole),
                   born_tail=[bt.real, bt.imag], born_tail_abs=abs(bt), F_sf_at_T=Fabs,
                   resid_abs=abs(Vsf-pole), resid_over_born=abs(Vsf-pole)/abs(bt),
                   resid_after_born=abs(Vsf-pole-bt)/abs(pole))
        rows.append(rec); print({k: (('%.4g' % v) if isinstance(v, float) else v) for k, v in rec.items() if k not in ('Vsf', 'born_tail')}, '%.0fs' % (time.time()-t0), flush=True)
        json.dump(rows, open('../data/testB_v03_m2.json', 'w'), indent=1)

if part == '2':
    gs, gf = 1.0, 3.0; m = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    lams = [float(a) for a in sys.argv[3].split(',')] if len(sys.argv) > 3 else [0.2, 0.14, 0.1, 0.07, 0.05]
    b = L.bath_poly(m, gs=gs, gf=gf); D = np.zeros((2, 2)); dg = gf-gs
    kap = b.kappa[m]; h = -(-gs)**m*kap/dg
    print('poly bath m=%d: eta %.5g A %.3g p %s mu %s kappa_m %s h_m %s min eig %.3g' % (m, b.eta, b.A, b.pcoef, b.mu, kap, h, b.min_eig(8001)), flush=True)
    C0 = np.array([[0.2, 0], [0, 0.2]], complex); bb = np.array([0, 0.1], complex)
    css, cf, bf2 = 0.2, 0.2, 0.01; d0 = float(L.dplus_scalar(cf, bf2))
    cD, cB = 2*m/dg, 4*m/dg
    cs = sorted(set([round(0.25*i, 2) for i in range(2, int(4*(cB+1.5))+1)]))
    mem = L.Memory(b, 1100.0, dtype=np.longdouble)
    rows = []
    for lam in lams:
        eps = lam*lam; t0 = time.time()
        Tmax = max(cs)*np.log(1/lam)/eps
        N = int((Tmax + 45/eps)/np.pi) + 300
        fr = L.Friedrichs(b, D, lam, N, dtype=np.longdouble)
        fr.moments(int(fr.R*Tmax + 12*(fr.R*Tmax)**(1/3) + 150))
        Lc = min(550.0, 5.0/(eps*(gf+0.8)))
        Ts = [c*np.log(1/lam)/eps for c in cs]
        S, P, A, Zs = L.contour_Z_hp(mem, D, eps, Lc, -(gs+gf)/2*eps, (dg/2+0.8)*eps, 200, Ts, dps=80)
        ev = sorted([e/eps for e in mp.eig(A, left=False, right=False)], key=lambda e: -mp.re(e)); zs, zf = ev
        print('lam', lam, 'zeta_s', mp.nstr(zs, 10), 'zeta_f', mp.nstr(zf, 10), 'N', N, 'K', fr.K, 'L', Lc, flush=True)
        for c, T, Z in zip(cs, Ts, Zs):
            tau = eps*T
            V = L.mp_from_cld_matrix(fr.V(T, dps=40)); mp.mp.dps = 80
            E = mp.matrix([[mp.e**(-gs*tau), 0], [0, mp.e**(-gf*tau)]])
            Dm, Bm, Fm = L.entropy_and_defect(V, C0, bb); Dz, Bz, Fz = L.entropy_and_defect(Z, C0, bb); De, Be, Fe = L.entropy_and_defect(E, C0, bb)
            h2 = abs(h)**2
            ef, es = mp.e**(-2*gf*tau), mp.e**(-2*gs*tau)
            rf, rs = mp.e**(2*mp.re(zf)*tau), mp.e**(2*mp.re(zs)*tau)
            D_lead = 2*tau*bf2*(gf*ef + gs*h2*eps**(2*m)*es)
            D_ren = bf2*rf*(-mp.log(rf) - mp.log(cf) + d0/bf2) + bf2*h2*eps**(2*m)*rs*(-mp.log(rs) - mp.log(css))
            F_asym = d0*ef + bf2**2*h2**2*eps**(4*m)*es/(2*css)
            F_ren = d0*rf + bf2**2*h2**2*eps**(4*m)*rs/(2*css)
            pole_sf = eps**m*kap*((-gf)**m*mp.e**(-gf*tau)-(-gs)**m*mp.e**(-gs*tau))/dg
            rows.append(dict(m=m, lam=lam, c=c, tau=tau, T=T, N=N, K=fr.K, L=Lc, cD=cD, cB=cB,
                 err_pole=float(L.mp_norm2(V-Z)), err_davies=float(L.mp_norm2(V-E)),
                 V_sf_abs=float(abs(V[0, 1])), V_ff_abs=float(abs(V[1, 1])), V_sf_pred=float(abs(pole_sf)),
                 D_micro=float(Dm), F_micro=float(Fm), D_pole=float(Dz), F_pole=float(Fz), D_davies=float(De), F_davies=float(Fe),
                 relD_pole=float(abs(Dm/Dz-1)), relF_pole=float(abs(Fm/Fz-1)),
                 D_lead=float(D_lead), D_ren=float(D_ren), F_asym=float(F_asym), F_ren=float(F_ren),
                 D_fast=float(2*tau*bf2*gf*ef), D_slow=float(2*tau*bf2*gs*h2*eps**(2*m)*es),
                 F_fast=float(d0*ef), F_slow=float(bf2**2*h2**2*eps**(4*m)*es/(2*css))))
            json.dump(rows, open('../data/testB_poly_m%d.json' % m, 'w'), indent=1)
        print('  done', lam, '%.0fs' % (time.time()-t0), flush=True)
