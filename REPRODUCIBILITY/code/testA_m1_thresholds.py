"""Test A: fixed physical bath, m=1 (V0.3 Example). Direct microscopic propagation (long double Chebyshev)
vs (i) pole reconstruction exp(TA)S (long double/mp), (ii) bare Davies, (iii) two-term asymptotic formulas."""
import numpy as np, mpmath as mp, json, sys, time, p030lib as L
mp.mp.dps = 80
lams = [float(a) for a in sys.argv[1].split(',')] if len(sys.argv) > 1 else [0.2, 0.14, 0.1, 0.07, 0.05]
out = sys.argv[2] if len(sys.argv) > 2 else '../data/testA_m1.json'
gs, gf, zeta = 1.0, 2.0, 0.4
b = L.bath_m1(gs, gf, zeta); D = np.zeros((2, 2))
kappa = b.kappa[1]; dg = gf-gs; m = 1
h = -(-gs)**m*kappa/dg
cs = [0.5+0.25*i for i in range(25)]            # c in [0.5, 6.5]
states = {
 'iso':  dict(C0=np.array([[0.2, 0], [0, 0.2]], complex), b=np.array([0, 0.1], complex)),
 'gen':  dict(C0=np.array([[0.25, 0.08+0.03j], [0.08-0.03j, 0.15]], complex), b=np.array([0, 0.12], complex)),
}
for st in states.values():
    C0, bb = st['C0'], st['b']; a0 = 1-np.trace(C0).real
    rho = np.block([[np.array([[a0]]), bb[None, :]], [bb.conj()[:, None], C0]])
    assert np.linalg.eigvalsh(rho).min() > 0
    Ci = np.linalg.inv(C0)
    st['css'] = C0[0, 0].real; st['cf'] = 1/Ci[1, 1].real; st['bf2'] = abs(bb[1])**2
    st['d0'] = float(L.dplus_scalar(st['cf'], st['bf2']))
mem = L.Memory(b, 600.0, dtype=np.longdouble)
rows = []
for lam in lams:
    eps = lam*lam; t0 = time.time()
    Tmax = max(cs)*np.log(1/lam)/eps
    N = int((Tmax + 45/eps)/np.pi) + 300
    fr = L.Friedrichs(b, D, lam, N, dtype=np.longdouble)
    fr.moments(int(fr.R*Tmax + 12*(fr.R*Tmax)**(1/3) + 150))
    Lc = min(300.0, 2.0/eps)
    Ts = [c*np.log(1/lam)/eps for c in cs]
    S, P, A, Zs = L.contour_Z_hp(mem, D, eps, Lc, -1.5*eps, 1.2*eps, 160, Ts, dps=80)
    ev = mp.eig(A, left=False, right=False); ev = sorted([e/eps for e in ev], key=lambda e: -mp.re(e)); zs, zf = ev
    print('lam', lam, 'zeta_s', mp.nstr(zs, 12), 'zeta_f', mp.nstr(zf, 12), 'N', N, 'K', fr.K, 'L', Lc, 'setup %.1fs' % (time.time()-t0), flush=True)
    for c, T, Z in zip(cs, Ts, Zs):
        tau = eps*T
        V = L.mp_from_cld_matrix(fr.V(T, dps=40)); mp.mp.dps = 80
        E = mp.matrix([[mp.e**(-gs*tau), 0], [0, mp.e**(-gf*tau)]])
        rec = dict(zeta_s=[float(mp.re(zs)), float(mp.im(zs))], zeta_f=[float(mp.re(zf)), float(mp.im(zf))], lam=lam, c=c, tau=tau, T=T, N=N, K=fr.K, L=Lc,
                   err_pole=float(L.mp_norm2(V-Z)), err_davies=float(L.mp_norm2(V-E)),
                   V_sf_abs=float(abs(V[0, 1])), V_ff_abs=float(abs(V[1, 1])),
                   V_sf_pred=float(abs(eps**m*kappa*((-gf)**m*mp.e**(-gf*tau)-(-gs)**m*mp.e**(-gs*tau))/dg)))
        for name, st in states.items():
            Dm, Bm, Fm = L.entropy_and_defect(V, st['C0'], st['b'])
            Dz, Bz, Fz = L.entropy_and_defect(Z, st['C0'], st['b'])
            De, Be, Fe = L.entropy_and_defect(E, st['C0'], st['b'])
            bf2, css, cf, d0 = st['bf2'], st['css'], st['cf'], st['d0']
            ef, es = mp.e**(-2*gf*tau), mp.e**(-2*gs*tau)
            h2 = abs(h)**2
            D_lead = 2*tau*bf2*(gf*ef + gs*h2*eps**(2*m)*es)
            D_ref = bf2*ef*(2*gf*tau - mp.log(cf) + d0/bf2) + bf2*h2*eps**(2*m)*es*(2*gs*tau - mp.log(css))
            F_asym = d0*ef + bf2**2*h2**2*eps**(4*m)*es/(2*css)
            rf, rs = mp.e**(2*mp.re(zf)*tau), mp.e**(2*mp.re(zs)*tau)      # renormalized resonance decay factors
            D_ren = bf2*rf*(-mp.log(rf) - mp.log(cf) + d0/bf2) + bf2*h2*eps**(2*m)*rs*(-mp.log(rs) - mp.log(css))
            F_ren = d0*rf + bf2**2*h2**2*eps**(4*m)*rs/(2*css)
            rec[name] = dict(D_micro=float(Dm), F_micro=float(Fm), D_pole=float(Dz), F_pole=float(Fz),
                             D_davies=float(De), F_davies=float(Fe),
                             relD_pole=float(abs(Dm/Dz-1)), relF_pole=float(abs(Fm/Fz-1)),
                             D_lead=float(D_lead), D_ref=float(D_ref), F_asym=float(F_asym), D_ren=float(D_ren), F_ren=float(F_ren),
                             D_fast=float(2*tau*bf2*gf*ef), D_slow=float(2*tau*bf2*gs*h2*eps**(2*m)*es),
                             F_fast=float(d0*ef), F_slow=float(bf2**2*h2**2*eps**(4*m)*es/(2*css)))
        rows.append(rec)
        json.dump(rows, open(out, 'w'), indent=1)
    print('  done lam', lam, '%.1fs' % (time.time()-t0), flush=True)
