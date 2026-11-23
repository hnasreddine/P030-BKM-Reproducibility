"""Test C: reconstruction exactly AT an exceptional point of the resonance pair (not only on either side).
Bath J = phi(x) [[2,1],[1,2]]/(2 pi), D = eps*diag(delta/2,-delta/2).  For each eps the detuning delta_*(eps) at which
the two zeros of the cutoff pencil coalesce (double zero on the real axis) is located by Newton iteration."""
import numpy as np, mpmath as mp, json, time, p030lib as L
mp.mp.dps = 60
b = L.bath_ep()
mem = L.Memory(b, 600.0, dtype=np.longdouble)
C0 = np.array([[0.2, 0], [0, 0.2]], complex); bb = np.array([0, 0.1], complex)
rows = []
for lam in [0.2, 0.14, 0.1, 0.07, 0.05]:
    eps = lam*lam; Lc = min(300.0, 2.0/eps); t0 = time.time()
    def kap(zeta, deriv=0):          # K_11(eps*zeta) and derivatives w.r.t. zeta (real zeta), long double -> mp
        v = mem.Khat(np.array([np.longdouble(float(zeta))*np.longdouble(eps)]), Lc, deriv=deriv)[0][0, 0]
        return L.ld_to_mp(np.real(v))*mp.mpf(eps)**deriv
    # f = (zeta+k)^2 + delta^2/4 - k^2/4 ; double zero: f=0, f_zeta=0
    def F(zeta, delta):
        k, k1 = kap(zeta), kap(zeta, 1)
        return ((zeta+k)**2 + delta**2/4 - k**2/4, 2*(zeta+k)*(1+k1) - k*k1/2)
    mp.mp.dps = 30
    zeta_star, delta_star = mp.findroot(lambda z, d: F(z, d), (mp.mpf(-1), mp.mpf(1)), tol=1e-17, maxsteps=60)
    mp.mp.dps = 60
    for tag, delta in (('below', 0.9*float(delta_star)), ('EP', float(delta_star)), ('above', 1.1*float(delta_star))):
        D = eps*np.diag([delta/2, -delta/2])
        cs = [1.0, 2.2, 4.0]; Ts = [c*np.log(1/lam)/eps for c in cs]
        N = int((max(Ts) + 45/eps)/np.pi) + 300
        fr = L.Friedrichs(b, D, lam, N, dtype=np.longdouble); fr.moments(int(fr.R*max(Ts) + 12*(fr.R*max(Ts))**(1/3) + 150))
        S, P, A, Zs = L.contour_Z_hp(mem, D, eps, Lc, -1.0*eps, 1.0*eps, 200, Ts, dps=60)
        ev, ER = mp.eig(A/eps)
        split = abs(ev[0]-ev[1]); condV = mp.cond(ER)
        K0 = np.array([[1, .5], [.5, 1.]])
        for c, T, Z in zip(cs, Ts, Zs):
            tau = float(eps*T)
            V = L.mp_from_cld_matrix(fr.V(T, dps=40)); mp.mp.dps = 60
            E = mp.expm(-mp.mpf(float(tau))*(mp.matrix(K0.tolist()) + 1j*mp.matrix(np.diag([delta/2, -delta/2]).tolist())))
            Dm, Bm, Fm = L.entropy_and_defect(V, C0, bb, dps=60); Dz, Bz, Fz = L.entropy_and_defect(Z, C0, bb, dps=60)
            rec = dict(lam=lam, tag=tag, delta=delta, delta_star=float(delta_star), zeta_star=float(zeta_star), c=c, tau=tau, L=Lc,
                       eig_split=float(split), eigvec_cond=float(condV), err_pole=float(L.mp_norm2(V-Z)), err_davies=float(L.mp_norm2(V-E)),
                       normV=float(L.mp_norm2(V)), S_minus_I_over_eps=float(L.mp_norm2(S-mp.eye(2))/eps),
                       relD=float(abs(Dm/Dz-1)), relF=float(abs(Fm/Fz-1)))
            rows.append(rec)
            print({k: ('%.4g' % v if isinstance(v, float) else v) for k, v in rec.items()}, flush=True)
    json.dump(rows, open('../data/testC_EP.json', 'w'), indent=1)
    print('lam', lam, '%.0fs' % (time.time()-t0), flush=True)
