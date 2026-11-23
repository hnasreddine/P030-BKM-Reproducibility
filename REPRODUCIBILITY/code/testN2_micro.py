"""Test N2: genuinely nonnormal Davies generator realized by a microscopic bath.
 J(x) = Gamma*phi(x)/(2 pi), Gamma=[[2,-q],[-q,4]], D = eps*H, H = i(M0-M0^dag)/2, M0=[[-1,q],[0,-2]], q=1/2.
 w = b0(-q,1)^T is annihilated by the slow spectral projector P_s = u l_s (not an orthogonal projection).
Compares: direct microscopic propagation (80-bit), grouped-resonance reconstruction, Davies semigroup,
the nonnormal two-term formulas, and the diagonal formula transplanted by similarity (which is wrong)."""
import numpy as np, mpmath as mp, json, sys, time, p030lib as L, nonnormal as N
mp.mp.dps = 80
lams = [float(a) for a in sys.argv[1].split(',')] if len(sys.argv) > 1 else [0.2, 0.14, 0.1, 0.07, 0.05, 0.035]
q = 0.5; b0 = 0.1
M0 = np.array([[-1, q], [0, -2.]], complex)
Gam = -(M0+M0.conj().T).real; H = 1j*(M0-M0.conj().T)/2
assert np.allclose(M0, -Gam/2-1j*H) and np.linalg.eigvalsh(Gam).min() > 0
assert np.linalg.norm(Gam@H-H@Gam) > 0.1                       # Gamma and H do not commute
def Jf(x):
    ph = L.bump(x); return ph[..., None, None]*(Gam/(2*np.pi)).astype(ph.dtype)
bath = L.Bath(Jf, 2, 'nonnormal')
zs0, zf0, u, ls, vf, lf = N.davies_data(M0)
w = b0*np.array([-q, 1.0], complex); assert abs(ls@w) < 1e-15
mem = L.Memory(bath, 600.0, dtype=np.longdouble)
G1 = np.asarray(mem.jets(300.0, 1)[1], complex)               # first jet (cutoff L=300; equals the physical jet to ~1e-10)
k1 = (G1[0, 0]/Gam[0, 0]).real; b1 = -k1
alpha = -(zs0/(zs0-zf0))*(ls@G1@w)                            # general m=1 coefficient
alpha_reference = -b1*b0*(q+q**3)
print('k1 = K_scalar\'(0) = %.8f, b1 = %.8f, alpha = %s, reference -b1 b0 (q+q^3) = %.8f' % (k1, b1, alpha, alpha_reference), flush=True)
states = {'pI': 0.2*np.eye(2, dtype=complex), 'general': np.array([[0.25, 0.08+0.03j], [0.08-0.03j, 0.15]])}
Ks = {k: N.constants(u, ls, w, C0) for k, C0 in states.items()}
for k, C0 in states.items():
    rho = np.block([[np.array([[1-np.trace(C0).real]]), w.conj()[None, :]], [w[:, None], C0]]); assert np.linalg.eigvalsh(rho).min() > 0
cs = [0.5+0.25*i for i in range(25)]
rows = []; meta = dict(q=q, b0=b0, k1=k1, b1=b1, alpha=[alpha.real, alpha.imag], alpha_proposed=alpha_reference, u=u.real.tolist(), l_s=ls.real.tolist(),
                       consts={k: {a: (float(abs(v)) if a == 'wpar' else float(v)) for a, v in K.items()} for k, K in Ks.items()})
for lam in lams:
    eps = lam*lam; t0 = time.time(); D = eps*H
    Tmax = max(cs)*np.log(1/lam)/eps
    Nn = int((Tmax + 45/eps)/np.pi) + 300
    fr = L.Friedrichs(bath, D, lam, Nn, dtype=np.longdouble); fr.moments(int(fr.R*Tmax + 12*(fr.R*Tmax)**(1/3) + 150))
    Lc = min(300.0, 2.0/eps); Ts = [c*np.log(1/lam)/eps for c in cs]
    S, P, A, Zs = L.contour_Z_hp(mem, D, eps, Lc, -1.5*eps, 1.2*eps, 160, Ts, dps=80)
    ev, ER = mp.eig(A/eps); idx = sorted(range(2), key=lambda i: -mp.re(ev[i])); zs, zf = ev[idx[0]], ev[idx[1]]
    # numerical slow residue amplitude: Q_s w = a_s * u_eps ; Q_s from the spectral decomposition of A and S
    ERi = mp.inverse(ER); us = ER[:, idx[0]]; Qs = us*ERi[idx[0], :]*S
    wm = mp.matrix([[mp.mpc(complex(x))] for x in w]); qsw = Qs*wm
    a_s = mp.sqrt(abs(qsw[0])**2+abs(qsw[1])**2)
    print('lam', lam, 'zeta_s', mp.nstr(zs, 10), 'zeta_f', mp.nstr(zf, 10), '|Q_s w|/(eps|alpha|) =', mp.nstr(a_s/(eps*abs(alpha)), 8), 'S-I/eps', mp.nstr(L.mp_norm2(S-mp.eye(2))/eps, 5), flush=True)
    for c, T, Z in zip(cs, Ts, Zs):
        tau = float(eps*T)
        V = L.mp_from_cld_matrix(fr.V(T, dps=40)); mp.mp.dps = 80
        E = mp.expm(mp.mpf(tau)*mp.matrix(M0.tolist()))
        rec = dict(lam=lam, c=c, tau=tau, T=T, N=Nn, K=fr.K, L=Lc, zeta_s=float(mp.re(zs)), zeta_f=float(mp.re(zf)), Qs_w_over_eps_alpha=float(a_s/(eps*abs(alpha))),
                   err_pole=float(L.mp_norm2(V-Z)), err_davies=float(L.mp_norm2(V-E)), normV=float(L.mp_norm2(V)), fast_amp=float(mp.e**(-2*tau)))
        for name, C0 in states.items():
            Dm, Bm, Fm = L.entropy_and_defect(V, C0, w.conj()); Dz, Bz, Fz = L.entropy_and_defect(Z, C0, w.conj()); De, Be, Fe = L.entropy_and_defect(E, C0, w.conj())
            f = N.formulas(Ks[name], alpha, 1, eps, zs, zf, tau)
            f0 = N.formulas(Ks[name], alpha, 1, eps, zs0, zf0, tau)
            # diagonal theorem transplanted by similarity: |b_f|^2 -> |l_f w|^2 ||v_f||^2 = ||w||^2, no slow-direction share
            nw2 = float(np.linalg.norm(w)**2); ef, es = mp.e**(2*mp.re(zf)*tau), mp.e**(2*mp.re(zs)*tau)
            D_naive = nw2*ef*(-mp.log(ef)) + abs(alpha)**2*eps**2*es*(-mp.log(es))
            rec[name] = dict(D_micro=float(Dm), F_micro=float(Fm), D_pole=float(Dz), F_pole=float(Fz), D_davies=float(De), F_davies=float(Fe),
                             relD_pole=float(abs(Dm/Dz-1)), relF_pole=float(abs(Fm/Fz-1)), D_form=float(f['D']), F_form=float(f['F']), D_form_bare=float(f0['D']), F_form_bare=float(f0['F']),
                             D_fastdir=float(f['D_fastdir']), D_slowdir=float(f['D_slowdir']), F_fast=float(f['F_fast']), F_slow=float(f['F_slow']), D_naive=float(D_naive),
                             slow_part_memory=float(abs(alpha)**2*eps**2*es), slow_part_nonorth=float(abs(Ks[name]['wpar'])**2*ef))
        rows.append(rec)
        json.dump(dict(meta=meta, rows=rows), open('../data/testN2_micro.json', 'w'), indent=1)
    print('  done', lam, '%.0fs' % (time.time()-t0), flush=True)
