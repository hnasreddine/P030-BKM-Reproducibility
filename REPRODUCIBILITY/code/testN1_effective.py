"""Test N1 (no bath): the mission's counterexample M0=[[-1,q],[0,-2]], q=1/2, w=b0(-q,1)^T, P_s w=0.
(a) pure Davies semigroup: entropy coefficient is NOT the diagonal one after a similarity transformation.
(b) spectral projector vs orthogonal projection: w orthogonal to v_s but P_s w != 0 => no threshold at all.
(c) synthetic two-pole family with slow amplitude alpha*eps^m: formulas vs exact, general C0, complex phases."""
import numpy as np, mpmath as mp, json, p030lib as L, nonnormal as N
mp.mp.dps = 500
q = 0.5; M0 = np.array([[-1, q], [0, -2.]], complex); b0 = 0.1
zs, zf, u, ls, vf, lf = N.davies_data(M0)
out = dict(P_s=np.outer(u, ls).real.tolist(), P_s_is_orthogonal=bool(np.allclose(np.outer(u, ls), np.outer(u, ls).conj().T)))
rows = []
for name, C0 in (('pI', 0.2*np.eye(2, dtype=complex)), ('general', np.array([[0.25, 0.08+0.03j], [0.08-0.03j, 0.15]]))):
    w = b0*np.array([-q, 1.0], complex); K = N.constants(u, ls, w, C0)
    for tau in (5, 10, 20, 40):
        Z = mp.expm(mp.mpf(tau)*mp.matrix(M0.tolist()))
        D, B, F = L.entropy_and_defect(Z, C0, w.conj(), dps=500)
        f = N.formulas(K, 0.0, 1, 1e-3, zs, zf, tau)
        naive = 2*tau*2.0*np.linalg.norm(w)**2*mp.e**(-4*tau)          # diagonal formula transplanted: 2 tau gamma_f |b|^2 e^{-2 gamma_f tau}
        lead = 2*tau*mp.e**(-4*tau)*(1.0*abs(K['wpar'])**2 + 2.0*K['wperp2'])
        rows.append(dict(C0=name, tau=tau, D_over_formula=float(D/f['D']), D_over_leading=float(D/lead), D_over_naive_diagonal=float(D/naive), F_over_formula=float(F/f['F'])))
out['a_davies_semigroup'] = rows
# (b)
w_orth = b0*np.array([0.0, 1.0], complex)       # orthogonal to v_s=(1,0) but l_s w = q b0 != 0
rows = []
for tau in (5, 10, 20):
    Z = mp.expm(mp.mpf(tau)*mp.matrix(M0.tolist())); D, B, F = L.entropy_and_defect(Z, 0.2*np.eye(2), w_orth.conj(), dps=500)
    rows.append(dict(tau=tau, l_s_w=float((ls@w_orth).real), logD_over_tau=float(mp.log(D)/tau), fast_rate_would_be=-4.0, slow_rate=-2.0))
out['b_orthogonal_but_not_dark'] = rows
# (c) synthetic nonsymmetric two-pole family with phases
rng = np.random.default_rng(7); rows = []
for trial in range(4):
    m = [1, 1, 2, 3][trial]
    gs, gf = 1.0, [2.0, 3.0, 2.5, 3.0][trial]; os_, of_ = rng.normal(size=2)
    Vm = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)); Vm[:, 0] /= np.linalg.norm(Vm[:, 0]); Lm = np.linalg.inv(Vm)
    u, ls, vf, lf = Vm[:, 0], Lm[0, :], Vm[:, 1], Lm[1, :]
    w = 0.1*vf/np.linalg.norm(vf); alpha = (0.7-0.4j)
    A0 = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)); C0 = A0@A0.conj().T; C0 = 0.4*C0/np.trace(C0).real
    K = N.constants(u, ls, w, C0); assert K['rho0'] < 0.6
    cD, cB = 2*m/(gf-gs), 4*m/(gf-gs)
    for eps in (1e-2, 1e-4, 1e-6):
        if eps**m < 1e-13: continue    # the float-built w is dark only to ~1e-17; keep alpha*eps^m well above that
        lam = np.sqrt(eps)
        for c in (0.5*cD, cD, 0.5*(cD+cB), cB, 1.5*cB):
            tau = c*np.log(1/lam)
            zs, zf = -gs+1j*os_, -gf+1j*of_
            # exact two-pole matrix: Q_s = (u + eps*du)(l_s + alpha eps^m/(l_f w) l_f), Q_f = v_f l_f  (bounded O(eps) deformations)
            um = mp.matrix([[mp.mpc(complex(x))] for x in u]); lsm = mp.matrix([[mp.mpc(complex(x)) for x in ls]]); vfm = mp.matrix([[mp.mpc(complex(x))] for x in vf]); lfm = mp.matrix([[mp.mpc(complex(x)) for x in lf]])
            rs = lsm + mp.mpc(alpha)*mp.mpf(eps)**m/mp.mpc(complex(lf@w))*lfm
            Phi = mp.e**(mp.mpc(zs)*tau)*um*rs + mp.e**(mp.mpc(zf)*tau)*vfm*lfm
            if L.mp_norm2(Phi) > 1: continue      # synthetic nonnormal family is not a contraction at short times; skip unphysical points
            D, B, F = L.entropy_and_defect(Phi, C0, w.conj(), dps=500)
            f = N.formulas(K, alpha, m, eps, zs, zf, tau)
            rows.append(dict(trial=trial, m=m, eps=eps, c_over_cD=float(c/cD), c_over_cB=float(c/cB), D_over_formula=float(D/f['D']), F_over_formula=float(F/f['F']),
                             logD_over_loglam=float(mp.log(D)/mp.log(lam)), pred_D=float(min(2*c*gf, 4*m+2*c*gs)), logF_over_loglam=float(mp.log(F)/mp.log(lam)), pred_F=float(min(2*c*gf, 8*m+2*c*gs))))
out['c_synthetic_two_pole'] = rows
json.dump(out, open('../data/testN1_effective.json', 'w'), indent=1)
print('P_s', out['P_s'], 'orthogonal?', out['P_s_is_orthogonal'])
for r in out['a_davies_semigroup']: print(r)
for r in out['b_orthogonal_but_not_dark']: print(r)
for r in out['c_synthetic_two_pole']:
    if r['eps'] in (1e-4, 1e-6): print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
