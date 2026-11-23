"""Test N3: generic nonnormal two-channel baths at the level of the reconstructed two-pole dynamics.
 J(x) = phi(x) B(x) B(x)^T / pi with B(x)=B0+x B1 (real 2x2), D = eps*Dt with a complex Hermitian detuning Dt:
 complex resonances (omega_s != omega_f), nonorthogonal eigenvectors, generic first jet.  w := fast right eigenvector.
 Checks (i) alpha = -(zeta_s/(zeta_s-zeta_f)) l_s G_1 w against the numerically computed slow residue,
 (ii) the two formulas against the exact two-pole values across both horizons, (iii) the phase-sensitive interference at c_D."""
import numpy as np, mpmath as mp, json, p030lib as L, nonnormal as N
mp.mp.dps = 160
ld, cld = np.longdouble, np.clongdouble; PI = np.arccos(ld(-1))
def residue(mem, Dt, eps, Lc, z0, rad, n=128):
    k = np.arange(n, dtype=ld)+ld(0.5); th = 2*PI*k/n; e = np.cos(th)+cld(1j)*np.sin(th)
    zeta = cld(complex(z0)) + ld(rad)*e; dzeta = cld(1j)*ld(rad)*e*(2*PI/n)
    K = mem.Khat(zeta*ld(eps), Lc)
    M = zeta[:, None, None]*np.eye(2, dtype=cld)[None] + cld(1j)*Dt.astype(cld)[None] + K
    det = M[:, 0, 0]*M[:, 1, 1]-M[:, 0, 1]*M[:, 1, 0]
    R = np.empty_like(M); R[:, 0, 0] = M[:, 1, 1]/det; R[:, 1, 1] = M[:, 0, 0]/det; R[:, 0, 1] = -M[:, 0, 1]/det; R[:, 1, 0] = -M[:, 1, 0]/det
    Q = np.einsum('z,zij->ij', dzeta, R)/(cld(2j)*PI); P1 = np.einsum('z,zij->ij', dzeta*zeta, R)/(cld(2j)*PI)
    return np.trace(P1)/np.trace(Q), Q
def rank1(Q):
    """exact rank-one completion from the dominant column and row (entries below 80-bit resolution are not used)."""
    Qm = L.mp_from_cld_matrix(Q); i, j = max(((a, b) for a in range(2) for b in range(2)), key=lambda t: abs(Qm[t[0], t[1]]))
    col = Qm[:, j]; row = Qm[i, :]
    return col*row/Qm[i, j]
rng = np.random.default_rng(11); rows = []; res = []
for case in range(3):
  while True:
    B0 = rng.normal(size=(2, 2)); B1 = 0.5*rng.normal(size=(2, 2))
    K0 = 0.5*B0@B0.T; Dt_try = np.array([[0.4, 0.3-0.5j], [0.3+0.5j, -0.6]])*[1.0, 0.5, 1.5][case]
    e0 = np.linalg.eigvals(-(1j*Dt_try+K0))          # Lamb shifts vanish for this even phi only if B1=0; this is a pre-screen
    if abs(e0[0].real-e0[1].real) > 0.5 and max(e0.real) < -0.15: break
  if True:
    def Jf(x, B0=B0, B1=B1):
        ph = L.bump(x); Bx = B0[None].astype(ph.dtype) + x[..., None, None]*B1[None].astype(ph.dtype)
        return (0.5/np.pi)*ph[..., None, None]*np.einsum('...ik,...jk->...ij', Bx, Bx)
    bath = L.Bath(Jf, 2, 'gen')
    Dt = np.array([[0.4, 0.3-0.5j], [0.3+0.5j, -0.6]])*[1.0, 0.5, 1.5][case]
    mem = L.Memory(bath, 3000.0, dtype=np.longdouble)
    G = [np.asarray(g, complex) for g in mem.jets(1500.0, 1)]
    M0 = -(1j*Dt + G[0]); zs0, zf0, u, ls, vf, lf = N.davies_data(M0)
    assert zs0.real - zf0.real > 0.2 and zs0.real < -0.05, (zs0, zf0)
    w = 0.1*vf/np.linalg.norm(vf)
    alpha = -(zs0/(zs0-zf0))*(ls@G[1]@w)
    C0 = np.array([[0.25, 0.08+0.03j], [0.08-0.03j, 0.15]]); K = N.constants(u, ls, w, C0)
    gs, gf = -zs0.real, -zf0.real; cD, cB = 2/(gf-gs), 4/(gf-gs)
    nonorth = abs(np.vdot(u, vf))/np.linalg.norm(vf)
    for eps in (1e-2, 1e-3, 1e-4, 1e-5):
        Lc = min(1500.0, 3.0/eps)
        zs, Qs = residue(mem, Dt, eps, Lc, zs0, 0.25*abs(zs0-zf0)); zf, Qf = residue(mem, Dt, eps, Lc, zf0, 0.25*abs(zs0-zf0))
        Qs_m, Qf_m = rank1(Qs), rank1(Qf); zs_m, zf_m = L.cld_to_mp(zs), L.cld_to_mp(zf)
        wm = mp.matrix([[mp.mpc(complex(x))] for x in w]); um = mp.matrix([[mp.mpc(complex(x)) for x in u.conj()]])
        a_num = (um*(Qs_m*wm))[0]                        # <u, Q_s w>
        res.append(dict(case=case, eps=eps, zs0=[zs0.real, zs0.imag], zf0=[zf0.real, zf0.imag], eigvec_overlap=float(nonorth), alpha=[alpha.real, alpha.imag],
                        a_num_over_alpha_eps=[float(mp.re(a_num/(mp.mpc(alpha)*eps))), float(mp.im(a_num/(mp.mpc(alpha)*eps)))]))
        print(res[-1], flush=True)
        lam = mp.sqrt(mp.mpf(eps))
        for c in [0.5*cD, 0.9*cD, cD, 1.1*cD, 0.5*(cD+cB), cB, 1.5*cB]:
            tau = c*mp.log(1/lam)
            Phi = mp.e**(zs_m*tau)*Qs_m + mp.e**(zf_m*tau)*Qf_m
            Dz, Bz, Fz = L.entropy_and_defect(Phi, C0, w.conj(), dps=160)
            f = N.formulas(K, alpha, 1, eps, complex(zs_m), complex(zf_m), tau)
            # same formula with the interference cross term dropped (incoherent sum) -- to show the cross term is needed at c_D
            Es, Ef = mp.e**(zs_m*tau), mp.e**(zf_m*tau)
            incoh = f['D'] - (f['slow_amp2'] - abs(mp.mpc(alpha))**2*mp.mpf(eps)**2*abs(Es)**2 - abs(mp.mpc(complex(K['wpar'])))**2*abs(Ef)**2)*(-mp.log(abs(Es)**2) - mp.log(K['cs']))
            rows.append(dict(case=case, eps=eps, c_over_cD=float(c/cD), c_over_cB=float(c/cB), D_over_formula=float(Dz/f['D']), F_over_formula=float(Fz/f['F']),
                             D_over_incoherent=float(Dz/incoh), logD_over_loglam=float(mp.log(Dz)/mp.log(lam)), pred_D=float(min(2*c*gf, 4+2*c*gs)),
                             logF_over_loglam=float(mp.log(Fz)/mp.log(lam)), pred_F=float(min(2*c*gf, 8+2*c*gs))))
        json.dump(dict(residues=res, asymptotics=rows), open('../data/testN3_generic.json', 'w'), indent=1)
