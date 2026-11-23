"""Slow-pole residue asymptotics: (Q_s)_{sf} = h_m eps^m (1+O(eps)),  h_m = -kappa_m (zeta_s^0)^m/(zeta_s^0-zeta_f^0),
for m=1,2,3, with and without Lamb shifts / diagonal detunings (zeta_j^0 = -(gamma_j + i omega_j + i d_j))."""
import numpy as np, mpmath as mp, json, p030lib as L
mp.mp.dps = 40
ld, cld = np.longdouble, np.clongdouble
def residue(mem, Dt, eps, Lc, zc, rad, n=96):
    """(1/2 pi i) oint R(z) dz on a small circle (macroscopic center zc, radius rad), long double."""
    k = np.arange(n, dtype=ld)+ld(0.5); th = 2*np.arccos(ld(-1))*k/n
    e = np.cos(th)+cld(1j)*np.sin(th)
    z = (cld(complex(zc)) + ld(rad)*e)*ld(eps); dz = cld(1j)*ld(rad)*ld(eps)*e*(2*np.arccos(ld(-1))/n)
    K = mem.Khat(z, Lc)
    T = z[:, None, None]*np.eye(2, dtype=cld)[None] + cld(1j)*ld(eps)*np.asarray(Dt).astype(cld)[None] + ld(eps)*K
    det = T[:, 0, 0]*T[:, 1, 1]-T[:, 0, 1]*T[:, 1, 0]
    R = np.empty_like(T); R[:, 0, 0] = T[:, 1, 1]/det; R[:, 1, 1] = T[:, 0, 0]/det; R[:, 0, 1] = -T[:, 0, 1]/det; R[:, 1, 0] = -T[:, 1, 0]/det
    return np.einsum('z,zij->ij', dz, R)/(cld(2j)*np.arccos(ld(-1)))
rows = []
for m in (1, 2, 3):
    for variant in ('plain', 'lamb+detuning'):
        gs, gf = 1.0, 2.5
        base = L.bath_poly(m, gs=gs, gf=gf)
        if variant == 'plain':
            bath = base; Dt = np.zeros((2, 2)); ws = wf = 0.0
        else:
            asym_s, asym_f = 0.6, -0.4
            eta = base.eta
            def Jf(x, base=base):
                J = base.Jfun(x); ph = L.bump(x)
                J[..., 0, 0] = J[..., 0, 0] + gs/np.pi*ph*asym_s*x; J[..., 1, 1] = J[..., 1, 1] + gf/np.pi*ph*asym_f*x
                return J
            bath = L.Bath(Jf, 2, 'lamb'); assert bath.min_eig(20001) >= -1e-15
            xg, wg = np.polynomial.legendre.leggauss(600); iphi = float(np.sum(wg*L.bump(xg)))
            ws, wf = -gs/np.pi*asym_s*iphi, -gf/np.pi*asym_f*iphi          # Lamb shifts: Ghat_jj(0) = gamma_j + i omega_j
            Dt = np.diag([0.3, -0.2])
            bath.kappa = {m: base.kappa[m]*eta/base.eta}
        kap = bath.kappa[m]
        zs0 = -(gs + 1j*(ws+Dt[0, 0])); zf0 = -(gf + 1j*(wf+Dt[1, 1]))
        h = -kap*zs0**m/(zs0-zf0)
        mem = L.Memory(bath, 3000.0, dtype=np.longdouble)
        for eps in (1e-2, 3e-3, 1e-3, 3e-4):
            Lc = min(1500.0, 4.0/eps)
            Q = residue(mem, Dt, eps, Lc, zs0, 0.45)
            q_sf = L.cld_to_mp(Q[0, 1]); q_ss = L.cld_to_mp(Q[0, 0]); q_ff = L.cld_to_mp(Q[1, 1])
            rec = dict(m=m, variant=variant, eps=eps, h_m=[h.real, h.imag], ratio=[float(mp.re(q_sf/(h*eps**m))), float(mp.im(q_sf/(h*eps**m)))],
                       ratio_minus_1_over_eps=float(abs(q_sf/(h*eps**m)-1)/eps), Qss_minus_1_over_eps=float(abs(q_ss-1)/eps),
                       rank_one_defect=float(abs(q_ss*q_ff-q_sf**2)/abs(q_sf)**2))
            rec['L'] = Lc; rows.append(rec); print(rec, flush=True)
json.dump(rows, open('../data/testG_residues.json', 'w'), indent=1)
