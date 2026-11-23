"""Test H: Theorem 6.4 at the level of the exact two-pole model, deep in the asymptotic regime.
Resonances zeta_j and residues Q_j of the cutoff pencil are computed in long double (small circles, large L);
Phi(tau) = e^{zeta_s tau} Q_s + e^{zeta_f tau} Q_f is then evaluated in multiprecision, so both thresholds can be
crossed for m = 1, 2, 3 without the floating-point limitation of direct propagation.  Also reports the residue
asymptotics (Q_s)_{sf} / (h_m eps^m) -> 1."""
import numpy as np, mpmath as mp, json, sys, p030lib as L
mp.mp.dps = 120
ld, cld = np.longdouble, np.clongdouble
PI = np.arccos(ld(-1))
def pole_and_residue(mem, eps, Lc, z0, rad=0.4, n=96):
    """zero (first moment / zeroth moment of 1/det) and residue matrix on a circle around the macroscopic point z0."""
    k = np.arange(n, dtype=ld)+ld(0.5); th = 2*PI*k/n; e = np.cos(th)+cld(1j)*np.sin(th)
    zeta = cld(complex(z0)) + ld(rad)*e; z = zeta*ld(eps); dzeta = cld(1j)*ld(rad)*e*(2*PI/n)
    K = mem.Khat(z, Lc)
    M = zeta[:, None, None]*np.eye(2, dtype=cld)[None] + K
    det = M[:, 0, 0]*M[:, 1, 1]-M[:, 0, 1]*M[:, 1, 0]
    R = np.empty_like(M); R[:, 0, 0] = M[:, 1, 1]/det; R[:, 1, 1] = M[:, 0, 0]/det; R[:, 0, 1] = -M[:, 0, 1]/det; R[:, 1, 0] = -M[:, 1, 0]/det
    Q = np.einsum('z,zij->ij', dzeta, R)/(cld(2j)*PI)
    P1 = np.einsum('z,zij->ij', dzeta*zeta, R)/(cld(2j)*PI)
    zj = np.trace(P1)/np.trace(Q)           # exact for a single simple pole inside (rank-one residue)
    return zj, Q
rows = []; res_rows = []
C0 = np.array([[0.25, 0.08+0.03j], [0.08-0.03j, 0.15]], complex); bb = np.array([0, 0.12], complex)
Ci = np.linalg.inv(C0); css = C0[0, 0].real; cf = 1/Ci[1, 1].real; bf2 = abs(bb[1])**2; d0 = L.dplus_scalar(cf, bf2)
gs, gf = 1.0, 3.0; dg = gf-gs
for m in (1, 2, 3):
    bath = L.bath_poly(m, gs=gs, gf=gf); kap = bath.kappa[m]; h = -(-gs)**m*kap/dg
    mem = L.Memory(bath, 3000.0, dtype=np.longdouble)
    cD, cB = 2*m/dg, 4*m/dg
    for eps in (1e-2, 1e-3, 1e-4):
        Lc = min(1500.0, 3.0/eps)
        zs, Qs = pole_and_residue(mem, eps, Lc, -gs); zf, Qf = pole_and_residue(mem, eps, Lc, -gf)
        zs_m, zf_m = L.cld_to_mp(zs), L.cld_to_mp(zf)
        Qs_m, Qf_m = L.mp_from_cld_matrix(Qs), L.mp_from_cld_matrix(Qf)
        # The residues are exactly symmetric of rank one (simple zeros, Prop. 6.2).  Their smallest entries,
        # (Q_s)_ff = (Q_s)_sf^2/(Q_s)_ss ~ eps^{2m}, are below 80-bit resolution, so impose the rank-one
        # structure from the two accurately computed entries instead of using the noisy small entry.
        rank1_defect = float(abs(Qs_m[0, 0]*Qs_m[1, 1]-Qs_m[0, 1]**2)/abs(Qs_m[0, 1])**2)
        Qs_m = mp.matrix([[Qs_m[0, 0], Qs_m[0, 1]], [Qs_m[0, 1], Qs_m[0, 1]**2/Qs_m[0, 0]]])
        Qf_m = mp.matrix([[Qf_m[0, 1]**2/Qf_m[1, 1], Qf_m[0, 1]], [Qf_m[0, 1], Qf_m[1, 1]]])
        ratio = Qs_m[0, 1]/(mp.mpc(h)*mp.mpf(eps)**m)
        res_rows.append(dict(m=m, eps=eps, L=Lc, zeta_s=[float(mp.re(zs_m)), float(mp.im(zs_m))], zeta_f=[float(mp.re(zf_m)), float(mp.im(zf_m))],
                             Qsf_over_h_eps_m=[float(mp.re(ratio)), float(mp.im(ratio))], abs_ratio_minus_1_over_eps=float(abs(ratio-1)/eps),
                             Qss_minus_1_over_eps=float(abs(Qs_m[0, 0]-1)/eps), rank_one_defect_raw=rank1_defect, S_minus_I_over_eps=float(L.mp_norm2(Qs_m+Qf_m-mp.eye(2))/eps)))
        print(res_rows[-1], flush=True)
        lam = mp.sqrt(mp.mpf(eps))
        for c in [0.5*cD, 0.9*cD, cD, 1.1*cD, 0.5*(cD+cB), 0.9*cB, cB, 1.1*cB, 1.5*cB]:
            tau = c*mp.log(1/lam)
            Phi = mp.e**(zs_m*tau)*Qs_m + mp.e**(zf_m*tau)*Qf_m
            Dz, Bz, Fz = L.entropy_and_defect(Phi, C0, bb, dps=120)
            h2 = abs(h)**2; rf, rs = mp.e**(2*mp.re(zf_m)*tau), mp.e**(2*mp.re(zs_m)*tau)
            ef, es = mp.e**(-2*gf*tau), mp.e**(-2*gs*tau)
            D_ren = bf2*rf*(-mp.log(rf)-mp.log(cf)+d0/bf2) + bf2*h2*mp.mpf(eps)**(2*m)*rs*(-mp.log(rs)-mp.log(css))
            F_ren = d0*rf + bf2**2*h2**2*mp.mpf(eps)**(4*m)*rs/(2*css)
            D_lead = 2*tau*bf2*(gf*ef + gs*h2*mp.mpf(eps)**(2*m)*es)
            F_lead = d0*ef + bf2**2*h2**2*mp.mpf(eps)**(4*m)*es/(2*css)
            rows.append(dict(m=m, eps=eps, c=float(c), c_over_cD=float(c/cD), c_over_cB=float(c/cB), tau=float(tau),
                             D_over_Dren=float(Dz/D_ren), F_over_Fren=float(Fz/F_ren), D_over_Dlead=float(Dz/D_lead), F_over_Flead=float(Fz/F_lead),
                             logD_over_loglam=float(mp.log(Dz)/mp.log(lam)), logF_over_loglam=float(mp.log(Fz)/mp.log(lam)),
                             exponent_D_pred=float(min(2*c*gf, 4*m+2*c*gs)), exponent_F_pred=float(min(2*c*gf, 8*m+2*c*gs)),
                             defect_over_D=float(Fz/Dz)))
        json.dump(dict(residues=res_rows, asymptotics=rows), open('../data/testH_pole_asymptotics.json', 'w'), indent=1)
