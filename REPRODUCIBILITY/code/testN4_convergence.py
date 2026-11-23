"""Test N4: convergence checks for the nonnormal microscopic test (lambda=0.1 and 0.07, c=2 and c=4)."""
import numpy as np, mpmath as mp, json, scipy.linalg as sla, p030lib as L, nonnormal as N
mp.mp.dps = 60
q = 0.5; M0 = np.array([[-1, q], [0, -2.]], complex); Gam = -(M0+M0.conj().T).real; H = 1j*(M0-M0.conj().T)/2
bath = L.Bath(lambda x: L.bump(x)[..., None, None]*(Gam/(2*np.pi)).astype(L.bump(x).dtype), 2, 'nonnormal')
w = 0.1*np.array([-q, 1.0], complex); C0 = 0.2*np.eye(2, dtype=complex)
mem = L.Memory(bath, 600.0, dtype=np.longdouble); memd = L.Memory(bath, 600.0)
nrm = lambda M: float(L.mp_norm2(M)); out = {}
for lam in (0.1, 0.07):
    eps = lam*lam; D = eps*H; res = {}
    for c in (2.0, 4.0):
        tau = c*np.log(1/lam); T = tau/eps; r = {}
        def direct(Nn, dtype=np.longdouble, Kfac=1.0):
            fr = L.Friedrichs(bath, D, lam, Nn, dtype=dtype); fr.moments(int(Kfac*(fr.R*T + 12*(fr.R*T)**(1/3) + 150)))
            V = fr.V(T, dps=40); mp.mp.dps = 60
            return L.mp_from_cld_matrix(V) if dtype == np.longdouble else L.mp_mat(np.asarray(V, complex))
        Nref = int((T + 60/eps)/np.pi) + 2000; Vref = direct(Nref)
        Dref, Bref, Fref = L.entropy_and_defect(Vref, C0, w.conj(), dps=60)
        rows = []
        for fac in (1.1, 1.5, 2.0, 3.0):
            V = direct(int(fac*T/np.pi)); Dm, Bm, Fm = L.entropy_and_defect(V, C0, w.conj(), dps=60)
            rows.append(dict(recurrence_over_T=fac, dV=nrm(V-Vref), relD=float(abs(Dm/Dref-1)), relF=float(abs(Fm/Fref-1))))
        r['bath_discretization'] = rows
        Vd = direct(int((T+45/eps)/np.pi)+300, dtype=np.float64); Dm, Bm, Fm = L.entropy_and_defect(Vd, C0, w.conj(), dps=60)
        r['precision_direct_double_vs_80bit'] = dict(dV=nrm(Vd-Vref), relD=float(abs(Dm/Dref-1)), relF=float(abs(Fm/Fref-1)))
        r['chebyshev_order_x1.5'] = nrm(direct(int((T+45/eps)/np.pi)+300, Kfac=1.5)-Vref)
        Lc = min(300.0, 2.0/eps)
        Zref = L.contour_Z_hp(mem, D, eps, Lc, -1.5*eps, 1.2*eps, 320, [T], dps=60)[3][0]
        Dz, Bz, Fz = L.entropy_and_defect(Zref, C0, w.conj(), dps=60)
        r['contour_nodes'] = [dict(n=n, dZ=nrm(L.contour_Z_hp(mem, D, eps, Lc, -1.5*eps, 1.2*eps, n, [T], dps=60)[3][0]-Zref)) for n in (16, 24, 32, 48, 96)]
        rows = []
        for Lcut in (25, 50, 100, 150, 200, 300):
            if 2*Lcut*2.7*eps > 12: continue
            Z = L.contour_Z_hp(mem, D, eps, float(Lcut), -1.5*eps, 1.2*eps, 160, [T], dps=60)[3][0]; Dm, Bm, Fm = L.entropy_and_defect(Z, C0, w.conj(), dps=60)
            rows.append(dict(L=Lcut, err_vs_direct=nrm(Z-Vref), relD_vs_direct=float(abs(Dm/Dref-1)), relF_vs_direct=float(abs(Fm/Fref-1))))
        r['memory_cutoff'] = rows
        zc, dzc = L.circle_contour(-1.5*eps, 1.2*eps, 240); S, P, A = L.contour_moments(memd, D, eps, Lc, zc, dzc)
        r['precision_reconstruction_double_vs_80bit'] = nrm(L.mp_mat(sla.expm(T*A)@S)-Zref)
        r['final'] = dict(T=T, err_pole=nrm(Zref-Vref), relD=float(abs(Dz/Dref-1)), relF=float(abs(Fz/Fref-1)), D=float(Dref), defect=float(Fref))
        res['c=%g' % c] = r
    out[str(lam)] = res; json.dump(out, open('../data/testN4_convergence.json', 'w'), indent=1); print(lam, 'done', flush=True)
