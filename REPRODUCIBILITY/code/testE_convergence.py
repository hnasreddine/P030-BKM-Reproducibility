"""Test E: numerical convergence studies (m=1 bath, lambda=0.1 and 0.07, tau = 2.2 log(1/lambda))."""
import numpy as np, mpmath as mp, json, scipy.linalg as sla, p030lib as L
mp.mp.dps = 50
b = L.bath_m1(); D = np.zeros((2, 2)); out = {}
def nrm(M): return float(L.mp_norm2(M))
for lam in (0.1, 0.07):
    eps = lam*lam; tau = 2.2*np.log(1/lam); T = tau/eps; res = {}
    Kneed = int(1.48*T + 12*(1.48*T)**(1/3) + 150)
    # reference: long double, fine grid
    Nref = int((T + 60/eps)/np.pi) + 2000
    fr = L.Friedrichs(b, D, lam, Nref, dtype=np.longdouble); fr.moments(Kneed)
    Vref = L.mp_from_cld_matrix(fr.V(T, dps=40)); mp.mp.dps = 50
    # (1) bath-frequency discretization
    rows = []
    for N in [int(f*T/np.pi) for f in (1.02, 1.1, 1.3, 1.6, 2.0, 3.0)] + [int((T+45/eps)/np.pi)+300]:
        f2 = L.Friedrichs(b, D, lam, N, dtype=np.longdouble); f2.moments(Kneed)
        V = L.mp_from_cld_matrix(f2.V(T, dps=40)); mp.mp.dps = 50
        rows.append(dict(N=N, recurrence_over_T=float(f2.recurrence/T), diff_vs_ref=nrm(V-Vref)))
    res['bath_discretization'] = rows
    # (2) frequency window: grid extended beyond the support of J (J=0 there) -- must not change anything
    rows = []
    for xlim in (1.0, 1.25, 1.6):
        N = int(xlim*((T + 45/eps)/np.pi + 300))
        f2 = L.Friedrichs(b, D, lam, N, dtype=np.float64, xlim=xlim); f2.moments(int(f2.R*T + 12*(f2.R*T)**(1/3) + 150))
        V = L.mp_mat(np.asarray(f2.V(T), complex)); rows.append(dict(xlim=xlim, N=N, diff_vs_ref=nrm(V-Vref)))
    res['frequency_window'] = rows
    # (3) floating-point precision of the direct propagation
    N = int((T + 45/eps)/np.pi) + 300
    fd = L.Friedrichs(b, D, lam, N, dtype=np.float64); fd.moments(Kneed)
    res['precision_direct'] = dict(double_vs_longdouble=nrm(L.mp_mat(np.asarray(fd.V(T), complex))-Vref))
    # (4) Chebyshev order / spectral bound / dense diagonalization (long-time propagation)
    rows = []
    for Kc in (Kneed, int(1.5*Kneed)):
        f2 = L.Friedrichs(b, D, lam, N, dtype=np.longdouble); f2.moments(Kc)
        rows.append(dict(K=Kc, R=f2.R, diff_vs_ref=nrm(L.mp_from_cld_matrix(f2.V(T, dps=40))-Vref))); mp.mp.dps = 50
    f2 = L.Friedrichs(b, D, lam, N, dtype=np.longdouble); f2.R = 1.3*f2.R; f2.moments(int(1.3*Kneed)+200)
    rows.append(dict(K=f2.K, R=f2.R, diff_vs_ref=nrm(L.mp_from_cld_matrix(f2.V(T, dps=40))-Vref))); mp.mp.dps = 50
    if lam == 0.1:
        Nd = 1400; f3 = L.Friedrichs(b, D, lam, Nd, dtype=np.float64); f3.moments(Kneed)
        H = np.zeros((2+2*Nd, 2+2*Nd)); G = np.zeros((2, 2*Nd))
        for n in range(Nd): G[:, 2*n:2*n+2] = f3.G[n].T
        H[:2, 2:] = lam*G; H[2:, :2] = lam*G.T; H[2:, 2:] = np.diag(np.repeat(f3.x, 2))
        w, U = np.linalg.eigh(H); Pm = U[:2, :]
        Tshort = 0.9*f3.recurrence/1.0 if 0.9*f3.recurrence < T else T
        Vd = (Pm*np.exp(-1j*w*T))@Pm.conj().T
        rows.append(dict(check='dense eigh vs Chebyshev (same N=%d, same T; both include recurrence)' % Nd, diff=float(np.abs(Vd-np.asarray(f3.V(T), complex)).max())))
    res['long_time_propagation'] = rows
    # (5) contour quadrature, (6) memory cutoff, (7) contour shape, (8) precision of the reconstruction
    mem = L.Memory(b, 600.0, dtype=np.longdouble); memd = L.Memory(b, 600.0)
    Lc = min(300.0, 2.0/eps)
    Zref = L.contour_Z_hp(mem, D, eps, Lc, -1.5*eps, 1.2*eps, 320, [T], dps=50)[3][0]
    res['contour_quadrature'] = [dict(n=n, diff_vs_n320=nrm(L.contour_Z_hp(mem, D, eps, Lc, -1.5*eps, 1.2*eps, n, [T], dps=50)[3][0]-Zref)) for n in (16, 24, 32, 48, 64, 96, 128, 160, 240)]
    rows = []
    for Lcut in (25, 50, 75, 100, 150, 200, 250, 300):
        if 2*Lcut*2.7*eps > 12: continue
        Z = L.contour_Z_hp(mem, D, eps, float(Lcut), -1.5*eps, 1.2*eps, 160, [T], dps=50)[3][0]
        rows.append(dict(L=Lcut, L1_tail=mem.L1_tail(float(Lcut)), eps_T_tail=float(eps*T*mem.L1_tail(float(Lcut))), left_edge_times_L=float(2.7*eps*Lcut), err_vs_direct=nrm(Z-Vref), diff_vs_Lmax=nrm(Z-Zref)))
    res['memory_cutoff'] = rows
    rows = []
    for (cen, rad) in ((-1.5, 1.2), (-1.5, 0.9), (-1.2, 1.5), (-1.6, 1.0)):
        Z = L.contour_Z_hp(mem, D, eps, Lc, cen*eps, rad*eps, 240, [T], dps=50)[3][0]
        rows.append(dict(center_over_eps=cen, radius_over_eps=rad, diff=nrm(Z-Zref)))
    zc, dzc = L.circle_contour(-1.5*eps, 1.2*eps, 240)
    S, P, A = L.contour_moments(memd, D, eps, Lc, zc, dzc); Zd = sla.expm(T*A)@S
    res['contour_shape'] = rows
    res['precision_reconstruction'] = dict(double_vs_longdouble=nrm(L.mp_mat(Zd)-Zref))
    res['final'] = dict(T=T, tau=tau, err_pole_vs_direct=nrm(Zref-Vref), normV=nrm(Vref))
    out[str(lam)] = res
    json.dump(out, open('../data/testE_convergence.json', 'w'), indent=1)
    print(lam, json.dumps(res, indent=0)[:3000], flush=True)
