"""Test D: several Bohr-frequency groups coupled through one common bath (d=4), gaps from O(1) down to 0.
Direct propagation (double) vs grouped contour reconstruction on ONE contour enclosing all 4 resonances."""
import numpy as np, scipy.linalg as sla, json, time, p030lib as L
rng = np.random.default_rng(2030)
B0 = rng.normal(size=(4, 4)); B1 = 0.6*rng.normal(size=(4, 4))
def Jf(x):
    ph = L.bump(x).astype(float)
    Bx = B0[None] + np.asarray(x, float)[..., None, None]*B1[None]
    return 0.1/np.pi*ph[..., None, None]*np.einsum('...ik,...jk->...ij', Bx, Bx)
bath = L.Bath(Jf, 4, 'cluster4')
print('min eig J', bath.min_eig(4001), 'Davies K0 eig', np.linalg.eigvalsh(np.pi*Jf(np.array([0.0]))[0]))
mem = L.Memory(bath, 520.0)
rows = []
for lam in [0.1, 0.07, 0.05, 0.035]:
    eps = lam*lam
    for gap_label, gap in [('0.60', 0.60), ('0.10', 0.10), ('0.02', 0.02), ('3eps', 3*eps), ('0.3eps', 0.3*eps), ('0', 0.0)]:
        t0 = time.time()
        Dd = np.array([-gap/2, -gap/2+0.7*eps, gap/2, gap/2-1.3*eps])
        D = np.diag(Dd)
        c = 2.2; tau = c*np.log(1/lam); T = tau/eps
        N = int((T + 40/eps)/np.pi) + 300
        fr = L.Friedrichs(bath, D, lam, N); fr.moments(int(fr.R*T + 12*(fr.R*T)**(1/3) + 150))
        V = np.asarray(fr.V(T), complex)
        # one ellipse around the segment -i[min d, max d], shifted into the left half plane by the mean damping
        K0 = np.pi*Jf(np.array([0.0]))[0]; gbar = np.trace(K0)/4; gmax = np.linalg.eigvalsh(K0).max()
        cx = -eps*gmax/2; rx = eps*(gmax/2+0.5); ry = (Dd.max()-Dd.min())/2 + rx
        Lc = min(260.0, 1.5/(rx-cx))
        # thin rectangle [cx-rx, cx+rx] x [-dmax-2rx, -dmin+2rx], composite Gauss-Legendre (panels of length ~rx/2)
        g16, w16 = np.polynomial.legendre.leggauss(16)
        def side(p, q):
            npan = max(1, int(np.ceil(abs(q-p)/(rx/2)))); ed = p + (q-p)*np.linspace(0, 1, npan+1)
            zz = ((ed[:-1, None]+ed[1:, None])/2 + (ed[1:, None]-ed[:-1, None])/2*g16[None, :]).ravel()
            ww = ((ed[1:, None]-ed[:-1, None])/2*w16[None, :]).ravel(); return zz, ww
        lo, hi = -Dd.max()-2*rx, -Dd.min()+2*rx
        corners = [cx-rx+1j*lo, cx+rx+1j*lo, cx+rx+1j*hi, cx-rx+1j*hi, cx-rx+1j*lo]
        parts = [side(p, q) for p, q in zip(corners[:-1], corners[1:])]
        z = np.concatenate([p[0] for p in parts]); dz = np.concatenate([p[1] for p in parts]); nq = len(z)
        S = np.zeros((4, 4), complex); P = np.zeros((4, 4), complex)
        for i in range(0, nq, 1500):
            zz, dd = z[i:i+1500], dz[i:i+1500]
            Kz = mem.Khat(zz, Lc)
            R = np.linalg.inv(zz[:, None, None]*np.eye(4)[None] + 1j*D[None] + eps*Kz)
            S += np.einsum('z,zij->ij', dd, R)/(2j*np.pi); P += np.einsum('z,zij->ij', dd*zz, R)/(2j*np.pi)
        A = P@np.linalg.inv(S)
        Z = sla.expm(T*A)@S
        E = sla.expm(T*(-1j*D - eps*(K0 + 0j)))   # Davies-level comparison without Lamb shifts (for scale only)
        ev = np.linalg.eigvals(A)
        rec = dict(lam=lam, gap=gap_label, gap_value=gap, tau=tau, T=T, N=N, nq=nq, L=Lc, left_edge_times_L=float((rx-cx)*Lc),
                   err_pole=float(np.linalg.norm(V-Z, 2)), normV=float(np.linalg.norm(V, 2)), rel_err=float(np.linalg.norm(V-Z, 2)/np.linalg.norm(V, 2)),
                   S_minus_I_over_eps=float(np.linalg.norm(S-np.eye(4), 2)/eps), A_plus_iD_over_eps=float(np.linalg.norm(A+1j*D, 2)/eps),
                   max_re_eig_over_eps=float(ev.real.max()/eps), min_re_eig_over_eps=float(ev.real.min()/eps))
        rows.append(rec); print({k: ('%.4g' % v if isinstance(v, float) else v) for k, v in rec.items()}, '%.0fs' % (time.time()-t0), flush=True)
        json.dump(rows, open('../data/testD_clusters.json', 'w'), indent=1)
