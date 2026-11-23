import numpy as np, time, scipy.linalg as sla, p030lib as L
from scipy.sparse.linalg import expm_multiply
from scipy.sparse import csc_matrix
b = L.bath_m1(); lam = 0.2; D = np.zeros((2,2))
for dt in (np.float64, np.longdouble):
    fr = L.Friedrichs(b, D, lam, 1500, dtype=dt)
    t=time.time(); fr.moments(900); print(dt.__name__, 'moments time', time.time()-t, 'R', fr.R, 'recurrence', fr.recurrence)
    for T in (10., 100., 400.):
        V = fr.V(T)
        # dense reference
        N=fr.N; H = np.zeros((2+2*N, 2+2*N))
        G = np.asarray(fr.G, float).transpose(1,0,2)  # not used
        Gm = np.zeros((2, 2*N))
        for n in range(N): Gm[:, 2*n:2*n+2] = np.asarray(fr.G[n], float).T
        H[:2, 2:] = lam*Gm; H[2:, :2] = lam*Gm.T; H[2:, 2:] = np.diag(np.repeat(np.asarray(fr.x, float), 2))
        if T == 10.:
            w, U = np.linalg.eigh(H); P = U[:2, :]
        Vref = (P*np.exp(-1j*w*T)) @ P.conj().T
        print(' T', T, 'max diff vs eigh', np.abs(np.asarray(V, complex)-Vref).max(), 'V_sf', complex(V[0,1]))
