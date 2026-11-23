"""Nonnormal two-channel theorem: shared definitions (clean-room derivation, v0.5).

Convention: Phi(tau) = e^{zeta_s tau} Q_s + e^{zeta_f tau} Q_f,  w = b^dagger (column),  P_s w = 0.
 u      : unit right eigenvector of the slow resonance (limit eps->0)
 l_s    : left (row) eigenvector normalized by l_s u = 1          (P_s = u l_s, NOT an orthogonal projection)
 w_par  : <u, w>            (component of the coherence vector along the slow eigenvector)
 w_perp : || w - u w_par ||
 c_s    : l_s C0 l_s^dagger ;  c_f = det C0 / c_s ;  rho0 = b C0^{-1} b^dagger ;  d0 = c_f * gtilde(rho0)
 alpha  : Q_s(eps) w = alpha eps^m u (1+O(eps));  for m=1: alpha = -(zeta_s/(zeta_s-zeta_f)) l_s G_1 w
Formulas (e_j = |e^{zeta_j tau}|^2, E_j = e^{zeta_j tau}):
 D      ~ w_perp^2 e_f [log(1/e_f) - log c_f + d0/w_perp^2] + |alpha eps^m E_s + w_par E_f|^2 [log(1/e_s) - log c_s]
 D - B  ~ d0 e_f + |alpha|^4 eps^{4m} e_s / (2 c_s)
"""
import numpy as np, mpmath as mp

def gtilde(r): return r + (1-r)*mp.log(1-r)

def davies_data(M0):
    """eigen-data of a 2x2 Davies matrix with distinct real parts: returns zs, zf, u, l_s, v_f, l_f (numpy complex)."""
    ev, V = np.linalg.eig(M0)
    order = np.argsort(-ev.real); ev = ev[order]; V = V[:, order]
    Linv = np.linalg.inv(V)
    u = V[:, 0]/np.linalg.norm(V[:, 0]); ls = Linv[0, :]*np.linalg.norm(V[:, 0])
    vf = V[:, 1]; lf = Linv[1, :]
    return ev[0], ev[1], u, ls, vf, lf

def constants(u, ls, w, C0):
    w = np.asarray(w, complex).reshape(-1); b = w.conj()
    wpar = np.vdot(u, w); wperp2 = float(np.linalg.norm(w - u*wpar)**2)
    cs = float((ls @ C0 @ ls.conj()).real); cf = float(np.linalg.det(C0).real/cs)
    rho0 = float((b @ np.linalg.inv(C0) @ b.conj()).real)
    d0 = cf*gtilde(mp.mpf(rho0))
    return dict(wpar=wpar, wperp2=wperp2, cs=cs, cf=cf, rho0=rho0, d0=d0)

def formulas(K, alpha, m, eps, zs, zf, tau):
    """two-term nonnormal formulas with given (complex) resonances zs, zf (macroscopic units)."""
    tau = mp.mpf(tau); eps = mp.mpf(eps)
    Es, Ef = mp.e**(mp.mpc(zs)*tau), mp.e**(mp.mpc(zf)*tau); es, ef = abs(Es)**2, abs(Ef)**2
    slow_amp2 = abs(mp.mpc(alpha)*eps**m*Es + mp.mpc(complex(K['wpar']))*Ef)**2
    D = K['wperp2']*ef*(-mp.log(ef) - mp.log(K['cf']) + K['d0']/K['wperp2']) + slow_amp2*(-mp.log(es) - mp.log(K['cs']))
    F = K['d0']*ef + abs(mp.mpc(alpha))**4*eps**(4*m)*es/(2*K['cs'])
    D_fastdir = K['wperp2']*ef*(-mp.log(ef)); D_slowdir = slow_amp2*(-mp.log(es))
    return dict(D=D, F=F, D_fastdir=D_fastdir, D_slowdir=D_slowdir, F_fast=K['d0']*ef, F_slow=abs(mp.mpc(alpha))**4*eps**(4*m)*es/(2*K['cs']), slow_amp2=slow_amp2)
