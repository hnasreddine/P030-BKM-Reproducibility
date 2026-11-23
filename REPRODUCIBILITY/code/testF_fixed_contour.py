"""Counterexample to V0.3's FIXED contour (fixed positive distance from the imaginary axis).
The cutoff pencil T(z)=z+iD+eps*K_L(z), L=eps^{-alpha}, is entire; on Re z=-delta<0 fixed the term eps*K_L(z)
is NOT O(eps): it is astronomically large, and det T has many spurious zeros inside any fixed contour.
The corrected contour (left edge at Re z = -1/L) contains exactly d zeros and gives S = I + O(eps)."""
import numpy as np, json, p030lib as L
b = L.bath_m1(); D = np.zeros((2, 2)); d = 2
mem = L.Memory(b, 420.0)
rows = []
for lam, alpha in [(0.2, 0.5), (0.1, 0.5), (0.05, 0.5), (0.05, 0.75), (0.03, 0.5), (0.03, 0.75)]:
    eps = lam*lam; Lc = eps**(-alpha)
    # (a) size of eps*K on the left edge of a fixed contour (delta = 0.3) vs on the thin contour (1/L)
    y = np.linspace(-1.2, 1.2, 241)
    sizes = {}
    for name, delta in (('fixed_delta_0.3', 0.3), ('thin_1_over_L', 1/Lc)):
        K = mem.Khat(-delta+1j*y, Lc)
        sizes[name] = float(eps*np.max(np.linalg.norm(K, ord=2, axis=(1, 2))))
    # (b) zero count inside rectangles by the argument principle
    def count(re_lo, nleft=30000, nother=3000):
        re_hi, im = 0.5, 1.3
        tl = np.linspace(0, 1, nleft, endpoint=False); to = np.linspace(0, 1, nother, endpoint=False)
        z = np.concatenate([re_lo+(re_hi-re_lo)*to - 1j*im, re_hi + 1j*im*(2*to-1), re_hi-(re_hi-re_lo)*to + 1j*im, re_lo - 1j*im*(2*tl-1), [re_lo-1j*im]])
        wn, det = L.winding_number(mem, D, eps, Lc, z)
        ph = np.angle(det); jump = float(np.max(np.abs(np.angle(det[1:]/det[:-1]))))   # must stay well below pi
        return wn, jump
    wn_thin, j1 = count(-1/Lc)
    wn_fixed, j2 = count(-0.3)
    # (c) moments on the two contours
    def SP(re_lo, ns):
        Z, W = [], []
        cs = [re_lo-1.3j, 0.5-1.3j, 0.5+1.3j, re_lo+1.3j, re_lo-1.3j]
        for a, c in zip(cs[:-1], cs[1:]):
            g, w = np.polynomial.legendre.leggauss(ns); Z.append((a+c)/2+(c-a)/2*g); W.append((c-a)/2*w)
        z, dz = np.concatenate(Z), np.concatenate(W)
        S, P, A = L.contour_moments(mem, D, eps, Lc, z, dz); return S
    # thin rectangle needs fine quadrature near the poles: use the equivalent small circle inside it instead
    z, dz = L.circle_contour(-1.5*eps, 1.2*eps, 400)
    S_thin, _, _ = L.contour_moments(mem, D, eps, Lc, z, dz)
    S_fixed = SP(-0.3, 4000)
    rec = dict(lam=lam, eps=eps, alpha=alpha, L=Lc, eps_K_fixed=sizes['fixed_delta_0.3'], eps_K_thin=sizes['thin_1_over_L'],
               zeros_in_thin=round(wn_thin, 3), zeros_in_fixed=round(wn_fixed, 3), max_phase_step=max(j1, j2),
               S_thin_minus_I=float(np.linalg.norm(S_thin-np.eye(2), 2)), S_thin_minus_I_over_eps=float(np.linalg.norm(S_thin-np.eye(2), 2)/eps),
               S_fixed_minus_I=float(np.linalg.norm(S_fixed-np.eye(2), 2)))
    print(rec, flush=True); rows.append(rec)
json.dump(rows, open('../data/testF_fixed_contour.json', 'w'), indent=1)
