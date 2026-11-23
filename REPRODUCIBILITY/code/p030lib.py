"""P030 independent verification library (clean-room; does not import the V0.3 scripts).

Conventions (rotating frame, x = omega - Delta, eps = lambda^2):
  dV/dT = -i D V - eps * int_0^T F(s) V(T-s) ds,   F(s) = int J(x) exp(-i x s) dx,
  Laplace:  Vhat(z) = [z + iD + eps*Fhat(z)]^{-1},  Fhat(z) = int J(x)/(z+ix) dx  (Re z>0).
Macroscopic time tau = eps*T.
"""
import numpy as np, mpmath as mp
from numpy.polynomial.legendre import leggauss

def mp_to_ld(v):
    hi = float(v); lo = float(v - mp.mpf(hi))
    return np.longdouble(hi) + np.longdouble(lo)

def ld_to_mp(x):
    x = np.longdouble(x); hi = float(x); lo = float(x - np.longdouble(hi))
    return mp.mpf(hi) + mp.mpf(lo)

def cld_to_mp(z):
    return mp.mpc(ld_to_mp(np.real(z)), ld_to_mp(np.imag(z)))

def leggauss_ld(n):
    """Gauss-Legendre nodes/weights refined to long-double accuracy (Newton on P_n)."""
    x0, _ = leggauss(n)
    x = x0.astype(np.longdouble)
    for _ in range(3):
        p0 = np.ones_like(x); p1 = x.copy()
        for k in range(2, n+1):
            p0, p1 = p1, ((2*k-1)*x*p1 - (k-1)*p0)/k
        dp = n*(x*p1 - p0)/(x*x - 1)
        x = x - p1/dp
    p0 = np.ones_like(x); p1 = x.copy()
    for k in range(2, n+1):
        p0, p1 = p1, ((2*k-1)*x*p1 - (k-1)*p0)/k
    dp = n*(x*p1 - p0)/(x*x - 1)
    w = 2/((1 - x*x)*dp*dp)
    return x, w

# ---------------------------------------------------------------- bump functions
def bump(x):
    x = np.asarray(x); 
    if x.dtype != np.longdouble: x = x.astype(float)
    y = np.zeros_like(x)
    m = np.abs(x) < 1
    y[m] = np.exp(1.0 - 1.0/(1.0 - x[m]**2))
    return y

def bump_derivs(t, order):
    """bump and its derivatives up to `order` at points t (|t|<1), via exact recursion.
    bump = exp(1 - 1/(1-t^2)).  Write bump^(k) = P_k(t)/(1-t^2)^(2k) * bump with polynomials P_k:
    P_{k+1} = P_k' (1-t^2)^2 + P_k * (4k t (1-t^2) - 2 t)."""
    t = np.asarray(t)
    if t.dtype != np.longdouble: t = t.astype(float)
    from numpy.polynomial import polynomial as Pl
    one_m = np.array([1.0, 0.0, -1.0])
    P = [np.array([1.0])]
    for k in range(order):
        Pk = P[-1]
        term1 = Pl.polymul(Pl.polyder(Pk) if len(Pk) > 1 else np.array([0.0]), Pl.polymul(one_m, one_m))
        term2 = Pl.polymul(Pk, Pl.polyadd(4*k*Pl.polymul(np.array([0.0, 1.0]), one_m), np.array([0.0, -2.0])))
        P.append(Pl.polyadd(term1, term2))
    out = []
    m = np.abs(t) < 1
    b = bump(t)
    for k in range(order+1):
        y = np.zeros_like(t)
        pv = np.zeros_like(t[m])
        for cf in P[k][::-1]: pv = pv*t[m] + t.dtype.type(cf)
        y[m] = pv / (1 - t[m]**2)**(2*k) * b[m]
        out.append(y)
    return out

# ---------------------------------------------------------------- baths
class Bath:
    """Matrix spectral density J(x) (d x d, PSD), supported in [-1,1]."""
    def __init__(self, Jfun, d, name=""):
        self.Jfun, self.d, self.name = Jfun, d, name
    def J(self, x):
        x = np.asarray(x)
        return self.Jfun(x if x.dtype == np.longdouble else x.astype(float))
    def min_eig(self, n=4001):
        x = np.linspace(-1, 1, n)
        return float(np.min(np.linalg.eigvalsh(self.J(x))))

A_M1 = None
def a_m1():
    global A_M1
    if A_M1 is None:
        x, w = leggauss(400)
        A_M1 = float(np.sum(w*x*x*bump(x))/np.sum(w*bump(x)))
    return A_M1

def bath_m1(gs=1.0, gf=2.0, zeta=0.4):
    """V0.3 Example: J = phi(x) [[gs/pi, zeta x(x^2-a)],[., gf/pi]] ; first off-diagonal jet order m=1."""
    a = a_m1()
    def Jf(x):
        ph = bump(x); J = np.zeros(x.shape+(2, 2), dtype=ph.dtype)
        J[..., 0, 0] = gs/np.pi*ph; J[..., 1, 1] = gf/np.pi*ph
        J[..., 0, 1] = J[..., 1, 0] = zeta*x*(x*x-a)*ph
        return J
    b = Bath(Jf, 2, "m1"); b.kappa = {1: -1j*np.pi*zeta*a}; b.gs, b.gf, b.m = gs, gf, 1
    return b

def q_jet(x, m, x0=0.45, w=0.35):
    """q_m(x) = x^m theta^{(m)}(x), theta(x)=bump((x-x0)/w): int q x^{-k-1}=0 (k<m), = m! int theta x^{-m-1} (k=m)."""
    x = np.asarray(x)
    if x.dtype != np.longdouble: x = x.astype(float)
    t = (x-x0)/w
    return x**m * bump_derivs(t, m)[m] / w**m

def bath_jet(m, eta=None, gs=1.0, gf=2.0, x0=0.45, w=0.35, safety=0.9):
    """Diagonal Davies generator, off-diagonal causal jet of exact order m (explicit construction)."""
    xx = np.linspace(x0-w, x0+w, 200001)[1:-1]
    ratio = np.max(np.abs(q_jet(xx, m, x0, w))/bump(xx))
    eta_max = np.sqrt(gs*gf)/np.pi/ratio
    if eta is None: eta = safety*eta_max
    def Jf(x):
        ph = bump(x); J = np.zeros(x.shape+(2, 2), dtype=ph.dtype)
        J[..., 0, 0] = gs/np.pi*ph; J[..., 1, 1] = gf/np.pi*ph
        J[..., 0, 1] = J[..., 1, 0] = eta*q_jet(x, m, x0, w)
        return J
    b = Bath(Jf, 2, f"jet{m}"); b.eta, b.eta_max, b.gs, b.gf, b.m = eta, eta_max, gs, gf, m
    xg, wg = leggauss(4000); xg = x0 + w*xg; wg = w*wg
    th = bump((xg-x0)/w)
    mu = {k: float(np.sum(wg*q_jet(xg, m, x0, w)/xg**(k+1))) for k in range(m+3)}
    b.mu = mu
    b.mu_exact_m = float(np.math.factorial(m)*np.sum(wg*th/xg**(m+1))) if hasattr(np, "math") else None
    b.kappa = {k: eta*(-1)**k/(1j)**(k+1)*mu[k] for k in mu}
    return b

def bath_jet2(m, gs=1.0, gf=3.0, A=2.0, x0=0.5, w=0.3, wpsi=0.4, safety=0.8):
    """Diagonal Davies generator with rates gs<gf, zero Lamb shifts, off-diagonal causal jet of exact order m:
       J_jj(x) = (g_j/pi) phi(x) + A psi(x),  psi(x) = bump((|x|-x0)/wpsi)   (even, vanishing near x=0),
       J_sf(x) = eta * x^m theta^{(m)}(x),    theta(x) = bump((x-x0)/w),  w < wpsi.
    The off-resonant diagonal weight A*psi leaves gamma_j and the (zero) Lamb shifts unchanged and
    allows an O(1) off-diagonal jet coefficient kappa_m."""
    assert x0 - wpsi > 0 and x0 + wpsi < 1 and w < wpsi
    def psi(x): return bump((np.abs(x)-x0)/wpsi)
    xx = np.linspace(x0-w, x0+w, 400001)[1:-1]
    ratio = np.max(np.abs(q_jet(xx, m, x0, w))/(A*psi(xx)))
    eta_max = 1.0/ratio; eta = safety*eta_max
    def Jf(x):
        ph = bump(x); ps = psi(x); J = np.zeros(x.shape+(2, 2), dtype=ph.dtype)
        J[..., 0, 0] = gs/np.pi*ph + A*ps; J[..., 1, 1] = gf/np.pi*ph + A*ps
        J[..., 0, 1] = J[..., 1, 0] = eta*q_jet(x, m, x0, w)
        return J
    b = Bath(Jf, 2, f"jet2_{m}"); b.eta, b.eta_max, b.gs, b.gf, b.m, b.A = eta, eta_max, gs, gf, m, A
    xg, wg = leggauss(6000); xg = x0 + w*xg; wg = w*wg
    b.mu = {k: float(np.sum(wg*q_jet(xg, m, x0, w)/xg**(k+1))) for k in range(m+3)}
    import math
    b.mu_exact_m = float(math.factorial(m)*np.sum(wg*bump((xg-x0)/w)/xg**(m+1)))
    b.kappa = {k: eta*(-1)**k/(1j)**(k+1)*b.mu[k] for k in b.mu}
    return b

def poly_jet(m, x0=0.5, w=0.4, nq=8000):
    """Real polynomial p of degree <= m, max|p|=1 on supp theta, with int theta p x^{-k-1} dx = 0 (k<m).
    Then mu_m = int theta p x^{-m-1} dx != 0 automatically.  Returns (coefficients low->high, moments)."""
    xg, wg = leggauss(nq); xg = x0 + w*xg; wg = w*wg
    th = bump((xg-x0)/w)
    M = np.array([[np.sum(wg*th*xg**(j-k-1)) for j in range(m+1)] for k in range(m)]) if m > 0 else np.zeros((0, 1))
    if m > 0:
        _, _, vt = np.linalg.svd(M); a = vt[-1]
    else:
        a = np.array([1.0])
    xs = np.linspace(x0-w, x0+w, 20001); pv = np.polynomial.polynomial.polyval(xs, a)
    a = a/np.max(np.abs(pv))
    if np.sum(wg*th*np.polynomial.polynomial.polyval(xg, a)/xg**(m+1)) < 0: a = -a
    mu = [float(np.sum(wg*th*np.polynomial.polynomial.polyval(xg, a)/xg**(k+1))) for k in range(m+3)]
    return a, mu

def bath_poly(m, gs=1.0, gf=3.0, A=2.0, x0=0.5, w=0.4, safety=0.9):
    """Diagonal Davies generator (rates gs<gf, zero Lamb shifts), off-diagonal causal jet of exact order m:
       J_jj(x) = (g_j/pi) phi(x) + A theta(|x|),   theta(x) = bump((x-x0)/w) supported in (x0-w, x0+w) subset (0,1),
       J_sf(x) = eta * theta(x) p_m(x)  for x>0 (zero for x<=0),  eta = safety*A  (so |J_sf| <= safety*A*theta: PSD)."""
    assert x0-w > 0 and x0+w < 1
    a, mu = poly_jet(m, x0, w); eta = safety*A
    def theta(x): return bump((x-x0)/w)
    def Jf(x):
        ph = bump(x); te = bump((np.abs(x)-x0)/w); J = np.zeros(x.shape+(2, 2), dtype=ph.dtype)
        pv = np.zeros_like(x)
        for cf in a[::-1]: pv = pv*x + x.dtype.type(cf)
        J[..., 0, 0] = gs/np.pi*ph + A*te; J[..., 1, 1] = gf/np.pi*ph + A*te
        J[..., 0, 1] = J[..., 1, 0] = eta*theta(x)*pv
        return J
    b = Bath(Jf, 2, f"poly{m}"); b.eta, b.gs, b.gf, b.m, b.A, b.pcoef, b.mu = eta, gs, gf, m, A, a, mu
    b.kappa = {k: eta*(-1)**k/(1j)**(k+1)*mu[k] for k in range(len(mu))}
    b.qfun = lambda x: eta*theta(x)*np.polynomial.polynomial.polyval(np.asarray(x, float), a)
    return b

def bath_ep():
    """J = phi(x) [[2,1],[1,2]]/(2 pi): Davies matrix [[1,.5],[.5,1]] -> EP under detuning."""
    M = np.array([[2., 1.], [1., 2.]])/(2*np.pi)
    b = Bath(lambda x: bump(x)[..., None, None]*M, 2, "ep"); return b

# ---------------------------------------------------------------- memory kernel and transforms
class Memory:
    """F(s) on composite Gauss-Legendre panels in s; cutoff transforms."""
    def __init__(self, bath, Smax, panel=1.0, npan=10, nx=None, dtype=np.float64):
        self.bath = bath; d = bath.d; self.dtype = dtype
        hp = dtype == np.longdouble
        cplx = np.clongdouble if hp else np.complex128
        # trapezoid rule in x: spectrally accurate for C_c^infty densities flat at the band edges
        # (aliasing error ~ |F(2*pi/dx - s)|); Gauss-Legendre is NOT adequate for Gevrey bumps.
        if nx is None: nx = int(2*(Smax + 2500.0)/(2*np.pi)) + 2
        x = np.linspace(-1, 1, nx).astype(dtype) if not hp else (np.arange(nx, dtype=dtype)*(dtype(2)/(nx-1)) - dtype(1))
        wx = np.full(nx, x[1]-x[0], dtype=dtype)
        J = bath.J(x)
        npanels = int(np.ceil(Smax/panel)); edges = np.linspace(0, Smax, npanels+1).astype(dtype)
        g, wg = leggauss_ld(npan) if hp else leggauss(npan)
        s = (edges[:-1, None] + (g[None, :]+1)/2*(edges[1:, None]-edges[:-1, None])).ravel()
        ws = (wg[None, :]/2*(edges[1:, None]-edges[:-1, None])).ravel()
        F = np.empty((s.size, d, d), cplx)
        Jw = (J*wx[:, None, None]).reshape(nx, d*d).astype(cplx)
        ch = 500 if hp else 2000
        for i in range(0, s.size, ch):
            F[i:i+ch] = (np.exp(cplx(-1j)*np.outer(s[i:i+ch], x)) @ Jw).reshape(-1, d, d)
        self.s, self.ws, self.F, self.Smax = s, ws, F, Smax
    @staticmethod
    def chi(t):
        """C-infinity cutoff: 1 on [0,1], 0 on [2,inf)."""
        t = np.asarray(t)
        if t.dtype != np.longdouble: t = t.astype(float)
        out = np.ones_like(t)
        out[t >= 2] = 0
        m = (t > 1) & (t < 2)
        u = t[m]-1
        h1 = np.exp(-1/(1-u)); h2 = np.exp(-1/u)
        out[m] = h1/(h1+h2)
        return out
    def weights(self, L):
        assert 2*L <= self.Smax*(1+1e-12)
        return self.ws*self.chi(self.s/L)
    def Khat(self, z, L, deriv=0):
        """K(z)=int chi(s/L) F(s) e^{-zs} ds  (physical variable z); array of z allowed."""
        z = np.atleast_1d(np.asarray(z)); z = z.astype(self.F.dtype); w = self.weights(L)
        act = w != 0
        E = np.exp(-np.outer(z, self.s[act]))*(w[act]*(-self.s[act])**deriv)[None, :]
        return np.tensordot(E, self.F[act], axes=(1, 0))
    def L1_tail(self, L):
        nrm = np.linalg.norm(np.asarray(self.F, dtype=complex), ord=2, axis=(1, 2))
        return float(np.sum(self.ws*(1-self.chi(self.s/L))*nrm))
    def jets(self, L, order):
        """Taylor coefficients G_k = K^{(k)}(0)/k! of the cutoff transform."""
        return [self.Khat(0.0, L, deriv=k)[0]/self.dtype(int(mp.factorial(k))) for k in range(order+1)]

def rect_contour(re_lo, re_hi, im_lo, im_hi, n_per_unit=None, n_side=None):
    """Positively oriented rectangle, Gauss-Legendre on each side. Returns nodes z and weights dz."""
    corners = [re_lo+1j*im_lo, re_hi+1j*im_lo, re_hi+1j*im_hi, re_lo+1j*im_hi, re_lo+1j*im_lo]
    Z, W = [], []
    for a, b in zip(corners[:-1], corners[1:]):
        n = n_side if n_side else 64
        g, w = leggauss(n)
        Z.append((a+b)/2 + (b-a)/2*g); W.append((b-a)/2*w)
    return np.concatenate(Z), np.concatenate(W)

def circle_contour(center, radius, n):
    th = 2*np.pi*(np.arange(n)+0.5)/n
    z = center + radius*np.exp(1j*th)
    return z, 1j*radius*np.exp(1j*th)*(2*np.pi/n)

def contour_moments(mem, D, eps, L, z, dz):
    """S, P, A = P S^{-1} for the cutoff pencil  T(z)= z + iD + eps*K(z)  (physical variable)."""
    d = mem.bath.d
    K = np.asarray(mem.Khat(z, L), dtype=complex)
    T = z[:, None, None]*np.eye(d)[None] + 1j*D[None] + eps*K
    R = np.linalg.inv(T)
    S = np.einsum('z,zij->ij', dz, R)/(2j*np.pi)
    P = np.einsum('z,zij->ij', dz*z, R)/(2j*np.pi)
    return S, P, P @ np.linalg.inv(S)

def winding_number(mem, D, eps, L, z_closed):
    """Number of zeros of det T inside a closed polygonal/dense curve (argument principle by phase unwrapping)."""
    d = mem.bath.d
    out = []
    for i in range(0, len(z_closed), 400):
        zz = z_closed[i:i+400]
        K = mem.Khat(zz, L)
        T = zz[:, None, None]*np.eye(d)[None] + 1j*D[None] + eps*K
        out.append(np.linalg.det(T))
    det = np.concatenate(out)
    ph = np.unwrap(np.angle(det))
    return (ph[-1]-ph[0])/(2*np.pi), det

# ---------------------------------------------------------------- direct microscopic propagation
class Friedrichs:
    """One-excitation sector of H = H_S + H_B + lambda V with trapezoid-discretized continuum on [-1,1].
    H = [[D, lam*G],[lam*G^dag, X]].  Chebyshev moments m_k = P T_k(H/R) P allow V(T) for all T<Tmax."""
    def __init__(self, bath, D, lam, N, dtype=np.float64, xlim=1.0):
        self.d = d = bath.d; self.lam = lam; self.N = N
        cplx = np.clongdouble if dtype == np.longdouble else np.complex128
        self.cplx = cplx; self.rdt = dtype
        if dtype == np.longdouble and d == 2:
            # long-double grid and closed-form 2x2 square root: the discretized model is then
            # consistent to ~1e-19 (a double-precision grid leaves an N-dependent 1e-14 floor).
            ld = np.longdouble
            x = np.arange(N, dtype=ld)*(ld(2*xlim)/(N-1)) - ld(xlim); dx = x[1]-x[0]
            J = bath.J(x)
            det = np.maximum(J[:, 0, 0]*J[:, 1, 1] - J[:, 0, 1]*J[:, 1, 0], 0); sd = np.sqrt(det)
            den = np.sqrt(J[:, 0, 0] + J[:, 1, 1] + 2*sd); den[den == 0] = 1
            sq = (J + sd[:, None, None]*np.eye(2, dtype=ld))/den[:, None, None]
        else:
            x = np.linspace(-xlim, xlim, N); dx = x[1]-x[0]
            J = bath.J(x)
            w, U = np.linalg.eigh(J)
            w = np.maximum(w, 0)
            sq = np.einsum('nij,nj,nkj->nik', U, np.sqrt(w), U)      # J^{1/2}
        J = np.asarray(J, dtype=float)
        self.x = x.astype(dtype)
        self.G = (np.sqrt(dtype(dx))*sq).astype(dtype)                  # (N,d,d)
        self.D = np.asarray(D).astype(cplx)
        self.dx = float(dx)
        dx = float(dx)
        nG = np.sqrt(np.sum(np.trace(J, axis1=1, axis2=2))*dx)
        self.R = float(xlim + np.linalg.norm(D, 2) + 2*lam*nG + 0.05)
        self.recurrence = 2*np.pi/dx
    def matvec(self, phi, f):
        """(phi: d x d cols; f: N x d x dcols)."""
        lam = self.lam
        phi2 = self.D @ phi + lam*np.einsum('naj,nac->jc', self.G, f)
        f2 = self.x[:, None, None]*f + lam*np.einsum('naj,jc->nac', self.G, phi)
        return phi2, f2
    def moments(self, K):
        d, N = self.d, self.N; R = self.rdt(self.R)
        phi0 = np.eye(d, dtype=self.cplx); f0 = np.zeros((N, d, d), self.cplx)
        m = np.empty((K+1, d, d), self.cplx); m[0] = phi0
        p1, f1 = self.matvec(phi0, f0); p1 = p1/R; f1 = f1/R
        m[1] = p1
        pa, fa, pb, fb = phi0, f0, p1, f1
        for k in range(2, K+1):
            pc, fc = self.matvec(pb, fb)
            pc = 2*pc/R - pa; fc = 2*fc/R - fa
            m[k] = pc
            pa, fa, pb, fb = pb, fb, pc, fc
        self.m = m; self.K = K
        return m
    def V(self, T, dps=30):
        """V(T) = sum_k c_k m_k, c_k=(2-delta_k0)(-i)^k J_k(R T); Bessel functions by mpmath (downward recurrence)."""
        mp.mp.dps = dps
        xarg = mp.mpf(self.R)*mp.mpf(T)
        K = self.K
        need = int(float(xarg) + 12*float(xarg)**(1/3) + 60)
        if need > K: raise ValueError(f"Chebyshev order {K} too small for R*T={float(xarg):.1f} (need {need})")
        # Miller backward recurrence
        start = need + 40 + (need % 2)
        jp, jc = mp.mpf(0), mp.mpf('1e-300')
        vals = [None]*(need+1); norm = mp.mpf(0)
        for k in range(start, 0, -1):
            jm = 2*k/xarg*jc - jp
            jp, jc = jc, jm
            if k-1 <= need: vals[k-1] = jc
            if (k-1) % 2 == 0: norm += (2 if k-1 > 0 else 1)*jc
        Jk = np.array([mp_to_ld(v/norm) for v in vals], dtype=np.longdouble)
        k = np.arange(need+1)
        ph = np.array([1, -1j, -1, 1j], dtype=self.cplx)[k % 4]
        c = (2*Jk*ph).astype(self.cplx); c[0] = c[0]/2
        return np.tensordot(c, self.m[:need+1], axes=(0, 0))

# ---------------------------------------------------------------- entropy functionals (mpmath, exact given V)
def mp_mat(V):
    V = np.asarray(V)
    return mp.matrix([[mp.mpc(mp.mpf(float(np.real(V[i, j]))), mp.mpf(float(np.imag(V[i, j])))) for j in range(V.shape[1])] for i in range(V.shape[0])])

def entropy_and_defect(V, C0, b, dps=80):
    """rho(V) = [[a, b V^dag],[V b^dag, V C0 V^dag]], a = 1 - Tr(V C0 V^dag).
    Returns (D, B, D-B) as mp numbers: D = D(rho||Pi rho), B = BKM midpoint term."""
    mp.mp.dps = dps
    Vm = V if isinstance(V, mp.matrix) else mp_mat(V)
    C0m = C0 if isinstance(C0, mp.matrix) else mp_mat(C0)
    bm = b if isinstance(b, mp.matrix) else mp_mat(np.asarray(b).reshape(1, -1))
    d = Vm.rows
    C = Vm*C0m*Vm.H
    C = (C + C.H)/2
    v = Vm*bm.H
    a = 1 - mp.re(sum(C[i, i] for i in range(d)))
    R = mp.matrix(d+1)
    R[0, 0] = a
    for i in range(d):
        R[0, i+1] = mp.conj(v[i]); R[i+1, 0] = v[i]
        for j in range(d): R[i+1, j+1] = C[i, j]
    ev = mp.eighe(R, eigvals_only=True)
    cs, Q = mp.eighe(C)
    h = lambda t: t*mp.log(t) if t > 0 else mp.mpf(0)
    Dval = mp.fsum(h(mp.re(e)) for e in ev) - h(a) - mp.fsum(h(mp.re(c)) for c in cs)
    Bval = mp.mpf(0)
    for j in range(d):
        pr = mp.fsum(mp.conj(Q[k, j])*v[k] for k in range(d))
        cj = mp.re(cs[j])
        Bval += abs(pr)**2*(mp.log(a)-mp.log(cj))/(a-cj)
    return Dval, Bval, Dval-Bval

def dplus_scalar(c, k):
    c, k = mp.mpf(c), mp.mpf(k)
    return (c-k)*mp.log((c-k)/c) + k


# ---------------------------------------------------------------- high-precision 2x2 contour reconstruction
def contour_Z_hp(mem, D, eps, L, center, radius, n, T_list, dps=40):
    """Long-double contour moments for d=2 (explicit adjugate inverse), then mp for A=P S^{-1} and exp(T A) S.
    center/radius in the physical variable. Returns (S, P, A, [Z(T)]) as mp matrices."""
    assert mem.bath.d == 2
    mp.mp.dps = dps
    ld, cld = np.longdouble, np.clongdouble
    k = np.arange(n, dtype=ld) + ld(0.5)
    th = 2*ld(np.pi)*k/n if False else (2*np.arccos(ld(-1))*k/n)
    e = np.cos(th) + cld(1j)*np.sin(th)
    z = cld(center) + ld(radius)*e
    dz = cld(1j)*ld(radius)*e*(2*np.arccos(ld(-1))/n)
    K = mem.Khat(z, L)
    Tm = z[:, None, None]*np.eye(2, dtype=cld)[None] + cld(1j)*np.asarray(D).astype(cld)[None] + ld(eps)*K
    det = Tm[:, 0, 0]*Tm[:, 1, 1] - Tm[:, 0, 1]*Tm[:, 1, 0]
    R = np.empty_like(Tm)
    R[:, 0, 0] = Tm[:, 1, 1]/det; R[:, 1, 1] = Tm[:, 0, 0]/det
    R[:, 0, 1] = -Tm[:, 0, 1]/det; R[:, 1, 0] = -Tm[:, 1, 0]/det
    twopi_i = cld(2j)*np.arccos(ld(-1))
    S = np.einsum('z,zij->ij', dz, R)/twopi_i
    P = np.einsum('z,zij->ij', dz*z, R)/twopi_i
    Sm = mp.matrix([[cld_to_mp(S[i, j]) for j in range(2)] for i in range(2)])
    Pm = mp.matrix([[cld_to_mp(P[i, j]) for j in range(2)] for i in range(2)])
    Am = Pm*mp.inverse(Sm)
    Zs = [mp.expm(mp.mpf(T)*Am)*Sm for T in T_list]
    return Sm, Pm, Am, Zs

def mp_from_cld_matrix(V):
    return mp.matrix([[cld_to_mp(V[i, j]) for j in range(V.shape[1])] for i in range(V.shape[0])])

def mp_norm2(M):
    return mp.sqrt(max(mp.re(e) for e in mp.eighe(M.H*M, eigvals_only=True)))
