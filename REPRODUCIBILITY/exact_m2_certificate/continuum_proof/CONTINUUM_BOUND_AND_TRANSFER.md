# P030 — Certified continuum-to-bath-discretization bound for the m=2 example

**Scope:** one-excitation, zero-temperature, rotating-wave Friedrichs Hamiltonian; smooth compactly supported matrix-valued spectral density; \(\lambda=0.05\), \(c=5\), \(N=8500\), \(\gamma_s=1\), \(\gamma_f=3\). Date: 21 November 2026.

## Result 1 — a rigorous discretization bound for the survival propagator

For the **exact** continuous-spectrum survival propagator \(V(T)\) and the **exactly propagated** \(N\)-mode Hermitian bath discretization \(V_N(T)\),

\[
 T=\frac{5\log20}{0.05^2}=5991.464547107981986\ldots,\qquad
 P=\pi(N-1)=26700.39596285965\ldots,
\]

we establish

\[
 \boxed{\|V(T)-V_{8500}(T)\|_2<9.356\times10^{-23}<10^{-22}.}
\]

This is a **continuum-discretization theorem**; the finite-precision propagation is bounded separately by the accompanying numerical certificate. The proof is independent of any spectral gap, singular-value separation, or two-pole approximation. It applies both to the mathematically exact normalized \(m=2\) orthogonal polynomial and to the original **rounded** degree-two polynomial in P030 V0.8, regarded as a precisely defined real polynomial. The exact order-two moment conditions hold only for the first, not automatically for the historical rounded polynomial.

### 1. Matrix Volterra perturbation

For \(\varepsilon=\lambda^2\) and correlations \(F(s)=\int J(x)e^{-ixs}\,dx\), \(F_N(s)=h\sum_{j=0}^{N-1}J(x_j)e^{-ix_js}\), the difference of the causal Volterra resolvents is the double convolution

\[
 V-V_N=-\varepsilon\,V*(F-F_N)*V_N.
\]

Both \(V\) and \(V_N\) are compressions of unitary one-excitation Hamiltonian groups, hence contractions. Consequently

\[
 \|V(T)-V_N(T)\|\le\varepsilon\int_0^T(T-s)\|F(s)-F_N(s)\|\,ds.
\]

### 2. Poisson aliasing and derivatives

The bump density and all its derivatives vanish at the endpoints \([-1,1]\), so trapezoidal sampling with spacing \(h=2/(N-1)\) is a full lattice sum with zero exterior terms. Poisson summation gives \(F_N(s)=\sum_{k\in\mathbb Z}F(s+kP)\), \(P=2\pi/h\), with the Fourier convention above. Repeated integration by parts yields

\[
 \|F(t)\|\le C_r/|t|^r,\quad C_r=\int_{-1}^1\|J^{(r)}(x)\|_2dx.
\]

For \(r>1\), \(P>T\),

\[
 \|V(T)-V_N(T)\|\le\varepsilon T^2\frac{\zeta(r)C_r}{(P-T)^r}.
\]

### 3. Certifying the bump derivative norms without numerical quadrature

Define \(b(t)=\exp[1-(1-t^2)^{-1}]\) for \(|t|<1\), and zero elsewhere. Its exact derivative numerator polynomials are \(Q_0(t)=1\) and

\[
 Q_{k+1}=Q_k'(1-t^2)^2+Q_k\bigl(4kt(1-t^2)-2t\bigr),\qquad
 b^{(k)}(t)=\frac{Q_k(t)}{(1-t^2)^{2k}}b(t).
\]

Every real root of \(Q_k\) inside \((-1,1)\) is isolated using SymPy's exact rational Sturm sequences. Between successive roots \(b^{(k-1)}\) is monotone. Therefore \(B_k:=\int_{-1}^1|b^{(k)}|\,dt\) equals the **total variation** of \(b^{(k-1)}\), evaluated at the isolated roots and at the flat endpoints (where derivatives vanish). Endpoint evaluations and differences are enclosed with `mpmath.iv` outward interval arithmetic, giving safe rational-decimal upper bounds:

| Derivative order | Isolated roots | Certified safe upper bound on \(B_k\) |
|---:|---:|---:|
| 12 | 20 | \(4.378\times10^{21}\) |
| 13 | 23 | \(2.412\times10^{24}\) |
| 14 | 24 | \(1.566\times10^{27}\) |

Reproduce with `full_certified_operator_bound.py`; it calculates the interval enclosures, checks that the stated decimal bounds dominate them, then computes the operator bound.

### 4. Rigorous polynomial bounds without depending on the SVD's coefficients

The original `poly_jet` normalizes its degree-two polynomial to maximum absolute value 1 over 20,001 equally spaced points of \([0.1,0.9]\). Markov's inequality on an interval of length 0.8 implies that its true supremum differs from the sampled supremum by at most a factor \((1-10\times0.00002)^{-1}<1.001\). Therefore, conservatively,

\[
 \|p\|_\infty\le1.001,\quad\|p'\|_\infty\le10.01,\quad\|p''\|_\infty\le25.025.
\]

These bounds hold regardless of the rounding in the three SVD coefficients; even the original rounded density is positive semidefinite because its off-diagonal term obeys \(|1.8\,\vartheta p|<2\vartheta\), while each diagonal spectral density contains at least \(2\vartheta\) on the positive bump.

Writing \(w=0.4\), \(\eta=1.8\), \(r=14\), the density entries obey

\[
 C_{14}\le \left(\frac4\pi+8w^{-13}\right)B_{14}
 +2\eta\left(1.001w^{-13}B_{14}
 +14(10.01)w^{-12}B_{13}
 +91(25.025)w^{-11}B_{12}\right).
\]

With the certified \(B_k\) ceilings this gives \(C_{14}<2.782\times10^{33}\). The calculation uses the elementary upper bound \(\zeta(14)\le1+2^{-14}+2^{-13}/13\), avoiding unverified special-function numerics. `mpmath.iv` evaluation yields

\[
 E_{8500}(T)<9.3550840\times10^{-23}.
\]

## Result 2 — an exact ideal bath of transfer order m=2 (without rounded moment coefficients)

Let \(\vartheta(x)=b((x-0.5)/0.4)\), supported in \([0.1,0.9]\). Define the positive measure \(d\mu(x)=\vartheta(x)x^{-2}dx\) and its unique **monic orthogonal polynomial of degree two**, \(p_2(x)=x^2+a x+b\), by

\[
 \int p_2\,d\mu=0,\qquad \int xp_2\,d\mu=0.
\]

Its coefficients are defined **exactly** by these integrals, not by the rounded decimal approximations. Set \(p=p_2/\max_{x\in[0.1,0.9]}|p_2(x)|\). Then

\[
 \int\vartheta(x)\frac{p(x)}{x}\,dx=0,\qquad
 \int\vartheta(x)\frac{p(x)}{x^2}\,dx=0
\]

**exactly**. The third moment is strictly positive: orthogonal polynomial zeros lie inside \((0.1,0.9)\), so \(p_2(0)>0\), and

\[
 \int\vartheta(x)\frac{p_2(x)}{x^3}\,dx
 =\int\frac{p_2(x)}x\,d\mu(x)
 =\frac1{p_2(0)}\int\frac{p_2(x)^2}{x}\,d\mu(x)>0.
\]

The last identity uses orthogonality to the linear polynomial \((p_2(x)-p_2(0))/x\). Thus the exact smooth positive-frequency bath has transfer order **exactly two**. Its numerical monic polynomial (informational, not the exact definition) is

\[
 p_2(x)\approx x^2-0.86595156053103764976333104627649x
 +0.16164064224041778253637867509165.
\]

This corrects the historical `poly_jet` implementation's two residual lower moments (~\(10^{-15}\)). It does not retroactively transform the old N-mode simulations into results for the exactly orthogonal polynomial; rerunning them with the exact coefficients would be needed to make that identification.

## Result 3 — a state-adaptive multiplicative BKM transfer theorem

For any invertible discretized survival matrix \(W=V_N(T)\), suppose \(\|V-W\|\le E<\sigma_{\min}(W)\) and let \(\delta=E/\sigma_{\min}(W)\), \(s=\operatorname{Tr}(W C_0W^\dagger)\), \(a=1-s\), \(\kappa=(w^\dagger C_0^{-1}w)/a<1\). Put \(a_\pm=(1\pm\delta)^2\) and

\[
 c_-=\frac a{1-a_-s}<1,\qquad c_+=\frac a{1-a_+s}>1,
 \qquad \Phi(c,\kappa)=\frac{c^2(1-\kappa)}{1-c\kappa}.
\]

Assume \(1-a_+s>0\), \(c_+\kappa<1\). Let \(\Delta(W)\) be the BKM defect obtained from the exact scalar Schur kernel. Then

\[
 \boxed{a_-\Phi(c_-,\kappa)\,\Delta(W)
 \ \le\ \Delta(V)\ \le\ a_+\Phi(c_+,\kappa)\,\Delta(W).}
\]

**Proof.** Set \(R=VW^{-1}\), \(\|R-I\|\le\delta\). Then \(C_V=RC_WR^\dagger\), \(v_V=Rv_W\) and \(R^{-1}R^{-\dagger}\) lies in the Loewner interval \([(1+\delta)^{-2}I,(1-\delta)^{-2}I]\). Accordingly \(q_W(u/a_-)\le q_V(u)\le q_W(u/a_+)\) for \(q_X(u)=v_X^\dagger(C_X+uI)^{-1}v_X\). The trace satisfies \(a_-s\le\operatorname{Tr}C_V\le a_+s\). Substitute \(u=a_\pm t\) in the positive integral \(\Delta=\int_0^\infty g(q(u)/(a_X+u))du\), \(g(x)=-\log(1-x)-x\), and compare denominators. Finally, using \(g(cx)=c^2\int_0^x t/(1-ct)dt\) and \(0\le x\le\kappa\), bound the ratio \(g(cx)/g(x)\) by \(\Phi(c,\kappa)\), from below if \(c\le1\) and above if \(c\ge1\). No spectral gap or small singular-value approximation enters this bound apart from the explicit \(\delta\).


## Combining the bounds for the exactly orthogonal m=2 bath

For the exact polynomial, the separate interval and rounding estimates are
documented in `../CERTIFICATE_IDEAL_ORDER2_2026-11-22.md`. With the conservative
propagator perturbation budget `E < 1.2e-22` and the smallest-singular-value
lower bound of the displayed matrix, the transfer theorem above gives

`3.95299e-41 < Delta_BKM,continuum < 3.98513e-41`.

This conclusion concerns precisely m=2, lambda=0.05, c=5, N=8500,
(gamma_s,gamma_f)=(1,3), the initial state described in the exact-bath
certificate and the documented IEEE 754 / interval-arithmetic assumptions.
The historical SVD-rounded simulations do not inherit this interval.
