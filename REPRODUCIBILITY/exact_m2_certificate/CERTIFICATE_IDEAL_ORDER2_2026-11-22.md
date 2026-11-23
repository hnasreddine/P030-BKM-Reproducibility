# P030 — Exactly orthogonal m=2 bath: continuum BKM enclosure

**Date:** 22 November 2026. **Model:** vacuum, positive-frequency rotating-wave one-excitation Friedrichs model at lambda=0.05, c=5, N=8500, gamma_s=1, gamma_f=3, T=2000 log(20). **Scope:** one specified bath and one specified time, **not** a universal numerical theorem for all lambda or all transfer orders.

**Result:** The numerical enclosure below is supported by the accompanying analytical estimates and executable interval computations for the single specified model and arithmetic environment. It is a computer-assisted result, not a proof-assistant-checked theorem or a statement about other parameter values.

## A. Exact order-two polynomial, not SVD rounding

Let `b(t)=exp(1-1/(1-t^2))` for `|t|<1`, zero otherwise, `theta(x)=b((x-.5)/.4)`. Define

`I_j = integral_{.1}^{.9} theta(x) x^j dx`, `j=-2,-1,0,1`.

Let `p_2(x)=x^2+a x+b` be the unique monic polynomial satisfying

`b I_{-2} + a I_{-1}+I_0 =0`, `b I_{-1}+a I_0+I_1=0`.

Set `p(x)=p_2(x)/max_{x in [.1,.9]} |p_2(x)|`. The maximum is attained at x=.9; this was proven with interval arithmetic by comparing both endpoints and the vertex. For the *mathematically exact* p:

`integral theta(x) p(x)/x dx = integral theta(x) p(x)/x^2 dx = 0`.

The next weighted moment is strictly positive. Indeed, with `dmu=theta(x)x^{-2}dx`, use orthogonality of `p_2` to the linear polynomial `(p_2(x)-p_2(0))/x` to get

`integral theta p_2(x)/x^3 dx = (1/p_2(0)) integral theta p_2(x)^2 / x^3 dx > 0`, since `p_2(0)>0`.

Thus the **ideal positive bath genuinely has transfer order exactly m=2** (the rounded-SVD bath does not enforce exact zero moments).

### Validated coefficient intervals

`certify_polynomial_moments.py` encloses the four integrals with `mpmath.iv` and **150,000** uniformly spaced panels. The analytic Poisson/trapezoid remainder uses `r=14`, the established Sturm bounds `B12<4.378e21`, `B13<2.412e24`, `B14<1.566e27`, and the elementary inequality `B_l<=2^(12-l) B12` for `l<=12`. For each integrand `theta(x) x^j`, the 14th derivative L1 norm is bounded by the product rule and elementary bounds on derivatives of x^j for x>=.1. The remaining quadrature error is bounded by `2*zeta(14)*C14*(h/(2*pi))^14`; rational conservative replacements `2*zeta(14)<2.002` and `2*pi>6.28` are used. All quadrature sums and coefficient algebra are interval-enclosed.

This proves that the sum of absolute coefficient differences between the exactly orthogonal polynomial and the recorded high-precision decimal polynomial is **less than 4.703e-30**, and in particular **less than 5e-30**, uniformly implying `sup_{x in [-1,1]} |p_exact(x)-p_decimal(x)|<5e-30`.

The recorded decimal polynomial coefficients, low-to-high, are:

```
0.840633866412294604753274979989083894100408184544464751...
-4.503497377672166729567865731729741282194883749015500736...
5.200634288262537594886177998231707728240724925394427051...
```

See `moments_certificate.log` for the full coefficient interval endpoints.

## B. Analytic comparison of the ideal and decimal-coefficient Hamiltonians

The spectral density has the form

`J_p(x) = diag(gamma_s phi(x)/pi, gamma_f phi(x)/pi) + 2 theta(|x|) I + 1.8 theta(x) p(x) [[0,1],[1,0]]`.

On the support of `theta`, both the exact and recorded polynomial have absolute value at most `1+5e-30`. Hence `lambda_min(J_p(x)) >= .19 theta(x)`, and the same for `J_q(x)`. Where `theta=0`, the two matrices coincide. From the Sylvester equation for the difference of positive matrix square roots,

`||sqrt(J_p)-sqrt(J_q)||_F <= ||J_p-J_q||_F / (sqrt(lambda_min(J_p))+sqrt(lambda_min(J_q)))`.

Since `||J_p-J_q||_F <= sqrt(2)*1.8 theta |p-q|`, and the discrete coupling is `G_j=sqrt(h)*sqrt(J(x_j))`, the grid sum `h sum_j theta(x_j) < .801` yields

`||G_p-G_q||_F <= (1.8 sqrt(2)/(2 sqrt(.19))) sqrt(.801) * (5e-30) < 1.307e-29`.

The Hermitian arrowhead Hamiltonian satisfies `||H_p-H_q||_2 <= lambda ||G_p-G_q||_F < 6.534e-31`. Duhamel's formula for the two unitary groups gives

**`||V_{N,p}(T)-V_{N,q}(T)||_2 < 4e-27`.**

This closes the *exact-moment polynomial versus rounded numerical representation* error at the specified lambda, c and N, without assuming a tiny BKM difference merely from close coefficients. See `sensitivity_to_exact_polynomial.py`.

## C. Direct Hamiltonian propagation and interval BKM integral

`prepare_exact_bath.py` generated all 8500 high-precision grid nodes and 2x2 matrix-square-root coupling blocks using the coefficient interval specified in Section A and the Bessel coefficients and the same lambda/T. The corrected full-Hamiltonian C++ Chebyshev propagation was run with `-O3 -std=c++17 -fno-fast-math -ffp-contract=off`.

A separately instrumented strict double-double run reached K=12865 on both columns and reproduced the non-instrumented output **byte-for-byte**. The DD input-parser verification tested all 46,867 tokens against exact rationals: zero violations of the specified parser budget. No disallowed DD nonoverlap was detected. The instrument reports:

- `GLOBAL_SUMRES 7.8585094074430482021e-26`
- `INITRES 3.0299503750250085138e-31`
- `OUTSUM 2.2659878295377634434e-29`
- `BAD_NONOVERLAP 0`, `MAX_OVERLAP 1.1102095020735415134e-16` (below `u=2^-53`)

The positivity integral for the exactly recorded decimal output matrix `W_hat` is

`Delta(W_hat) = integral_0^infty g(x(t)) dt`, with `g(x)=-log(1-x)-x` and

`x(t)=(alpha*t+kappa*d)/((t^2+s*t+d)*(a+t))`,

`C=0.2 W_hat W_hat^*`, `s=tr C`, `a=1-s`, `d=det C`, `alpha=0.01 ||W_hat e2||^2`, `kappa=.05`.

`check_BKM_monotonicity.py` proves with intervals that **all coefficients in the derivative numerator of x(t) are negative**, so `g(x(t))` is decreasing. `certified_BKM_tight.py` applied rigorous interval Darboux sums over 122,099 bins, with explicit tails, obtaining:

**`3.9670526841476009e-41 < Delta(W_hat) < 3.9710197867922034e-41`.**

This enclosure concerns the output matrix as exact displayed decimals, not the continuum dynamics before the error budget.

As an independent **numerical diagnostic**, 160-digit eigenvalue evaluation of the full 3x3 state and the dephased 2x2 excited block gives `Delta(W_hat)=3.96903554965836265433e-41`, inside the positive-kernel Darboux interval (`independent_BKM_diagonalization.py`). This ordinary multiprecision calculation is a cross-check, **not** an independent interval proof.

## D. Full operator error budget and continuum transfer

The continuum-to-exact-discrete bound applies to the exact polynomial as well as the rounded-SVD polynomial. The program `full_certified_operator_bound.py` gives `9.355083991611e-23` (the Poisson alias identity for the shifted grid includes `(-1)^k`; its absolute-value bound is unchanged).

The double-double arithmetic lemma proved in the accompanying note is used **under** the actual IEEE 754 binary64 nearest-rounding mode, normalized operands, correct FMA and no reassociation. The strict C++ execution recorded normalized double-double values and identical outputs. Conservative accounting of all components gives:

| Error source | Strict upper budget |
|---|---:|
| Continuum to exact 8500-mode Hamiltonian | `9.356e-23` |
| Arithmetic recurrence, including local DD residuals | `1.334e-23` |
| Printed decimal output | `1.3e-24` |
| Decimal input to double-double Hamiltonian | `1.3e-24` |
| Ideal polynomial vs recorded polynomial | `4e-27` |
| Bessel coefficient input conversion | `1.07e-28` |
| Bessel ODE residual | `1.717e-32` |
| Time and ideal Chebyshev truncation | negligible at this precision |
| **TOTAL** | **`<1.09505e-22 <1.2e-22`** |

The coefficient bound `|J_0|+2 sum_{k=1}^K |J_k|<135` was checked with exact decimal arithmetic; the binary64 printed-output discrepancy was checked with exact rational fractions. The safety margin on double-double residual accumulation is 0.5%. See `final_error_budget.py`.

The matrix `W_hat` has interval-lower-bounded smallest singular value `>6.76529344235146e-20`, allowing the already established *multiplicative BKM transfer theorem*. Using the deliberately loose bound `E=1.2e-22`, `delta < .001774`, and the interval Darboux bounds gives the **continuum enclosure for the exactly orthogonal m=2 bath**:

**`3.9529919747214352e-41 < Delta_BKM,continuum < 3.9851195445159386e-41`.**

One may conservatively quote the outward-rounded interval

**`3.95299e-41 < Delta_BKM,continuum < 3.98513e-41`**.

No small-defect subtraction was used. The transfer employs positivity of the Schur kernel and Loewner comparison of `R^* R` for `R=V W_hat^-1`.

## E. Exact scope and computational requirements

This enclosure is restricted to the exactly orthogonal order-two bath at the stated
parameter point, the stated initial state, the specified discretization and the
IEEE 754 non-contracted binary64 execution. The rounded-SVD polynomial
and its other numerical table entries are distinct experiments. The analytic
proofs and programs in this archive give the mathematical assumptions and
individually checkable bounds; they do not constitute proof-assistant verification
or a theorem for all parameter values.

The aggregate budget script also includes the 42-digit printed-grid/coupling
contribution, bounded by `2.2e-37` for the propagator, which is negligible within
the conservative `<1.2e-22` total.
