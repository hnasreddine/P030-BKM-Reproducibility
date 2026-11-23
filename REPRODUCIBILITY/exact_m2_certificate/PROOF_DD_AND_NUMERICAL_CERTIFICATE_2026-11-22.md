# P030 — Supporting numerical rounding and BKM transfer estimates

This mathematical note records the double-double operations, their error
estimates and the positive-kernel multiplicative transfer argument. Some
illustrative numerical values below concern the **rounded-SVD-polynomial
bath**, not the separate exactly orthogonal m=2 model. The exact-bath proof,
conservative total budget and final interval are in
`CERTIFICATE_IDEAL_ORDER2_2026-11-22.md`.

## 1. Double-double operations used in the C++ propagation

Let u=2^-53 be the IEEE 754 binary64 unit roundoff, ν=2^-1074 the minimum positive subnormal, and write each double-double number as X=xh+xl with |xl|≤u|xh|. All operations below are performed round-to-nearest, ties-to-even; explicit `std::fma` is correctly rounded and the C++ compiler must not apply reassociation or fast-math. No overflow or invalid nonfinite arithmetic is permitted. As usual, for a finite basic arithmetic operation, a conservative rounding model is `|fl(t)-t|≤u|t|+ν`, with the ν term covering underflow.

The actual source routines in `cheb_dd_lambda_corrected.cpp` and `strict_dd_verification.cpp` implement:

**Add:** `s=fl(ah+bh); v=fl(s-ah); t0=fl(fl(bh-v)+fl(ah-fl(s-v))); t=fl(fl(t0+al)+bl); z=fl(s+t); l=fl(t-fl(z-s))`. **Sub** is the same, with b negated. This is an exact `TwoSum` for the high components followed by two rounded corrections and a `FastTwoSum`-style output.

**Multiply:** `p=fl(ah*bh); r=fl_fma(ah,bh,-p); e=fl(fl(r+fl(ah*bl))+fl(al*bh)); z=fl(p+e); l=fl(e-fl(z-p))`. This is two-product/FMA followed by both cross products, omitting `al*bl` (a second-order term), then `FastTwoSum`-style output.

### Lemma: conservative operation bounds

For inputs satisfying |al|≤u|ah| and |bl|≤u|bh|, under the environment above, the resulting unevaluated sum Z=zh+zl obeys

`|Z-(A+B)| ≤ 100 u² (|A|+|B|) + 1e-300`;

`|Z-(A-B)| ≤ 100 u² (|A|+|B|) + 1e-300`;

`|Z-AB| ≤ 100 u² |AB| + 1e-300`.

The relative factors here use the magnitudes of the *exact double-double input values* A=ah+al, B=bh+bl, as in the instrument's `mg` upper bound. Because (1-u)|ah|≤|A|≤(1+u)|ah|, replacing |ah|+|bh| by |A|+|B| changes the constants by less than one part in 10^15; the generous value 100 absorbs this.

**Addition proof, including cancellation.** Put H=|ah|+|bh|. Standard exact TwoSum gives `s+t0=ah+bh` in the normal range; if an intermediate underflows, all departures from the identity are absorbed in the explicitly allowed absolute error `1e-300` (a finite fixed number of binary64 operations each has absolute underflow error ≤ν). Also |t0|≤uH+ν, and the normalized input assumption gives |al|+|bl|≤uH. The two rounded additions making `t=fl(fl(t0+al)+bl)` contribute at most 6u²H plus finitely many ν terms, and |t|≤3uH plus tiny absolute terms. If |s|≥|t|, the final `FastTwoSum` is an error-free transformation (outside underflow) and returns `z+l=s+t`. Otherwise |s|<|t|≤3uH, and direct rounding bounds on `z=fl(s+t)`, `w=fl(z-s)`, `l=fl(t-w)` give `|z+l-(s+t)|≤30u²H` plus finitely many ν terms: indeed |z|≤7uH, |z-s|≤10uH, |t-w|≤14uH, and the two subtractions introduce ≤(10+14)u²H. Therefore the total is less than 40u²H plus finitely many ν terms. The same argument, with exact sign reversal of bh and bl, proves subtraction. This is safely below 100u²(|A|+|B|)+1e-300. No lower bound on the *true result* |A+B| is needed: the bound is absolute and cancellation-safe.

**Multiplication proof.** Put H=|ah bh|. With an exact FMA residual, `p+r=ah bh`, |r|≤uH+ν. The products ah bl and al bh have magnitude at most uH each; their individual rounding errors are ≤u²H+ν, and the two subsequent additions each incur O(u²H) error because their operands are O(uH). A safe combined ceiling is 10u²H plus finitely many ν. Meanwhile |e|≤4uH, whereas |p|≥(1-u)H for normal products, so the final `FastTwoSum` is exact (up to underflow). The omitted product |al bl|≤u²H. Thus the total multiplication error is <12u²H plus finitely many ν, and hence <100u²|AB|+1e-300. Cases in which high product or intermediate operations underflow are covered by a fixed number of absolute ν terms; 1e-300 is many orders of magnitude larger than those terms. This proof presupposes *correctly rounded IEEE 754* arithmetic and no overflow, as checked in the runtime preflight.

The lemma follows from the standard TwoSum/FastTwoSum identities and FMA error-free product residual, **not** merely from randomized testing. It is specific to the implementation and its required normalized input condition; a complete machine-level formalization of the compiler and library semantics has not been carried out here.

## 2. Arithmetic preconditions and recorded residuals

- The verification executable `strict_dd_verification.cpp` was compiled with `g++ -O3 -std=c++17 -fno-fast-math -ffp-contract=off` and checked `std::numeric_limits<double>::is_iec559` and `fegetround()==FE_TONEAREST` before execution: **PASS**.
- It processed both columns with the recorded 8500 modes and K=12865. The recorded output was reproduced (`cmp` exit status 0).
- The strict verifier recorded `BAD_NONOVERLAP 0`, `MAX_OVERLAP 1.1102056881859389490e-16 < u=1.1102230246251565e-16`. The check is on the actual intermediate DD results in the propagation, not a claim about every possible input to the program.
- An additional exact-rational comparison of all 46,867 parsed decimal tokens against binary64 `high+low` found zero violations of `|DD(decimal)-decimal| ≤ 2^-100 |decimal| + 1e-300`. The maximal observed overlap among nonzero parsed DD inputs was 0.997592u. Three extremely tiny input tokens were below 1e-300 and were represented as exact binary zero, which is permitted by the absolute term.
- An independent, deterministic stress script checked 25,000 each of +, − and × against **exact rational** results, with error ratios to u² scale below 1.83, 1.87 and 2.96 respectively. This is a counterexample search, not the proof.
- The instrument reports `GLOBAL_SUMRES 7.8585094074430482021e-26`, `INITRES 3.0299503750250085138e-31`, `OUTSUM 2.2659878295377634434e-29`. The strict replay reproduces these values.

Caveat: the numerical overlap maxima use normal binary64 comparisons; the gap below u is much larger than a rounding ulp, so this check is decisive for the actual run, but independent formal compiler/runtime verification remains a separate level of assurance.

## 3. Chebyshev propagation and Bessel input error

The input certificate gives `||G||_F ≤ 2.323522` and `||H_N/2|| ≤ 0.558089 < 0.6`. For selfadjoint A=H_N/2, the second-kind polynomials obey the *uniform* spectral estimate

`||U_k(A)|| ≤ 1/sqrt(1-||A||²) < 1.25`, for all k≥0.

Thus the recurrence errors in the two-column instrument are propagated without a factor K. The coefficient ℓ¹ bound is `|J0|+2 sum_{k>=1}|Jk| < 134.051`. With a 0.5% safety margin for strictly nonnegative floating error-radius accumulations (much larger than γ_n for the source's <10^12 operations), one gets

`E_DD < 1.324e-23`.

The exact-rational parser verification justifies the missing term `E_Bessel_parse ≤ 134.051 × 2^-100 + (2K+1)×1e-300 < 1.058e-28` from converting the recorded decimal Bessel coefficients to DD, separately from the already validated ODE residual coefficient error <1.717e-32. The conversion of source decimal grid/coupling tokens into DD and the resulting Hamiltonian perturbation has budget <1.236e-24 (from the existing source bound). Conversion of the binary64 DD output to the printed decimal representation is <1.243e-24, checked against exact fractions. The series truncation estimate is <3e-97. The recorded time token is within 4.605e-49 of 2000 log(20), so a generous contribution of 1e-47 is harmless.

The continuum aliasing certificate (with the corrected Poisson shift phase (-1)^k) gives `<9.356e-23` between the continuous Hamiltonian and **exact** N-mode propagation. Summing the conservative budgets:

| Component | Conservative bound |
|--|--:|
| Continuum → exact 8500-mode model | 9.356e-23 |
| DD recurrence and accumulation | 1.324e-23 |
| Printed output conversion | 1.243e-24 |
| Decimal input → DD Hamiltonian | 1.236e-24 |
| Decimal Bessel coefficients → DD | 1.058e-28 |
| Bessel residual/ODE, Chebyshev tail, time | <2e-32 (plus negligible 1e-47) |
| **Total, with conservative round-ups** | **<1.093e-22** |

The conservative bound `E=1.2e-22` has a margin of approximately 1.07e-23 in this numerical environment. These arithmetic estimates require the stated IEEE/compiler environment; arbitrary-platform portability is not asserted.

## 4. Positive BKM integral and transfer

The prior Darboux interval script treats the *printed decimal* W_hat as an exact matrix and evaluates the positive integral

`Delta(W_hat) = ∫_0^∞ g(x(t)) dt`, `g(x)=-log(1-x)-x`,
`x(t) = (alpha t+ kappa0 d)/((t²+s t+d)(a+t))`, where `s=tr(C), a=1-s, d=det(C), alpha=0.01 ||W_hat e2||²`, `kappa0=.05`.

A missing explicit monotonicity check in the earlier narrative can be supplied: with `D(t)=(t²+s t+d)(t+a)`, the sign of x'(t) is that of

`-2 alpha t³ - [alpha(s+a)+3 kappa0 d] t² - 2 kappa0 d(s+a) t + d[alpha a - kappa0(d+s a)]`.

Here interval-evaluated alpha≈1.7314e-27, s≈3.89034e-14, a≈1, d≈3.56115e-53, kappa0=0.05; hence `alpha a < kappa0(d+s a)` with a margin ≈1.945e-15. Every coefficient is therefore strictly negative for t≥0. The Darboux monotonicity assumption is valid **for this saved matrix**.

### Proof of the multiplicative BKM transfer

Let `C=W C0 W*`, `b=W b0`, and `R=V W^{-1}` (well-defined because W is invertible). If `||V-W|| ≤ E`, with `delta=E/sigma_min(W)<1`, then `||R-I||≤delta`, hence `a_minus I ≤ R*R ≤ a_plus I` with `a_minus=(1-delta)^2` and `a_plus=(1+delta)^2`. The perturbed quantities are `C_V=R C R*`, `b_V=R b`, `s_V=tr(C_V)` and `a_V=1-s_V`. Thus `a_minus s ≤ s_V ≤ a_plus s`, and, for `q_W(t)=b*(C+tI)^{-1}b`, the exact change-of-basis identity gives

`q_V(t)= b* (C + t(R*R)^-1)^-1 b`.

Inverting the two Loewner inequalities yields `q_W(t/a_minus) ≤ q_V(t) ≤ q_W(t/a_plus)`. Set `x_W(t)=q_W(t)/(a+t)` and `x_V(t)=q_V(t)/(a_V+t)`. Using the bounds on a_V and the elementary facts `a_minus≤1≤a_plus`, one gets, for every t≥0,

`x_V(a_minus*t) ≥ c_minus*x_W(t)`, `x_V(a_plus*t) ≤ c_plus*x_W(t)`,

where `c_minus = a/(1-a_minus*s)≤1`, `c_plus=a/(1-a_plus*s)≥1`. In particular positivity of `1-a_plus*s` is necessary and is verified here.

For `g(x)=sum_{k=2}^∞ x^k/k`, let `f(x)=(1-x)g(x)/x²` for 0<x<1. Expanding, `f(x)=1/2+sum_{n≥1}(1/(n+2)-1/(n+1))*x^n`, so f is nonincreasing. Hence, for 0<c<1, `g(cx)/g(x)≥c²(1-x)/(1-cx)`; for c>1 with cx<1, the inequality reverses. Since `0≤x_W(t)≤kappa:=b*C^-1 b/a`, these imply

`g(c_minus*x_W)≥Phi(c_minus,kappa)*g(x_W)`,

`g(c_plus*x_W)≤Phi(c_plus,kappa)*g(x_W)`,

where `Phi(c,kappa)=c²(1-kappa)/(1-c*kappa)`. Changing variables `t=a_minus v` and `t=a_plus v` in the positive integral gives the fully explicit inequality

`a_minus Phi(c_minus,kappa) Delta(W) ≤ Delta(V) ≤ a_plus Phi(c_plus,kappa) Delta(W)`.

This proves the multiplicative transfer step directly, provided `C0>0`, W invertible, δ<1, `1-a_plus*s>0`, and `c_plus*kappa<1`, all of which are satisfied for the recorded parameter intervals. It does **not** depend on an empirical observation that the relative BKM error is small.

The 122,099-bin Darboux evaluation gives

`3.9670526733148689e-41 < Delta(W_hat) < 3.9710197759486385e-41`.

The interval singular-value bound `sigma_min(W_hat)>6.7652934423519234e-20`, together with `E=1.2e-22`, yields delta <0.001774. Under the multiplicative BKM transfer theorem independently proved in the preceding paragraph, and its positivity requirements, the final interval is

**`3.9529919639270994e-41 < Delta_BKM,continuum < 3.985119533633871e-41`.**

The rounding-operation and BKM-transfer arguments apply under the stated numerical assumptions. For the exactly orthogonal bath and the full continuum budget, use `CERTIFICATE_IDEAL_ORDER2_2026-11-22.md` and the executable `final_error_budget.py`, rather than the rounded-SVD-polynomial values in this auxiliary note.


## References for floating-point arithmetic identities

- M. Joldes, J.-M. Muller, and V. Popescu, *Tight and Rigorous Error Bounds for Basic Building Blocks of Double-Word Arithmetic*, ACM Transactions on Mathematical Software 44 (2017).
- S. Boldo, S. Graillat, and J.-M. Muller, *On the Robustness of the 2Sum and Fast2Sum Algorithms*, ACM Transactions on Mathematical Software 44 (2017).
