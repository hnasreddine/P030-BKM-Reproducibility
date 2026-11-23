# Mathematical proof and executable evidence map (P030, release v1.0.0)

**Exact scope**: the continuum BKM enclosure `3.95299e-41 < Δ < 3.98513e-41`
for the exactly orthogonal order-two spectral density with λ=0.05, c=5,
N=8500, (γs,γf)=(1,3), T=2000 log 20 and the explicitly specified
initial state. This is a **computer-assisted** one-point statement; the
historical rounded-SVD table entries do not inherit this enclosure.

| Mathematical obligation | Mathematical source | Reproduction script/data |
| --- | --- | --- |
| Moment orthogonality, order two, positivity | `REPRODUCIBILITY/exact_m2_certificate/CERTIFICATE_IDEAL_ORDER2_2026-11-22.md` | `certify_polynomial_moments.py`, `check_orthogonal.py`, `orthogonal_coefs.txt` |
| Exact polynomial vs recorded coefficients | same certificate | `sensitivity_to_exact_polynomial.py`, `validate_ideal_inputs.py` |
| Continuum vs N=8500 Fourier-grid/Volterra bound | `continuum_proof/CONTINUUM_BOUND_AND_TRANSFER.md` | `continuum_proof/full_certified_operator_bound.py`, certified roots and derivative bounds |
| Chebyshev/Bessel, parser, double-double error | `PROOF_DD_AND_NUMERICAL_CERTIFICATE_2026-11-22.md` | `cheb_dd_lambda_corrected.cpp`, `strict_dd_audit.cpp`, `certified_bessel_ode.py`, inputs, outputs and `final_error_budget.py` |
| Positive-kernel BKM enclosure and transfer | `CERTIFICATE_IDEAL_ORDER2_2026-11-22.md` | `certified_BKM_tight.py`, `transfer_ideal_interval.py`, `independent_BKM_diagonalization.py` |

The proof is not mechanically checked (no Lean/Coq/Isabelle); the stated
trusted base includes IEEE 754 round-to-nearest, correctly rounded FMA,
non-contracted arithmetic and the implementation of interval/exact libraries.
A completed numerical replay is not a proof of the underlying inequalities;
read the derivations as well as running the scripts.

**Fast recorded-value and budget checks**: `bash RUN_PRECOMPUTED.sh`.
**Full numerical reconstruction**: `bash RUN_EXACT.sh`.
They are separate routes to reconstruct the calculation; they are not two
independent formal proofs. For figures and uncertified historical simulations,
see `REPRODUCIBILITY/code/`, `REPRODUCIBILITY/data/` and `order2_high_precision/`.
