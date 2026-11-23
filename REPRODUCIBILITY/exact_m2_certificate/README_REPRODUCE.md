# Reproducing the exact-m=2 single-point computer-assisted enclosure

Scope: m=2, lambda=0.05, c=5, N=8500, gamma_s=1, gamma_f=3,
T=2000 ln 20, C0=0.2 I, b=(0,0.1), a0=0.6. No certification is asserted for
any other parameter point.

Requirements: Linux x86-64, IEEE 754 binary64 round-to-nearest-even, `g++`
C++17 with `libquadmath`, no floating-point contraction or reassociation,
no FTZ/DAZ and correct FMA. Python 3.13, mpmath 1.3.0, SymPy 1.14.0, NumPy and SciPy.
The supplied interval routines use `mpmath.iv` and SymPy; no other
interval-arithmetic package is required for the included workflows.

Run at release root:

```sh
python tools/verify_release.py
bash tools/reproduce_certified_point.sh
bash tools/reproduce_certified_point.sh full
```

For the underlying derivations and exact assumptions read
`CERTIFICATE_IDEAL_ORDER2_2026-11-22.md`,
`PROOF_DD_AND_NUMERICAL_CERTIFICATE_2026-11-22.md` and
`continuum_proof/CONTINUUM_BOUND_AND_TRANSFER.md`.

Only `input_005_5_8500_ideal_m2.txt` defines the exactly orthogonal bath
used in the certified enclosure. The rounded-SVD `input_005_5_8500.txt` is
included solely to reproduce the separately labeled SVD-rounded calculations.
The full C++ rerun requires the documented non-contracted compiler flags.
