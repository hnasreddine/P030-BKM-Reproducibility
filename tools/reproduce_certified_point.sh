#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CERT="$ROOT/REPRODUCIBILITY/exact_m2_certificate"
python "$ROOT/tools/verify_release.py"
cd "$CERT"
python final_error_budget.py
python transfer_ideal_interval.py
python independent_BKM_diagonalization.py
CXX="${CXX:-g++}"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
"$CXX" -O3 -std=c++17 -fno-fast-math -ffp-contract=off cheb_dd_lambda_corrected.cpp -o "$TMP/cheb_dd" -lquadmath
"$TMP/cheb_dd" input_005_5_8500_ideal_m2.txt "$TMP/out.txt"
cmp "$TMP/out.txt" out_ideal_m2_lean.txt
printf 'PASS: regenerated ideal-bath C++ output matches recorded values byte-for-byte\n'
if [[ "${1:-}" == "full" ]]; then
 python certify_polynomial_moments.py
 python certified_bessel_ode.py
 python certified_BKM_tight.py out_ideal_m2_lean.txt
 (cd continuum_proof && python full_certified_operator_bound.py)
fi
printf 'P030 exact m=2 single-point reproduction: PASS (scope H1-H5)\n'
