#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
bash RUN_PRECOMPUTED.sh
bash tools/reproduce_certified_point.sh full
# Recompute the instrumented floating-point rounding bounds from the canonical input.
CERT="$ROOT/REPRODUCIBILITY/exact_m2_certificate"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cd "$CERT"
"${CXX:-g++}" -O3 -std=c++17 -fno-fast-math -ffp-contract=off strict_dd_audit.cpp -o "$TMP/strict_dd" -lquadmath
"$TMP/strict_dd" input_005_5_8500_ideal_m2.txt "$TMP/strict_out.txt" 2> "$TMP/strict.log"
cmp "$TMP/strict_out.txt" out_ideal_m2_strict.txt
grep -q 'BAD_NONOVERLAP 0' "$TMP/strict.log"
grep -q 'NONFINITE 0 ZERO_HIGH_NONZERO_LOW 0' "$TMP/strict.log"
cat "$TMP/strict.log" | tail -3
printf '\nP030 full exact-bath computational replay: PASS (within H1-H5 assumptions)\n'
