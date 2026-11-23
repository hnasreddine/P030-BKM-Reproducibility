#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
python tools/verify_release.py
cd REPRODUCIBILITY/exact_m2_certificate
python final_error_budget.py
python transfer_ideal_interval.py
python independent_BKM_diagonalization.py
python check_orthogonal.py
printf '\nP030 precomputed verification: PASS (recorded evidence and one-point mathematical hypotheses)\n'
