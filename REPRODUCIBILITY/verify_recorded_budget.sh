#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/exact_m2_certificate"
python3 final_error_budget.py
python3 transfer_ideal_interval.py
printf '\nRecorded-budget checks: PASS.\n'
