#!/usr/bin/env bash
# Invoke only AFTER activating this repository in Zenodo's GitHub integration.
set -euo pipefail
OWNER="${P030_GITHUB_OWNER:-hnasreddine}"
REPO="${P030_GITHUB_REPO:-P030-BKM-Reproducibility}"
command -v gh >/dev/null || { echo 'Install/authenticate gh first.'; exit 2; }
gh repo view "$OWNER/$REPO" >/dev/null
read -r -p "Confirm Zenodo GitHub integration is enabled for $OWNER/$REPO. Type ZENODO: " response
[[ "$response" == "ZENODO" ]] || { echo CANCELLED; exit 1; }
if gh release view v1.0.0 -R "$OWNER/$REPO" >/dev/null 2>&1; then
  echo 'STOP: Release v1.0.0 already exists; no duplicate.' >&2;exit 3
fi
gh release create v1.0.0 -R "$OWNER/$REPO" --target main --title 'P030 BKM reproducibility v1.0.0' --notes 'Code, input data, recorded outputs and computer-assisted m=2 certificate; detailed numerical scope and mathematical assumptions in README.md and REPRODUCIBILITY/exact_m2_certificate/README_REPRODUCE.md.'
printf 'Release created. Check Zenodo integration for the actual DOI, then update JMP Data Availability.\n'
