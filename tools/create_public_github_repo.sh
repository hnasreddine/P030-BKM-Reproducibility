#!/usr/bin/env bash
# Run only with user-authorized, authenticated GitHub CLI. No Zenodo submission here.
set -euo pipefail
OWNER="${P030_GITHUB_OWNER:-hnasreddine}"
REPO="${P030_GITHUB_REPO:-P030-BKM-Reproducibility}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
python "$ROOT/tools/verify_release.py"
command -v gh >/dev/null || { echo 'STOP: Install/authenticate GitHub CLI (gh auth login) first.' >&2; exit 2; }
gh auth status >/dev/null
if gh repo view "$OWNER/$REPO" >/dev/null 2>&1; then
  echo "STOP: $OWNER/$REPO already exists; refusing to overwrite" >&2; exit 3
fi
read -r -p "Create public GitHub repository $OWNER/$REPO with code and data? Type PUBLISH: " response
[[ "$response" == "PUBLISH" ]] || { echo CANCELLED; exit 1; }
TEMP="$(mktemp -d)"
trap 'rm -rf "$TEMP"' EXIT
cp -a "$ROOT"/. "$TEMP"/
cd "$TEMP"
git init -b main
git add -A
git -c user.name="${GIT_AUTHOR_NAME:-Hassan Nasreddine}" -c user.email="${GIT_AUTHOR_EMAIL:-noreply@users.noreply.github.com}" commit -m 'P030 BKM reproducibility archive v1.0.0'
gh repo create "$OWNER/$REPO" --public --source=. --remote=origin --push
printf '\nGitHub public repository created: https://github.com/%s/%s\n' "$OWNER" "$REPO"
printf 'Next: enable this repo in Zenodo GitHub integrations BEFORE creating the release tag.\n'
