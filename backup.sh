#!/usr/bin/env bash
set -euo pipefail

# Automatically snapshot and push the entire repo.
# Assumes the remote "origin" and branch are already set up and credentials cached.

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_DIR}"

git add -A
if git diff --cached --quiet; then
  echo "Nothing to commit."
  exit 0
fi

STAMP=$(date +"%Y-%m-%d %H:%M:%S")
git commit -m "Automated backup ${STAMP}"
git push origin HEAD
