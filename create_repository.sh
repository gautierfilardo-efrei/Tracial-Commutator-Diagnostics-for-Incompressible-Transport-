#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_owner="gautierfilardo-efrei"
repo_name="tracial-transport-diagnostics"

command -v gh >/dev/null || { echo "GitHub CLI is required (gh)." >&2; exit 1; }
command -v git >/dev/null || { echo "Git is required." >&2; exit 1; }
active_login="$(gh api --hostname github.com user --jq .login)"
if [[ "$active_login" != "$repo_owner" ]]; then
  echo "Authenticate GitHub CLI as gautierfilardo-efrei before running this script." >&2
  exit 1
fi
git -C "$project_dir" config user.name >/dev/null || {
  echo "Configure your Git user.name first." >&2; exit 1;
}
git -C "$project_dir" config user.email >/dev/null || {
  echo "Configure your Git user.email first." >&2; exit 1;
}
if git -C "$project_dir" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "This directory is already inside a Git repository. No changes made." >&2
  exit 1
fi
if gh repo view "$repo_owner/$repo_name" >/dev/null 2>&1; then
  echo "The remote repository already exists. Inspect it before integration." >&2
  exit 1
fi

git -C "$project_dir" init -b main
git -C "$project_dir" add -- .
git -C "$project_dir" commit -m "Add corrected transport diagnostics and reproducible results"
gh repo create "$repo_owner/$repo_name" --private \
  --description "Tracial transport diagnostics: exact identities, Fourier saturation and reproducible numerical checks" \
  --source "$project_dir" --remote origin --push
gh repo view "$repo_owner/$repo_name" --json nameWithOwner,url,isPrivate
