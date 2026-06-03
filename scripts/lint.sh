#!/usr/bin/env bash
set -euo pipefail

# Resolve an optional target path BEFORE we cd, since it's relative to the
# caller's current directory, not the repo root.
target="."
if [ $# -gt 0 ]; then
  target="$1"
  case "$target" in
    /*) ;;                      # already absolute
    *) target="$(pwd)/$target" ;;
  esac
fi

# Run from the repo root so ruff finds its config.
cd "$(git rev-parse --show-toplevel)"

if [ "$target" != "." ] && [ ! -e "$target" ]; then
  echo "error: no such file or directory: $1" >&2
  exit 1
fi

# Locate ruff inside the project venv (handles Linux/macOS and Windows git-bash).
if [ -x ".venv/bin/ruff" ]; then
  RUFF=".venv/bin/ruff"
elif [ -x ".venv/Scripts/ruff.exe" ]; then
  RUFF=".venv/Scripts/ruff.exe"
else
  echo "error: ruff not found in .venv" >&2
  echo "  Set up the environment first:" >&2
  echo "    python -m venv .venv" >&2
  echo "    .venv/bin/pip install -r requirements.txt" >&2
  exit 1
fi

if [ "$target" = "." ]; then
  echo "Running ruff on all Python files in the repo..."
else
  echo "Running ruff on: $target"
fi

# Either command failing exits non-zero and blocks the commit.
"$RUFF" check "$target"          # lint
"$RUFF" format --check "$target" # verify formatting

echo "ruff: all checks passed."
