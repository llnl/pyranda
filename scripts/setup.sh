#!/usr/bin/env bash
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

python3 -m venv .venv
.venv/bin/pip3 install -r scripts/requirements.txt

echo "Setup complete."
