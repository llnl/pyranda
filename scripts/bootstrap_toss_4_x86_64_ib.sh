#!/bin/bash

set -euo pipefail

test -e INSTALL.md || { echo "this must run from the root of pyranda"; exit 1; }

env_dir="${1:-myEnv}"
python_bin="${PYTHON_FOR_VENV:-$(command -v python3 || true)}"
# The toss_4_x86_64 Python stack already carries numpy/scipy/matplotlib.
# Default to inheriting it so bootstrap does not depend on external index access.
venv_args_string="${BOOTSTRAP_VENV_ARGS:---system-site-packages}"
venv_args=()

if [[ -n "$venv_args_string" ]]; then
  read -r -a venv_args <<< "$venv_args_string"
fi

test -n "$python_bin" || {
  echo "unable to find python3 for virtual environment creation"
  exit 1
}

if [[ -e "$env_dir" ]]; then
  echo "environment path already exists: $env_dir"
  echo "remove it first or choose a different path"
  exit 1
fi

echo "creating virtual environment: $env_dir"
"$python_bin" -m venv "${venv_args[@]}" "$env_dir"

echo "bootstrapping pyranda into: $env_dir"
bash scripts/toss_4_x86_64_ib_cray.sh "$env_dir"

cat <<EOF

bootstrap complete
activate with:
  source $env_dir/bin/activate
EOF
