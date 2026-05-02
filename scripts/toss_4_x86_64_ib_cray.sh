#!/bin/bash

set -euo pipefail

test -e INSTALL.md || { echo "this must run from the root of pyranda"; exit 1; }

version="${MPI4PY_VERSION:-4.1.0}"
name="mpi4py-$version"
tarball="${name}.tar.gz"
env_dir="${1:-myEnv}"
python_bin="${PYTHON_BIN:-$env_dir/bin/python}"
mpi4py_mode="${PYRANDA_MPI4PY_MODE:-build}"

test -x "$python_bin" || {
  echo "python interpreter not found: $python_bin"
  echo "usage: $0 [path-to-venv]"
  exit 1
}

case "$mpi4py_mode" in
  build|reuse)
    ;;
  *)
    echo "invalid PYRANDA_MPI4PY_MODE: $mpi4py_mode"
    echo "expected one of: build, reuse"
    exit 1
    ;;
esac

prepend_cray_runtime_path() {
  [[ -n "${CRAY_LD_LIBRARY_PATH:-}" ]] || return 0

  case ":${LD_LIBRARY_PATH:-}:" in
    *":${CRAY_LD_LIBRARY_PATH}:"*)
      ;;
    *)
      export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
      ;;
  esac
}

install_activate_hook() {
  local activate_file="$env_dir/bin/activate"
  [[ -f "$activate_file" ]] || return 0
  grep -q "PYRANDA_CRAY_RUNTIME_HOOK" "$activate_file" && return 0

  cat >> "$activate_file" <<'EOF'

# PYRANDA_CRAY_RUNTIME_HOOK
if [ -n "${CRAY_LD_LIBRARY_PATH:-}" ]; then
    case ":${LD_LIBRARY_PATH:-}:" in
        *":${CRAY_LD_LIBRARY_PATH}:"*)
            ;;
        *)
            export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
            ;;
    esac
fi
EOF
}

prepend_cray_runtime_path
install_activate_hook

resolve_compiler() {
  local requested="${1:-}"
  if [[ -n "$requested" ]]; then
    echo "$requested"
    return 0
  fi
  command -v "$2" 2>/dev/null || true
}

pick_mpif90() {
  local mpicc_path="$1"
  local sibling_dir
  sibling_dir="$(dirname "$mpicc_path")"

  for candidate in \
    "${MPIF90:-}" \
    "${MPIFORT:-}" \
    "$sibling_dir/mpif90" \
    "$sibling_dir/mpifort" \
    "$(command -v mpif90 2>/dev/null || true)" \
    "$(command -v mpifort 2>/dev/null || true)"; do
    [[ -n "$candidate" && -x "$candidate" ]] && {
      echo "$candidate"
      return 0
    }
  done

  echo "unable to find mpif90/mpifort compatible with $mpicc_path" >&2
  exit 1
}

find_cython_site_packages() {
  "$python_bin" -c 'import pathlib; import Cython; print(pathlib.Path(Cython.__file__).resolve().parent.parent)' 2>/dev/null && return 0
  python3 -c 'import pathlib; import Cython; print(pathlib.Path(Cython.__file__).resolve().parent.parent)' 2>/dev/null && return 0
  return 1
}

find_build_support_site_packages() {
  "$python_bin" -c 'import pathlib; import setuptools; print(pathlib.Path(setuptools.__file__).resolve().parent.parent)' 2>/dev/null && return 0
  python3 -c 'import pathlib; import setuptools; print(pathlib.Path(setuptools.__file__).resolve().parent.parent)' 2>/dev/null && return 0
  return 1
}

has_importable_mpi4py() {
  "$python_bin" -c 'import mpi4py; print(mpi4py.__version__)' >/dev/null 2>&1
}

describe_mpi4py_source() {
  if [[ -f "$tarball" ]]; then
    echo "build from local tarball $tarball"
    return 0
  fi

  if [[ -d "$name" ]]; then
    echo "build from local source tree $name"
    return 0
  fi

  echo "download $tarball from GitHub and build from source"
}

show_mpi4py_details() {
  "$python_bin" - <<'EOF'
import mpi4py
from mpi4py import MPI

print("mpi4py", mpi4py.__version__)
print("mpi4py module", mpi4py.__file__)
print("mpi4py config", mpi4py.get_config())
print(MPI.Get_library_version().splitlines()[0])
EOF
}

prepare_source_tree() {
  local build_root
  build_root="$(mktemp -d "${TMPDIR:-/tmp}/mpi4py-build.XXXXXX")"

  if [[ -f "$tarball" ]]; then
    echo "using local tarball $tarball" >&2
    tar -xzf "$tarball" -C "$build_root"
    echo "$build_root/$name"
    return 0
  fi

  if [[ -d "$name" ]]; then
    echo "using local source tree $name" >&2
    echo "$name"
    return 0
  fi

  echo "downloading $tarball" >&2
  wget "https://github.com/mpi4py/mpi4py/releases/download/$version/$tarball"
  tar -xzf "$tarball" -C "$build_root"
  echo "$build_root/$name"
}

write_mpi4py_config() {
  local purelib
  purelib="$("$python_bin" -c 'import sysconfig; print(sysconfig.get_path("purelib"))')"
  cat > "$purelib/mpi4py/mpi.cfg" <<EOF
[mpi]
mpicc = $mpicc
mpif90 = $mpif90
mpifort = $mpif90
EOF
}

write_site_precedence_override() {
  local purelib
  purelib="$("$python_bin" -c 'import sysconfig; print(sysconfig.get_path("purelib"))')"
  cat > "$purelib/00-pyranda-local-site.pth" <<'EOF'
import sys, sysconfig; purelib = sysconfig.get_path("purelib"); sys.path.insert(0, sys.path.pop(sys.path.index(purelib))) if purelib in sys.path else None
EOF
}

mpicc="$(resolve_compiler "${MPICC:-}" mpicc)"
test -n "$mpicc" || { echo "unable to find mpicc in PATH"; exit 1; }
mpif90="$(pick_mpif90 "$mpicc")"
cython_site_packages="$(find_cython_site_packages)" || {
  echo "unable to find Cython for building mpi4py"
  echo "install Cython into $env_dir or make it importable from python3"
  exit 1
}
build_support_site_packages="$(find_build_support_site_packages)" || {
  echo "unable to find setuptools for source builds"
  echo "install setuptools into $env_dir or make it importable from python3"
  exit 1
}

cython_overlay="$(mktemp -d "${TMPDIR:-/tmp}/mpi4py-cython.XXXXXX")"
cleanup() {
  rm -rf "$cython_overlay"
}
trap cleanup EXIT

for package_name in Cython cython.py pyximport; do
  if [[ -e "$cython_site_packages/$package_name" ]]; then
    ln -s "$cython_site_packages/$package_name" "$cython_overlay/$package_name"
  fi
done

build_pythonpath="$cython_overlay:$build_support_site_packages"
if [[ "$cython_site_packages" != "$build_support_site_packages" ]]; then
  build_pythonpath="$build_pythonpath:$cython_site_packages"
fi
build_pythonpath="$build_pythonpath${PYTHONPATH:+:$PYTHONPATH}"

if has_importable_mpi4py; then
  importable_mpi4py=yes
else
  importable_mpi4py=no
fi

if [[ "$mpi4py_mode" == "reuse" && "$importable_mpi4py" == "yes" ]]; then
  echo "selected mpi4py action: reuse importable mpi4py from the selected python environment"
else
  action_description="$(describe_mpi4py_source)"
  if [[ "$mpi4py_mode" == "reuse" ]]; then
    echo "selected mpi4py action: requested reuse, but no importable mpi4py was found; $action_description"
  elif [[ "$importable_mpi4py" == "yes" ]]; then
    echo "selected mpi4py action: $action_description"
    echo "  ignoring importable mpi4py because PYRANDA_MPI4PY_MODE defaults to build"
    echo "  set PYRANDA_MPI4PY_MODE=reuse to reuse an existing mpi4py"
  else
    echo "selected mpi4py action: $action_description"
  fi
  srcdir="$(prepare_source_tree)"

  echo "building mpi4py $version with:"
  echo "  python: $python_bin"
  echo "  mpicc:  $mpicc"
  echo "  mpif90: $mpif90"
  echo "  source: $srcdir"

  env \
    PYTHONPATH="$build_pythonpath" \
    MPICC="$mpicc" \
    MPILD="$mpicc" \
    CC="$mpicc" \
    LDSHARED="$mpicc -shared" \
    "$python_bin" -m pip install \
      --no-build-isolation \
      --no-cache-dir \
      --force-reinstall \
      "$srcdir"

  write_mpi4py_config
fi

echo "mpi4py available in $env_dir:"
show_mpi4py_details

write_site_precedence_override

"$python_bin" -m pip install -r requirements.txt
make -C pyranda/parcop clean python="$python_bin"
env PYTHONPATH="$build_pythonpath" "$python_bin" -m pip install --no-build-isolation .
