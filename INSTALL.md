# Installing Pyranda

## Prerequisites
Current installs assume a Python 3 environment with the following available:

- `numpy` for `f2py` and array support during the Fortran extension build,
- `scipy` and `matplotlib` for runtime dependencies used throughout the package and examples,
- `mpi4py` built against the MPI toolchain you intend to use with pyranda, and
- an MPI-enabled Fortran compiler (`mpif90` or `mpifort`) visible in `PATH` or passed explicitly at build time.

Pyranda's build imports `numpy` and `mpi4py` from `setup.py`, so those packages must already
be importable before running `pip install --no-build-isolation .`.

### installing `mpi4py`
`mpi4py` should be installed with the `--no-cache-dir` option to avoid using an
existing build with a cached compiler.

```
# if your mpi*s are in your path
pip install mpi4py --no-cache-dir
# otherwise you can specify an environment variable
env MPICC=/path/to/your/mpi pip install mpi4py --no-cache-dir
```

On `toss_4_x86_64`, use [scripts/toss_4_x86_64_ib_cray.sh] instead of a plain `pip install mpi4py`. The script:

- prefers a local `mpi4py-<version>.tar.gz` before downloading anything,
- builds `mpi4py` from source with the active MPI wrapper,
- overrides Python's default linker wrapper so MPI extension linking matches the MPI compiler toolchain, and
- writes `mpi4py/mpi.cfg` so pyranda can recover the matching `mpif90`/`mpifort` compiler during its own build.

By default, the script rebuilds `mpi4py` from source even if the selected Python can already
import it. To reuse an existing `mpi4py` instead, set `PYRANDA_MPI4PY_MODE=reuse`.

For a full fresh-environment bootstrap on this host, use [scripts/bootstrap_toss_4_x86_64_ib.sh]. It creates a new virtual environment and then runs the full `mpi4py` + `pyranda` install flow into it.

```bash
BOOTSTRAP_VENV_ARGS=--system-site-packages bash scripts/bootstrap_toss_4_x86_64_ib.sh myEnv
```

On this host, the plain bootstrap still fails later if it cannot reach the configured package
index for `matplotlib`. Using `BOOTSTRAP_VENV_ARGS=--system-site-packages` avoids that by
reusing the base Python's `numpy`, `scipy`, and `matplotlib`, while still rebuilding `mpi4py`
from source by default.


## Installing with pip

### [optional] Using a virtualenv

The standard-library venv module is sufficient:

```bash
python3 -m venv my_venv
source my_venv/bin/activate
```

You can also verify that you're in your venv by checking your `$PATH`:

```
...> echo $PATH
/path/to/your/env/root/my_venv/bin:...
```

### install pyranda

Install the Python dependencies first so `numpy` and `mpi4py` are already importable when
`setup.py` runs:

```bash
python -m pip install -r scripts/requirements.txt
```

Then install pyranda itself:

```bash
python -m pip install --no-build-isolation . [--user]
```

## Installing without pip

```
[python setup.py build [extra_build_args]]
python setup.py install
```

## Legacy Instructions - Manual Install
These notes are kept for historical reference. They predate the current Python 3 install
scripts above and still use older `setup.py`-driven workflows.

This process should work on any system and will allow for an arbitrary compiler to be used for
the fortran and for the mpi4py.

### Step 1: Ensure python and numpy

#### Python
Use a Python 3 interpreter with `numpy` available.

#### numpy
As long as numpy is working with your version of python above, there will be no
compatibility issues.  This can be installed in a number of ways. http://www.numpy.org

### Step 2: Custom install of mpi4py
This python package provides MPI bindings to python and may or may not exist on your system
and python path.

#### Install mpi4py (this should work on most systems with a mpi compiler installed)
```
export version=4.1.0
wget https://github.com/mpi4py/mpi4py/releases/download/$version/mpi4py-$version.tar.gz
tar xvzf mpi4py-$version.tar.gz
cd mpi4py-$version
python setup.py build --mpicc=/where/you/have/mpicc
python setup.py install --prefix=install_location_mpi4py
```

** Add install_location_mpi4py/*/site_packages to PYTHONPATH **

### Step 3: Pyranda build/install
A fortran compiler compatible with the mpicc used in mpi4py is used by default.
2003 and above standards enforced and MPI libraries is required.
### Install pyranda
```
git clone https://github.com/LLNL/pyranda.git
cd pyranda
python setup.py build
python setup.py install --prefix=install_location_pyranda
```

** Add install_location_pyranda/*/site_packages to PYTHONPATH **

### Step 4: Run tests to check install
For a quick smoke test, navigate to `pyranda/examples` and run
```
MPLBACKEND=Agg MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/mpl-pyranda} python advection.py 32 1
```
