# Installing Pyranda

## Prerequisites
Before installing pyranda you'll need to install numpy and mpi4py. Numpy's f2py utility
is used to build the python-fortran interface and we try to validate that your choice of
compiler will work with the installed mpi4py's version of mpi.

### installing `mpi4py`
`mpi4py` should be installed with the `--no-cache-dir` option to avoid using an
existing build with a cached compiler.

```
# if your mpi*s are in your path
pip install mpi4py --no-cache-dir
# otherwise you can specify an environment variable
env MPICC=/path/to/your/mpi pip install mpi4py --no-cache-dir
```

On `toss_4_x86_64`, use [scripts/toss_4_x86_64_ib_cray.sh](/usr/WS2/olson45/SAND/pyranda-codex/scripts/toss_4_x86_64_ib_cray.sh) instead of a plain `pip install mpi4py`. The script:

- prefers a local `mpi4py-<version>.tar.gz` before downloading anything,
- builds `mpi4py` from source with the active MPI wrapper,
- overrides Python's default linker wrapper so MPI extension linking matches the MPI compiler toolchain, and
- writes `mpi4py/mpi.cfg` so pyranda can recover the matching `mpif90`/`mpifort` compiler during its own build.

For a full fresh-environment bootstrap on this host, use [scripts/bootstrap_toss_4_x86_64_ib.sh](/usr/WS2/olson45/SAND/pyranda-codex/scripts/bootstrap_toss_4_x86_64_ib.sh). It creates a new virtual environment and then runs the full `mpi4py` + `pyranda` install flow into it.

If your host cannot reach the configured Python package index but already has compatible `numpy`/`scipy`/`matplotlib` installed system-wide, you can create the venv with inherited site packages:

```bash
BOOTSTRAP_VENV_ARGS=--system-site-packages bash scripts/bootstrap_toss_4_x86_64_ib.sh myEnv
```

### installing `numpy`
`numpy` shouldn't have any special build steps, just install as normal:

```
pip install numpy
```

## Installing with pip

### [optional] Using a virtualenv

```
...> virtualenv [-p python] my_venv
...> source my_venv/bin/activate
(my_venv) ...>
```

You can also verify that you're in your venv by checking your `$PATH`:

```
...> echo $PATH
/path/to/your/env/root/my_venv/bin:...
```

### install pyranda

```
pip install --no-build-isolation . [--user]
```

###

## Installing without pip

```
[python setup.py build [extra_build_args]]
python setup.py install
```

## Legacy Instructions - Manual Install
This process should work on any system and will allow for an arbitrary compiler to be used for
the fortran and for the mpi4py.

### Step 1: Ensure python and numpy

#### Python
Though other versions of python may very well work, we recommend and support
python 2.7, 3.5, and 3.6 for pyranda.

#### numpy
As long as numpy is working with your version of python above, there will be no
compability issues.  This can be installed in a number of ways. http://www.numpy.org

### Step 2: Custom install of mpi4py
This python package provides MPI bindings to python and may or may not exists on your system
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
Trying navigating to pyranda/examples and running
```
python advection.py
```
