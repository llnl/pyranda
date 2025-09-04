#!/bin/bash

test -e INSTALL.md || { echo "this must run from the root of pyranda"; exit 1; }

version=4.1.0
name=mpi4py-$version
mpicc=/usr/tce/packages/cray-mpich/cray-mpich-9.0.1-cce-20.0.0-magic/bin/mpicc


if [[ ! -d $name ]]; then
  wget https://github.com/mpi4py/mpi4py/releases/download/$version/$name.tar.gz
  tar xvzf $name.tar.gz
fi

if [[ ! $(pip list | grep mpi4py) ]]; then
  pushd $name
  python3 setup.py build --mpicc=$mpicc
  python3 setup.py install
  popd
fi

pip install -r requirements.txt

pip install .
