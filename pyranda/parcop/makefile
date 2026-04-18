#===================================================================================================
# Copyright (c) 2018, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
#
# Written by: Britton J. Olson, olson45@llnl.gov
#===================================================================================================
#

# return the first argument if it's defined, otherwise return the second argument
# usage: $(call use-if,a,b)
use-if = $(if $1,$1,$2)

chaos_5_x86_64_ib.gnu.mpif90   = /collab/usr/gapps/python/build/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.4.7/mpich-3.0.4-evwtvwqv5i2ujysfnivcdtrnppmpcsrg/bin/mpif90
chaos_5_x86_64_ib.intel.mpif90 = /usr/local/tools/mvapich2-intel-2.2/bin/mpif90
chaos_5_x86_64_ib.python       = /usr/local/bin/python

toss_3_x86_64_ib.gnu.mpif90   = /usr/tce/packages/mvapich2/mvapich2-2.2-gcc-6.1.0/bin/mpif90
toss_3_x86_64_ib.intel.mpif90 = /usr/tce/packages/mvapich2/mvapich2-2.2-intel-16.0.3/bin/mpif90
toss_3_x86_64_ib.clang.mpif90 = /usr/tce/packages/mvapich2/mvapich2-2.2-clang-4.0.0/bin/mpif90
toss_3_x86_64_ib.pgi.mpif90   = /usr/tce/packages/mvapich2/mvapich2-2.2-pgi-18.1/bin/mpif90
toss_3_x86_64_ib.python       = /usr/tce/bin/python

blueos_3_ppc64le_ib.ibm.mpif90 = /usr/tce/packages/spectrum-mpi/spectrum-mpi-2018.04.27-xl-beta-2018.06.01/bin/mpif90
blueos_3_ppc64le_ib.gnu.mpif90 = /usr/tcetmp/bin/mpigfortran
blueos_3_ppc64le_ib.python     = /usr/tcetmp/bin/python

gopt = -funroll-loops -Ofast
gnu.fflags   = -fPIC -ffree-form -ffree-line-length-0 -fbackslash $(gopt) -std=legacy
intel.fflags = -fPIC -O3 -cpp
# it's likely that mpiclang is using gfortran
clang.fflags = $(gnu.fflags)
ibm.linkflags = -qoffload -W@,-v -qsmp=omp -qinfo=omperrtrace
ibm.fflags   = -fPIC -qxlf2003=polymorphic $(ibm.linkflags) -qsuffix=cpp=f90:f=f90
ibm.lopts    = -qmkshrobj
pgi.fflags   = -fPIC
cray.fflags  = -ffree -efF -fPIC

sys_type = $(call use-if,$(SYS_TYPE),$(shell uname))
python  ?= $(call use-if,$($(compiler).$(sys_type).python),$(shell which python))
mpif90  ?= $(call use-if,$($(sys_type).$(compiler).mpif90),$(shell which mpif90))
fcomp   ?= $(call use-if,$(compiler),$(shell ./vendor.sh $(mpif90)))
mpicc   ?= $(dir $(mpif90))mpicc
ccomp   ?= $(call use-if,$(compiler),$(shell ./vendor.sh $(mpicc)))
# make sure we use this python's f2py unless the user wants something else
f2py    ?= $(python) -c 'import numpy.f2py; numpy.f2py.main()'
fflags  ?= $($(fcomp).fflags)
lopts   ?= $(call use-if,$($(fcomp).lopts),-shared)

#
# Python and numpy directories for use when building the f2py wrapped module
#
python.Iflags  = -I$(shell $(python) -c "from sysconfig import get_config_var; print(get_config_var('INCLUDEPY'))")
python.libdir  = $(shell $(python) -c "from sysconfig import get_config_var; print(get_config_var('LIBDIR'))")
python.version = $(shell $(python) -c "from sysconfig import get_config_var; print(get_config_var('VERSION'))")
python.abiflag = $(shell $(python) -c "import sys; print(getattr(sys,'abiflags',''))")
# TODO: how should we do this so it's portable?
#python.Lflags  = -L$(python.libdir) -Wl,-R$(python.libdir)
python.Lflags  = -L$(python.libdir) -Wl,-rpath=$(python.libdir)
python.lflags  = -lpython$(python.version)$(python.abiflag)
mpi.rpathflags = $(shell $(mpif90) -show 2>/dev/null | tr ' ' '\n' | awk '/^-L/ {sub(/^-L/,""); printf "-Wl,-rpath=%s ", $$0}')

numpy.dir    = $(shell $(python) -c "import numpy; print(numpy.__path__[0])")
numpy.include = $(shell $(python) -c "import numpy; print(numpy.get_include())")
numpy.Iflags = -I$(numpy.include)
f2py.dir     = $(numpy.dir)/f2py/src
f2py.Iflags  = -I$(f2py.dir)

#
# source order is important because of mod file dependencies
#
source = blockmath.o ompsync.o stencils.o nrtype.o nrutil.o patch.o pentadiagonal.o \
	     comm.o compact_basetype.o compact_d1.o compact_r3.o compact_r4.o \
             compact.o mesh.o objects.o compact_operators.o operators.o

vpath %.c $(f2py.dir)

.PHONY: clean info

%.o: %.f90
	$(mpif90) $(fflags) -c $<

# TODO: is `-Wno-implicit...` portable?
%.o: %.c
	$(mpicc) -fPIC $(python.Iflags) $(numpy.Iflags) $(f2py.Iflags) -Wno-implicit-function-declaration -c $^

all: parcop.so test_pent

libparcop.a: $(source)
	@ar -ur $@ $(source)

#
# TODO: is the '2' in of the f2pywrappers2.f90 the major python version?
#
# f2py ended up not working well with certain compilers (ibm, pgi) so
# instead build the wrapper components using above flags
#
parcopmodule.c parcop-f2pywrappers2.f90: libparcop.a parcop.o
	@$(f2py) --lower --wrap-functions -m parcop parcop.f90 --verbose

parcop.so: parcopmodule.o fortranobject.o parcop-f2pywrappers2.o
	$(mpif90) $(lopts) $(python.Lflags) $(mpi.rpathflags) \
		parcopmodule.o fortranobject.o parcop.o parcop-f2pywrappers2.o libparcop.a $(python.lflags) \
		-o parcop.so

test_pent: libparcop.a parcop.o test_pent.o
	$(mpif90) test_pent.o parcop.o libparcop.a -o $@

miniApp: libparcop.a parcop.o miniApp.o
	$(mpif90) miniApp.o parcop.o libparcop.a -o $@

miniAppGPU: libparcop.a parcop.o miniApp.o
	$(mpif90) $(ibm.linkflags) miniApp.o parcop.o libparcop.a -o $@

info:
	@echo "sys_type:  $(sys_type)"
	@echo "python:    $(python)"
	@echo "mpif90:    $(mpif90)"
	@echo "mpicc:     $(mpicc)"
	@echo "f2py:      $(f2py)"
	@echo "fflags:    $(fflags)"
	@echo "mpi.rpath: $(mpi.rpathflags)"
	@echo "fcompiler: $(fcomp)"
	@echo "ccompiler: $(ccomp)"

clean:
	rm -rf *.o *.mod test_pent *.so *.pyf *.pyc libparcop.a *.dSYM parcopmodule.c parcop-f2pywrappers*.f90
