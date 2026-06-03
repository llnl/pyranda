################################################################################
# Copyright (c) 2018, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
################################################################################
import os
import sys
import shutil
import subprocess
from setuptools import setup
from distutils.command.build import build
from distutils.command.install import install

try:
    import mpi4py
    from numpy import f2py
except ImportError as e:
    print(e)
    print("some modules need to be installed before building pyranda")


def find_mpi4py_mpif90_compiler():
    # try to get the mpif90 compiler from mpi4py, it isn't always there..
    # if a user installs like `env MPICC=/path/to/mpicc pip install mpi4py ...`
    # it doesn't seem to have anything other than mpicc
    mpi4py_compilers = mpi4py.get_config()
    if "mpif90" in mpi4py_compilers:
        return mpi4py_compilers["mpif90"]
    elif "mpifort" in mpi4py_compilers:
        return mpi4py_compilers["mpifort"]

    # mpi4py >= 4 returns an empty config by default unless a packager ships
    # mpi.cfg, so fall back to the environment and PATH.
    mpicc = (
        mpi4py_compilers.get("mpicc")
        or os.environ.get("MPICC")
        or shutil.which("mpicc")
    )
    if mpicc:
        mpicc_dir = os.path.dirname(mpicc)
        for sibling in ("mpif90", "mpifort"):
            candidate = os.path.join(mpicc_dir, sibling)
            if os.path.exists(candidate):
                return candidate

    for envvar in ("MPIF90", "MPIFORT"):
        candidate = os.environ.get(envvar)
        if candidate:
            return candidate

    return shutil.which("mpif90") or shutil.which("mpifort")


def find_mpi4py_mpicc_compiler():
    mpi4py_compilers = mpi4py.get_config()
    return (
        mpi4py_compilers.get("mpicc")
        or os.environ.get("MPICC")
        or shutil.which("mpicc")
    )


def verify_mpi_compiler_compatibility(mpif90):
    mpi4py_compilers = mpi4py.get_config()
    configured_mpif90 = mpi4py_compilers.get("mpif90") or mpi4py_compilers.get(
        "mpifort"
    )
    configured_mpicc = find_mpi4py_mpicc_compiler()
    resolved_mpif90 = os.path.realpath(mpif90) if mpif90 else None

    if not resolved_mpif90:
        raise RuntimeError("unable to determine an mpif90/mpifort compiler for pyranda")

    if configured_mpif90:
        if resolved_mpif90 != os.path.realpath(configured_mpif90):
            raise RuntimeError(
                "mpi4py was built with {!r} but pyranda will build with {!r}".format(
                    configured_mpif90, mpif90
                )
            )
        return

    if configured_mpicc:
        mpicc_dir = os.path.dirname(os.path.realpath(configured_mpicc))
        for sibling in ("mpif90", "mpifort"):
            candidate = os.path.join(mpicc_dir, sibling)
            if os.path.exists(candidate):
                if resolved_mpif90 != os.path.realpath(candidate):
                    raise RuntimeError(
                        "mpi4py was built with {!r}, but pyranda will build with {!r}".format(
                            configured_mpicc, mpif90
                        )
                    )
                return

    print(
        "mpi4py did not advertise compiler metadata; using {!r} without additional verification".format(
            mpif90
        )
    )


distname = "pyranda"
fortran_module = "parcop"
fortran_package = "{}/{}".format(distname, fortran_module)
fortran_module_lib = "{}.so".format(fortran_module)
packages = [distname, fortran_package]
install_requires = ["matplotlib"]
description = """
A Python driven, Fortran powered Finite Difference solver for arbitrary
hyperbolic PDE systems. This is the mini-app for the Miranda code.
"""


class PyrandaMakeMixin:
    user_options = [
        ("mpif90=", None, "mpif90 compiler"),
        ("fflags=", None, "flags for mpif90"),
        (
            "mpiexec=",
            None,
            "mpi exec command used to verify when verifying the mpi compiler",
        ),
        (
            "numprocs-arg=",
            None,
            "mpi exec num procs arg used when verifying the mpi compiler",
        ),
        ("numprocs=", None, "number of procs used when verifying mpi"),
        (
            "check-mpi-compatibility",
            None,
            "check that the mpi compiler is compatible with mpi4py",
        ),
    ]

    def initialize_options(self):
        self.mpif90 = None
        self.fflags = None
        self.mpiexec = None
        self.numprocs_arg = None
        self.numprocs = None
        self.check_mpi_compatibility = None

    def finalize_options(self):
        pass

    def test(self):
        print("running regression tests for pyranda...")
        try:
            os.chdir("tests")
            subprocess.check_call(["python", "run_tests.py"])
        except Exception:
            print("Failed to run tests")
            raise

    def clean(self):
        print("cleaning up from {} build".format(fortran_module))
        try:
            subprocess.check_call(["make", "-C", fortran_package, "clean"])
        except subprocess.CalledProcessError:
            print("failed to clean {}".format(fortran_module))
            raise
        print("{} cleaned".format(fortran_module))

    def run(self):
        mpi4py_mpif90 = find_mpi4py_mpif90_compiler()
        python = sys.executable
        args = ["make", "-C", fortran_package, "python={}".format(python)]
        selected_mpif90 = None

        if self.mpif90 is not None:
            args.append("mpif90={}".format(self.mpif90))
            selected_mpif90 = self.mpif90
        elif mpi4py_mpif90 is not None:
            args.append("mpif90={}".format(mpi4py_mpif90))
            selected_mpif90 = mpi4py_mpif90

        # test mpi X mpi4py compatibility
        if self.check_mpi_compatibility:
            print("trying to verify that the mpi compiler is compatible with mpi4py")
            try:
                verify_mpi_compiler_compatibility(selected_mpif90)
            except RuntimeError:
                print("mpi verification failed")
                raise
            print("mpi verification successful")

        # build the .so module
        print("building {}".format(fortran_module))
        try:
            if self.fflags is not None:
                args.append("fflags={}".format(self.fflags))
            if self.numprocs_arg is not None:
                args.append("np_arg={}".format(self.numprocs_arg))
            if self.mpiexec is not None:
                args.append("mpirun={}".format(self.mpiexec))
            if self.numprocs is not None:
                args.append("numprocs={}".format(self.numprocs))

            args.append(fortran_module_lib)

            subprocess.check_call(args)
        except subprocess.CalledProcessError:
            print("failed to build {}".format(fortran_module))
            raise
        print("{} built".format(fortran_module))


class BuildPyranda(build, PyrandaMakeMixin):
    user_options = build.user_options + PyrandaMakeMixin.user_options

    def initialize_options(self):
        build.initialize_options(self)
        PyrandaMakeMixin.initialize_options(self)

    def finalize_options(self):
        build.finalize_options(self)
        PyrandaMakeMixin.finalize_options(self)

    def run(self):
        PyrandaMakeMixin.run(self)
        build.run(self)


class InstallPyranda(install, PyrandaMakeMixin):
    def initialize_options(self):
        install.initialize_options(self)

    def finalize_options(self):
        install.finalize_options(self)

    def run(self):
        install.run(self)
        # PyrandaMakeMixin.clean(self)


class CleanPyranda(install, PyrandaMakeMixin):
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        PyrandaMakeMixin.clean(self)


class TestPyranda(install, PyrandaMakeMixin):
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        PyrandaMakeMixin.test(self)


setup_args = dict(
    name=distname,
    description=description,
    packages=packages,
    package_data={distname: ["*.tex"], fortran_package: [fortran_module_lib]},
    install_requires=install_requires,
    cmdclass={
        "build": BuildPyranda,
        "install": InstallPyranda,
        "clean": CleanPyranda,
        "runtest": TestPyranda,
    },
)

setup(**setup_args)
