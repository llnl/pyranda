################################################################################
# Copyright (c) 2018, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
#
# Written by: Britton J. Olson, olson45@llnl.gov
################################################################################
import numpy
from .pyrandaUtils import fortran3d, splitLines
from . import parcop
from .pyrandaVar import pyrandaVar


class pyrandaMesh:
    def __init__(self):

        self.name = "base"
        self.options = {}
        self.dims = 0
        self.PyMPI = None
        self.x = None
        self.y = None
        self.z = None
        self.coordsys = None
        self.shape = None
        self.function = None
        # self.makeMesh()

    def makeMeshStr(self, str_mesh):

        meshStrLines = splitLines(str_mesh)  # Split into lines

        self.options = defaultMeshOptions()
        self.get_sMap()
        local_scope = {"self": self}
        for msl in meshStrLines:
            exec(fortran3d(msl, self.sMap), globals(), local_scope)

    def set_options(self, ind, x1, xn, nn, periodic=False):

        self.options["x1"][ind] = x1
        self.options["xn"][ind] = xn
        self.options["nn"][ind] = nn
        self.options["periodic"][ind] = periodic

    def get_sMap(self):
        sMap = {}
        sMap["xdom=("] = "self.set_options(0,"
        sMap["ydom=("] = "self.set_options(1,"
        sMap["zdom=("] = "self.set_options(2,"
        self.sMap = sMap

    def makeMesh(self, xr=None, yr=None, zr=None):

        # Per processor (local) mesh size
        ax, ay, az = self.PyMPI.ax, self.PyMPI.ay, self.PyMPI.az

        # Compute global index arrays
        self.indices = [
            pyrandaVar("iloc", "mesh", "scalar"),
            pyrandaVar("jloc", "mesh", "scalar"),
            pyrandaVar("kloc", "mesh", "scalar"),
        ]

        for ind in self.indices:
            ind.__allocate__(self.PyMPI)

        for i in range(ax):
            for j in range(ay):
                for k in range(az):
                    # Get global indices
                    ii = i + self.PyMPI.chunk_3d_lo[0]
                    jj = j + self.PyMPI.chunk_3d_lo[1]
                    kk = k + self.PyMPI.chunk_3d_lo[2]
                    self.indices[0].data[i, j, k] = ii
                    self.indices[1].data[i, j, k] = jj
                    self.indices[2].data[i, j, k] = kk

        options = self.options
        self.nn = options["nn"]

        if self.coordsys != 3:
            # Define the grid here (in Fortran)
            parcop.parcop.setup_mesh(self.PyMPI.patch, self.PyMPI.level)

        elif self.coordsys == 3:  # General/curvilinear - need grid from python
            x1 = options["x1"]
            xn = options["xn"]
            nn = options["nn"]

            chunk_lo = self.PyMPI.chunk_3d_lo
            chunk_hi = self.PyMPI.chunk_3d_hi

            x = numpy.linspace(x1[0], xn[0], num=nn[0])[chunk_lo[0] : chunk_hi[0] + 1]
            y = numpy.linspace(x1[1], xn[1], num=nn[1])[chunk_lo[1] : chunk_hi[1] + 1]
            z = numpy.linspace(x1[2], xn[2], num=nn[2])[chunk_lo[2] : chunk_hi[2] + 1]

            x, y, z = numpy.meshgrid(x, y, z, indexing="ij")

            x = numpy.asfortranarray(x)
            y = numpy.asfortranarray(y)
            z = numpy.asfortranarray(z)

            # if restart coords are given, use them
            if (
                isinstance(xr, numpy.ndarray)
                and isinstance(yr, numpy.ndarray)
                and isinstance(zr, numpy.ndarray)
            ):
                x = xr
                y = yr
                z = zr
            else:
                # Evaluate an ijk function for the mesh
                if "function" in self.options:
                    self.function = self.options["function"]
                    if self.function:
                        for i in range(ax):
                            for j in range(ay):
                                for k in range(az):
                                    # Get global indices
                                    ii = i + self.PyMPI.chunk_3d_lo[0]
                                    jj = j + self.PyMPI.chunk_3d_lo[1]
                                    kk = k + self.PyMPI.chunk_3d_lo[2]
                                    x[i, j, k], y[i, j, k], z[i, j, k] = self.function(
                                        ii, jj, kk
                                    )

            # Define the grid here (send to fortran
            mesh_per = options["periodicGrid"]
            parcop.parcop.setup_mesh_x3(
                self.PyMPI.patch, self.PyMPI.level, x, y, z, mesh_per
            )

        # Read in from the fortran and set to numpy arrays
        self.coords = [
            pyrandaVar("x", "mesh", "scalar"),
            pyrandaVar("y", "mesh", "scalar"),
            pyrandaVar("z", "mesh", "scalar"),
        ]
        for cor in self.coords:
            cor.__allocate__(self.PyMPI)

        self.coords[0].data = parcop.parcop.xgrid(ax, ay, az)
        self.coords[1].data = parcop.parcop.ygrid(ax, ay, az)
        self.coords[2].data = parcop.parcop.zgrid(ax, ay, az)

        # Write in x-y-z space
        if self.coordsys == 1:  # Cylindrical
            x = self.coords[0].data * numpy.cos(self.coords[1].data)
            y = self.coords[0].data * numpy.sin(self.coords[1].data)
            z = self.coords[2].data
        elif self.coordsys == 2:  # Spherical
            r = self.coords[0].data
            th = self.coords[1].data
            phi = self.coords[2].data

            x = r * numpy.sin(th) * numpy.cos(phi)
            y = r * numpy.sin(th) * numpy.sin(phi)
            z = r * numpy.cos(th)

        elif self.coordsys == 4:  # Spherical- ratio'd r-spacing
            r = self.coords[0].data
            th = self.coords[1].data
            phi = self.coords[2].data

            x = r * numpy.sin(th) * numpy.cos(phi)
            y = r * numpy.sin(th) * numpy.sin(phi)
            z = r * numpy.cos(th)

        if self.coordsys in [1, 2, 4]:
            # Save native coords
            self.coordsNative = [
                pyrandaVar("x", "mesh", "scalar"),
                pyrandaVar("y", "mesh", "scalar"),
                pyrandaVar("z", "mesh", "scalar"),
            ]
            for cor in self.coordsNative:
                cor.__allocate__(self.PyMPI)

            self.coordsNative[0].data[:, :, :] = self.coords[0].data
            self.coordsNative[1].data[:, :, :] = self.coords[1].data
            self.coordsNative[2].data[:, :, :] = self.coords[2].data

            self.coords[0].data[:, :, :] = x
            self.coords[1].data[:, :, :] = y
            self.coords[2].data[:, :, :] = z

        self.shape = list(self.coords[0].data.shape)

        self.d1 = parcop.parcop.dxgrid(ax, ay, az)
        self.d2 = parcop.parcop.dygrid(ax, ay, az)
        self.d3 = parcop.parcop.dzgrid(ax, ay, az)

        # Mesh data
        self.CellVol = parcop.parcop.mesh_getcellvol(
            self.PyMPI.ax, self.PyMPI.ay, self.PyMPI.az
        )
        self.GridLen = parcop.parcop.mesh_getgridlen(
            self.PyMPI.ax, self.PyMPI.ay, self.PyMPI.az
        )


def defaultMeshOptions():
    options = {}

    options["x1"] = [0.0, 0.0, 0.0]
    options["xn"] = [1, 1, 1]
    options["nn"] = [1, 1, 1]
    options["periodic"] = [False, False, False]
    options["periodicGrid"] = True  # Recompute
    options["symmetric"] = [[False, False], [False, False], [False, False]]
    options["dim"] = 1
    options["coordsys"] = 0
    options["function"] = None
    options["id"] = 0

    # Option to enable manual control of domain processor decomposition
    # Set to tell Pyranda to just use pn = [ px,py,pz ] as number of processors in x, y, and z
    # instead of attempting to decompose automatically
    options["pn"] = [0, 0, 0]  # Number of processors to use in x, y, and z

    return options
