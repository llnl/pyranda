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
from mpi4py import MPI
from .pyrandaPackage import pyrandaPackage
from .pyrandaMisc import ipoint


class pyrandaProbes(pyrandaPackage):
    def __init__(self, pysim, x=None, y=None, z=None):

        PackageName = "Probes"
        pyrandaPackage.__init__(self, PackageName, pysim)

        self.xpts = x
        self.ypts = y

        self.pysim = pysim
        # 3D not supported for now.
        self.zpts = None

        self.values = None

        self.points = []

        xmesh = pysim.PyMPI.ghost(pysim.x().data, 1)[:, :, 0]
        ymesh = pysim.PyMPI.ghost(pysim.y().data, 1)[:, :, 0]
        for ix, iy in zip(self.xpts, self.ypts):
            self.points.append(ipoint(ix, iy, xmesh, ymesh))

    def getOld(self, var, method="linear"):

        # Get the pysim data
        values = []
        vdata = self.pysim.PyMPI.ghost(self.pysim.variables[var].data, 1)
        vals = vdata[:, :, 0]

        # Get interpolated values
        for pt in self.points:
            values.append(pt.interp(vals, method))

        # Communicate
        Lvalues = numpy.array(values)
        Gvalues = self.pysim.PyMPI.comm.allreduce(Lvalues, op=MPI.SUM)

        self.values = Gvalues
        return self.values

    def get(self, var):

        # Get the pysim data
        values = []
        vdata = self.pysim.PyMPI.ghost(self.pysim.variables[var].data, 1)
        vals = vdata[:, :, 0]

        # Get interpolated values
        for pt in self.points:
            values.append(pt.interpFast(vals))

        # Communicate
        Lvalues = numpy.array(values)
        Gvalues = self.pysim.PyMPI.comm.allreduce(Lvalues, op=MPI.SUM)

        self.values = Gvalues
        return self.values

    def plot(self, var=None, style=None):

        # Plot the probes on the current figure
        # if not var:
        #    self.get(var)

        # vals = self.values
        x = self.xpts
        y = self.ypts

        if self.pysim.PyMPI.master:
            import matplotlib.pyplot as plt

            if style:
                plt.plot(x, y, style)
            else:
                plt.plot(x, y, "ko")

            plt.pause(0.1)
