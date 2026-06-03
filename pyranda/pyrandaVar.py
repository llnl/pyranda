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


class pyrandaVar:
    def __init__(self, name, kind, rank):

        self.name = name
        self.kind = kind
        self.rank = rank

        self.rhs = None  # RHS for algebraic closure
        self.dt = None  # dt for temporal integration

        self.data = None  # Need to allocate
        self.allocated = False
        self.PyMPI = None

    def __allocate__(self, pympi):  # Make this private

        if self.allocated:
            return

        self.PyMPI = pympi

        # Inherit the mesh size
        if self.rank == "scalar":
            self.data = pympi.emptyScalar()
        elif self.rank == "vector":
            self.data = pympi.emptyVector()
        else:
            raise ValueError("Error: rank: %s not recognized" % self.rank)

        self.allocated = True

    def sum(self, axis=None, idata=None):

        if type(idata) is not type(self.data):
            idata = self.data

        # Vector variables not  fully supported with axis
        if self.rank == "vector":
            return [
                self.PyMPI.sum3D(idata[:, :, :, 0]),
                self.PyMPI.sum3D(idata[:, :, :, 1]),
                self.PyMPI.sum3D(idata[:, :, :, 2]),
            ]

        # For no axis, return global sum
        if axis is None:
            return self.PyMPI.sum3D(idata)

        # Sum along 2 directions
        if type(axis) is type(()) or type(axis) is type([]):
            if len(axis) == 2:
                if (0 in axis) and (1 in axis):
                    return self.PyMPI.xysum(idata)
                elif (0 in axis) and (2 in axis):
                    return self.PyMPI.xzsum(idata)
                elif (1 in axis) and (2 in axis):
                    return self.PyMPI.yzsum(idata)
            else:
                self.PyMPI.iprint(
                    "Error: length of axis argument cant be larger than 2"
                )
                return None
        # Sum along 1 directions
        elif type(axis) is type(0):
            if axis == 0:
                return self.PyMPI.xsum(idata)
            if axis == 1:
                return self.PyMPI.ysum(idata)
            if axis == 2:
                return self.PyMPI.zsum(idata)
        else:
            self.PyMPI.iprint("Error: unknown type given as axis")
            return None

    def mean(self, axis=None, idata=None):

        npts = float(self.PyMPI.nx * self.PyMPI.ny * self.PyMPI.nz)
        rdata = self.sum(axis=axis, idata=idata)
        return rdata * rdata.size / npts

    def var(self, axis=None):

        bar = self.mean(axis=axis)
        bar2 = self.mean(axis=axis, idata=self.data**2)

        return bar2 - bar * bar

    # Allow global access via indices
    def __getitem__(self, given):

        # Check length of indices
        slen = len(given)

        # Only 3D arguments
        if slen != 3:
            self.PyMPI.iprint("Error: 3 index ranges must be given")
            return None

        indices = [None] * 3

        NX = [self.PyMPI.nx, self.PyMPI.ny, self.PyMPI.nz]
        for i in range(slen):
            if isinstance(given[i], slice):
                start = given[i].start
                stop = given[i].stop
                if start is None:
                    start = 0
                if stop is None:
                    stop = NX[i]
            else:
                start = given[i]
                stop = start + 1
            indices[i] = [start, stop]

        gdata = self.PyMPI.getIJK(self.data, indices[0], indices[1], indices[2])
        return numpy.squeeze(gdata)
