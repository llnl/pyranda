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
from mpi4py import MPI
import numpy
import sys
from functools import reduce

from .pyrandaVar import pyrandaVar
from . import parcop

# Global patch/levels to create new instances
GLOBAL_PATCH = 0
GLOBAL_LEVEL = 0
LAST_PATCH = -1
LAST_LEVEL = 0
PATCH_MAX = 100
LEVEL_MAX = 10


class pyrandaMPI:
    def __init__(self, mesh, comm=None):

        meshOptions = mesh.options
        self.nx = meshOptions["nn"][0]
        self.ny = meshOptions["nn"][1]
        self.nz = meshOptions["nn"][2]

        x1 = meshOptions["x1"][0]
        xn = meshOptions["xn"][0]
        y1 = meshOptions["x1"][1]
        yn = meshOptions["xn"][1]
        z1 = meshOptions["x1"][2]
        zn = meshOptions["xn"][2]

        dx = (meshOptions["xn"][0] - meshOptions["x1"][0]) / max(self.nx - 1, 1)
        dy = (meshOptions["xn"][1] - meshOptions["x1"][1]) / max(self.ny - 1, 1)
        dz = (meshOptions["xn"][2] - meshOptions["x1"][2]) / max(self.nz - 1, 1)
        self.dx = dx
        self.dy = dy
        self.dz = dz

        periodic = meshOptions["periodic"]
        symmetric = meshOptions["symmetric"]
        self.periodic = periodic

        self.coordsys = meshOptions["coordsys"]

        # Get the mesh function if it exists
        try:
            self.meshFunction = meshOptions["meshFunction"]
        except Exception:
            self.meshFunction = None

        # If explicit comm it given, use that communicator
        if comm:
            self.comm = comm
        else:
            self.comm = MPI.COMM_WORLD

        self.fcomm = self.comm.py2f()

        self.order = (10, 10, 10)
        self.filter_type = ("compact", "compact", "compact")

        # Check if we're doing a manual decomposition
        doManualDecompose = (
            ((meshOptions["pn"][0] >= 1) or (self.nx == 1))
            and ((meshOptions["pn"][1] >= 1) or (self.ny == 1))
            and ((meshOptions["pn"][2] >= 1) or (self.nz == 1))
        )

        if doManualDecompose == True:
            # Decompose based on user-supplied arguments
            px = max(meshOptions["pn"][0], 1)
            py = max(meshOptions["pn"][1], 1)
            pz = max(meshOptions["pn"][2], 1)
        else:
            # Compute the domain decomposition
            [px, py, pz] = decomp(self.comm.Get_size(), self.nx, self.ny, self.nz)

        if px * py * pz != self.comm.Get_size():
            if self.comm.Get_rank() == 0:
                print("Problem size = (%s,%s,%s)" % (self.nx, self.ny, self.nz))
                print("Tried (px,py,pz) = (%s,%s,%s)" % (px, py, pz))
                print("Nprocs = %s" % (px * py * pz))
                print("Available procs = %s" % self.comm.Get_size())
                print("Error: Processor under/over utilization not allowed")
                exit()
            else:
                print("Problem size = (%s,%s,%s)" % (self.nx, self.ny, self.nz))
                print("  (px,py,pz) = (%s,%s,%s)" % (px, py, pz))

        bx1 = "NONE"
        bxn = "NONE"
        by1 = "NONE"  # "PERI"
        byn = "NONE"  # "PERI"
        bz1 = "NONE"  # "PERI"
        bzn = "NONE"  # "PERI"

        # Operator dependent BCs need to be set at the solver level
        if periodic[0]:
            bx1 = "PERI"
            bxn = "PERI"
        if periodic[1]:
            by1 = "PERI"
            byn = "PERI"
        if periodic[2]:
            bz1 = "PERI"
            bzn = "PERI"

        symbc = "SYMM"
        if symmetric[0][0]:
            bx1 = symbc
        if symmetric[0][1]:
            bxn = symbc
        if symmetric[1][0]:
            by1 = symbc
        if symmetric[1][1]:
            byn = symbc
        if symmetric[2][0]:
            bz1 = symbc
        if symmetric[2][1]:
            bzn = symbc

        # Set at unique patch/level for each instance
        global GLOBAL_PATCH
        global LAST_PATCH
        global GLOBAL_LEVEL
        global LAST_LEVEL

        GLOBAL_PATCH = LAST_PATCH + 1
        if GLOBAL_PATCH >= PATCH_MAX:
            LAST_LEVEL = LAST_LEVEL + 1
            GLOBAL_PATCH = 0

        GLOBAL_LEVEL = int(LAST_LEVEL)
        self.patch = int(GLOBAL_PATCH)
        self.level = int(GLOBAL_LEVEL)
        LAST_LEVEL = self.level
        LAST_PATCH = self.patch

        parcop.parcop.setup(
            self.patch,
            self.level,
            self.fcomm,
            self.nx,
            self.ny,
            self.nz,
            px,
            py,
            pz,
            self.coordsys,
            x1,
            xn,
            y1,
            yn,
            z1,
            zn,
            bx1,
            bxn,
            by1,
            byn,
            bz1,
            bzn,
        )

        self.chunk_3d_size = numpy.zeros(3, dtype=numpy.int32, order="F")
        self.chunk_3d_lo = numpy.zeros(3, dtype=numpy.int32, order="F")
        self.chunk_3d_hi = numpy.zeros(3, dtype=numpy.int32, order="F")

        self.setPatch()
        self.xcom = MPI.Comm.f2py(parcop.parcop.commx())
        self.ycom = MPI.Comm.f2py(parcop.parcop.commy())
        self.zcom = MPI.Comm.f2py(parcop.parcop.commz())
        self.xycom = MPI.Comm.f2py(parcop.parcop.commxy())
        self.yzcom = MPI.Comm.f2py(parcop.parcop.commyz())
        self.xzcom = MPI.Comm.f2py(parcop.parcop.commxz())

        # Shifted comms
        (self.xcom_lo, self.xcom_hi) = self.xcom.Shift(0, 1)
        (self.ycom_lo, self.ycom_hi) = self.ycom.Shift(0, 1)
        (self.zcom_lo, self.zcom_hi) = self.zcom.Shift(0, 1)

        self.chunk_3d_size[0] = self.ax = int(self.nx / px)
        self.chunk_3d_size[1] = self.ay = int(self.ny / py)
        self.chunk_3d_size[2] = self.az = int(self.nz / pz)

        self.px = px
        self.py = py
        self.pz = pz

        self.chunk_3d_lo[0] = self.xcom.rank * self.chunk_3d_size[0] + 1
        self.chunk_3d_hi[0] = (self.xcom.rank + 1) * self.chunk_3d_size[0]

        self.chunk_3d_lo[1] = self.ycom.rank * self.chunk_3d_size[1] + 1
        self.chunk_3d_hi[1] = (self.ycom.rank + 1) * self.chunk_3d_size[1]

        self.chunk_3d_lo[2] = self.zcom.rank * self.chunk_3d_size[2] + 1
        self.chunk_3d_hi[2] = (self.zcom.rank + 1) * self.chunk_3d_size[2]

        self.chunk_3d_lo = self.chunk_3d_lo - 1  # Convert to 0 based indexing
        self.chunk_3d_hi = self.chunk_3d_hi - 1  # Convert to 0 based indexing

        procMap = {}
        rank = self.comm.rank
        procMap["%s-g1" % rank] = self.chunk_3d_lo
        procMap["%s-gn" % rank] = self.chunk_3d_hi

        # Custom comm operator for the processor map; make it, use it, free it (max Op count is 32)
        counterSumOp = MPI.Op.Create(addCounter, commute=True)
        self.procMap = self.comm.allreduce(procMap, op=counterSumOp)
        counterSumOp.Free()

        self.der = parcop_der()
        self.fil = parcop_sfil()
        self.gfil = parcop_gfil()

        # self.setPatch()

        # Communicators and proc boundaries
        self.master = False
        if self.comm.rank == 0:
            self.master = True

        self.x1proc = False
        if self.xcom.rank == 0:
            self.x1proc = True

        self.xnproc = False
        if self.xcom.rank == self.xcom.size - 1:
            self.xnproc = True

        self.y1proc = False
        if self.ycom.rank == 0:
            self.y1proc = True

        self.ynproc = False
        if self.ycom.rank == self.ycom.size - 1:
            self.ynproc = True

        self.z1proc = False
        if self.zcom.rank == 0:
            self.z1proc = True

        self.znproc = False
        if self.zcom.rank == self.zcom.size - 1:
            self.znproc = True

    def setPatch(self):
        parcop.parcop.set_patch(self.patch, self.level)

    def emptyScalar(self):
        return numpy.zeros(self.chunk_3d_size, dtype=numpy.float64, order="F") * 0.0

    def emptyVector(self, rank=3):
        blk_size = numpy.append(self.chunk_3d_size, rank)
        return numpy.zeros(blk_size, dtype=numpy.float64, order="F") * 0.0

    def gather2D(self, idata, com, n1, n2, g1, g2):
        """
        Give idata(a1,a2) return the assembled
          global array of gdata(n1,n2)
        """

        ldata = numpy.zeros((n1, n2))
        ldata[g1[0] : g1[1], g2[0] : g2[1]] = idata

        fldata = ldata.flatten()
        fldataT = fldata * 0.0
        com.Allreduce(fldata, fldataT, op=MPI.SUM)

        return numpy.reshape(fldataT, (n1, n2), order="C")

    def subsum3xz(self, data):
        icom = self.xzcom
        lsum = numpy.sum(data, (0, 2))
        gsum = self.xzcom.allreduce(lsum, op=MPI.SUM)
        Gsum = self.ycom.allgather(gsum)
        return numpy.concatenate(Gsum)

    def sum1D(self, data, com2, com1, indices):
        lsum = numpy.sum(data, indices)
        gsum = com2.allreduce(lsum, op=MPI.SUM)
        Gsum = com1.allgather(gsum)
        return numpy.concatenate(Gsum)

    def sum2D(self, data, com2, com1, index):
        lsum = numpy.sum(data, index)
        gsum = com1.allreduce(lsum, op=MPI.SUM)
        Gsum = com2.allgather(gsum)
        return numpy.concatenate(Gsum)

    def sum2Dold(self, data, com, n2, n3, index, g2, g3):

        a2 = g2[1] - g2[0]
        a3 = g3[1] - g3[0]
        gsum = numpy.zeros((a2, a3))

        # Get the local proc mean
        lsum = numpy.sum(data, index)

        tsum = self.gather2D(lsum, com, n2, n3, g2, g3)

        return tsum

    def sum3D(self, data):

        lsum = numpy.sum(data)
        tsum = self.comm.allreduce(lsum, op=MPI.SUM)

        return tsum

    def max3D(self, data):

        lmax = numpy.max(data)
        tmax = self.comm.allreduce(lmax, op=MPI.MAX)

        return tmax

    def min3D(self, data):

        lmin = numpy.min(data)
        tmin = self.comm.allreduce(lmin, op=MPI.MIN)

        return tmin

    # Deprecated functions
    def xbar(self, data):
        return self.xsum(data)

    def ybar(self, data):
        return self.ysum(data)

    def zbar(self, data):
        return self.zsum(data)

    # Reduction along 1 direction
    def xsum(self, data):
        return self.sum2D(data, self.yzcom, self.xcom, 0)

    def ysum(self, data):
        return self.sum2D(data, self.xzcom, self.ycom, 1)

    def zsum(self, data):
        return self.sum2D(data, self.xycom, self.zcom, 2)

    # Reduction along 2 directions
    def xysum(self, data):
        return self.sum1D(data, self.xycom, self.zcom, (0, 1))

    def xzsum(self, data):
        return self.sum1D(data, self.xzcom, self.ycom, (0, 2))

    def yzsum(self, data):
        return self.sum1D(data, self.yzcom, self.xcom, (1, 2))

    def mean3D(self, data, axis=None):
        """
        Compute the global mean of "data" and project it back onto the distributed
        3D data array
        """
        # Get the global mean
        tmpVar = pyrandaVar("tmp_mean3d_var", "mesh", "scalar")
        tmpVar.__allocate__(self)
        meanG = tmpVar.mean(axis=axis, idata=data)

        # Loop over local data and make it the mean
        if axis == None:
            tmpVar.data[:, :, :] = meanG
        else:
            if axis == 0:
                for k in range(self.az):
                    for j in range(self.ay):
                        jj = j + self.chunk_3d_lo[1]
                        kk = k + self.chunk_3d_lo[2]
                        tmpVar.data[:, j, k] = meanG[jj, kk]
            elif axis == 1:
                for k in range(self.az):
                    for i in range(self.ax):
                        ii = i + self.chunk_3d_lo[0]
                        kk = k + self.chunk_3d_lo[2]
                        tmpVar.data[i, :, k] = meanG[ii, kk]
            elif axis == 2:
                for j in range(self.ay):
                    for i in range(self.ax):
                        ii = i + self.chunk_3d_lo[0]
                        jj = j + self.chunk_3d_lo[1]
                        tmpVar.data[i, j, :] = meanG[ii, jj]
            elif axis == (0, 1):
                for k in range(self.az):
                    kk = k + self.chunk_3d_lo[2]
                    tmpVar.data[:, :, k] = meanG[kk]
            elif axis == (0, 2):
                for j in range(self.ay):
                    jj = j + self.chunk_3d_lo[1]
                    tmpVar.data[:, j, :] = meanG[jj]
            elif axis == (1, 2):
                for i in range(self.ax):
                    ii = i + self.chunk_3d_lo[0]
                    tmpVar.data[i, :, :] = meanG[ii]
            else:
                self.iprint("Error: unknown type given as axis")
                return None

        return tmpVar.data

    def iprint(self, sprnt):
        if self.master:
            print(sprnt)
            sys.stdout.flush()

    def getIJKpoo(self, data, I, J, K):
        """
        given [imin,imax, jmin,jmax, kmin,kmax] return the global data
        """
        gx = max(1, I[1] - I[0])
        gy = max(1, J[1] - J[0])
        gz = max(1, K[1] - K[0])
        gdata = numpy.zeros((gx, gy, gz))

        nx = data.shape[0]
        ny = data.shape[1]
        nz = data.shape[2]
        ax = max(1, self.chunk_3d_hi[0] - self.chunk_3d_lo[0] + 1)
        ay = max(1, self.chunk_3d_hi[1] - self.chunk_3d_lo[1] + 1)
        az = max(1, self.chunk_3d_hi[2] - self.chunk_3d_lo[2] + 1)
        for i in range(I[0], I[1]):
            for j in range(J[0], J[1]):
                for k in range(K[0], K[1]):
                    gi = i - I[0]
                    gj = j - J[0]
                    gk = k - K[0]
                    inI = (i >= self.chunk_3d_lo[0]) and (i <= self.chunk_3d_hi[0])
                    inJ = (j >= self.chunk_3d_lo[1]) and (j <= self.chunk_3d_hi[1])
                    inK = (k >= self.chunk_3d_lo[2]) and (k <= self.chunk_3d_hi[2])
                    if inI and inJ and inK:
                        gdata[gi, gj, gk] = data[i % ax, j % ay, k % az]

        gdata = self.comm.allreduce(gdata, op=MPI.SUM)

        return gdata

    def getIJK(self, data, I, J, K):
        """
        Same as readData but only reads in global range of data given by
        irange.
        """
        Rx = [0] * 2
        Ry = [0] * 2
        Rz = [0] * 2

        Rx[0] = I[0]
        Rx[1] = I[1]
        Ry[0] = J[0]
        Ry[1] = J[1]
        Rz[0] = K[0]
        Rz[1] = K[1]

        gx = max(1, I[1] - I[0])
        gy = max(1, J[1] - J[0])
        gz = max(1, K[1] - K[0])
        gdata = numpy.zeros((gx, gy, gz))

        g1 = self.chunk_3d_lo
        gn = self.chunk_3d_hi

        # Shift left point if node data
        iff = 0
        jff = 0
        kff = 0
        # try:
        c1 = Rx[1] in range(g1[0], gn[0] + 1)
        c2 = Rx[0] in range(g1[0], gn[0] + 1)
        c3 = (g1[0] and gn[0]) in range(Rx[0], Rx[1] + 1)
        CX = c1 or c2 or c3

        c1 = Ry[1] in range(g1[1], gn[1] + 1)
        c2 = Ry[0] in range(g1[1], gn[1] + 1)
        c3 = (g1[1] and gn[1]) in range(Ry[0], Ry[1] + 1)
        CY = c1 or c2 or c3

        c1 = Rz[1] in range(g1[2], gn[2] + 1)
        c2 = Rz[0] in range(g1[2], gn[2] + 1)
        c3 = (g1[2] and gn[2]) in range(Rz[0], Rz[1] + 1)
        CZ = c1 or c2 or c3

        if CX and CY and CZ:
            Li1 = numpy.max((0, Rx[0] - g1[0])) + iff
            Lif = numpy.min((Rx[1], gn[0] + 1)) - g1[0] + iff
            Ki1 = numpy.max((Rx[0], g1[0])) - Rx[0]
            Kif = Ki1 + (Lif - Li1)

            Lj1 = numpy.max((0, Ry[0] - g1[1])) + jff
            Ljf = numpy.min((Ry[1], gn[1] + 1)) - g1[1] + jff
            Kj1 = numpy.max((Ry[0], g1[1])) - Ry[0]
            Kjf = Kj1 + (Ljf - Lj1)

            Lk1 = numpy.max((0, Rz[0] - g1[2])) + kff
            Lkf = numpy.min((Rz[1], gn[2] + 1)) - g1[2] + kff
            Kk1 = numpy.max((Rz[0], g1[2])) - Rz[0]
            Kkf = Kk1 + (Lkf - Lk1)

            gdata[Ki1:Kif, Kj1:Kjf, Kk1:Kkf] = data[Li1:Lif, Lj1:Ljf, Lk1:Lkf]

        gdata = self.comm.allreduce(gdata, op=MPI.SUM)
        return gdata

    def ghost(self, data, np=1):

        bx = data.shape[0]
        by = data.shape[1]
        bz = data.shape[2]

        npx = npy = npz = 0
        if self.nx > 1:
            npx = np
        if self.ny > 1:
            npy = np
        if self.nz > 1:
            npz = np

        gdata = numpy.zeros((bx + 2 * npx, by + 2 * npy, bz + 2 * npz))
        gdata[npx : bx + npx, npy : by + npy, npz : bz + npz] = data

        if self.nx > 1:
            gdata = self.ghostx(gdata, npx)
        if self.ny > 1:
            gdata = self.ghosty(gdata, npy)
        if self.nz > 1:
            gdata = self.ghostz(gdata, npz)

        # For periodic data, clip sides
        if self.nx > 1:
            if self.x1proc:
                gdata = gdata[np:, :, :]
            if self.xnproc:
                gdata = gdata[:-np, :, :]

        if self.ny > 1:
            if self.y1proc:
                gdata = gdata[:, np:, :]
            if self.ynproc:
                gdata = gdata[:, :-np, :]

        if self.nz > 1:
            if self.z1proc:
                gdata = gdata[:, :, np:]
            if self.znproc:
                gdata = gdata[:, :, :-np]

        return gdata

    def ghostx(self, data, np):

        bx = data.shape[0]
        by = data.shape[1]
        bz = data.shape[2]

        i = numpy.ascontiguousarray(data[np : 2 * np, :, :])
        o = i * 0.0
        self.xcom.Sendrecv(i, self.xcom_lo, 0, o, self.xcom_hi, 0)
        data[bx - np : bx, :, :] = o

        i = numpy.ascontiguousarray(data[bx - 2 * np : bx - np, :, :])
        o = i * 0.0
        self.xcom.Sendrecv(i, self.xcom_hi, 1, o, self.xcom_lo, 1)
        data[0:np, :, :] = o

        if self.periodic[0]:
            return data

        if self.xcom.Get_rank() == 0:
            for i in range(np):
                data[i, :, :] = data[np, :, :]

        if self.xcom.Get_rank() == self.px - 1:
            for i in range(np):
                data[bx - np + i, :, :] = data[bx - np - 1, :, :]

        return data

    def ghosty(self, data, np):

        bx = data.shape[0]
        by = data.shape[1]
        bz = data.shape[2]

        i = numpy.ascontiguousarray(data[:, np : 2 * np, :])
        o = i * 0.0
        self.ycom.Sendrecv(i, self.ycom_lo, 0, o, self.ycom_hi, 0)
        data[:, by - np : by, :] = o

        i = numpy.ascontiguousarray(data[:, by - 2 * np : by - np, :])
        o = i * 0.0
        self.ycom.Sendrecv(i, self.ycom_hi, 1, o, self.ycom_lo, 1)
        data[:, 0:np, :] = o

        if self.periodic[1]:
            return data

        if self.ycom.Get_rank() == 0:
            for i in range(np):
                data[:, i, :] = data[:, np, :]

        if self.ycom.Get_rank() == self.py - 1:
            for i in range(np):
                data[:, by - np + i, :] = data[:, by - np - 1, :]

        return data

    def ghostz(self, data, np):

        bx = data.shape[0]
        by = data.shape[1]
        bz = data.shape[2]

        i = numpy.ascontiguousarray(data[:, :, np : 2 * np])
        o = i * 0.0
        self.zcom.Sendrecv(i, self.zcom_lo, 0, o, self.zcom_hi, 0)
        data[:, :, bz - np : bz] = o

        i = numpy.ascontiguousarray(data[:, :, bz - 2 * np : bz - np])
        o = i * 0.0
        self.zcom.Sendrecv(i, self.zcom_hi, 1, o, self.zcom_lo, 1)
        data[:, :, 0:np] = o

        if self.periodic[2]:
            return data

        if self.zcom.Get_rank() == 0:
            for i in range(np):
                data[:, :, i] = data[:, :, np]

        if self.zcom.Get_rank() == self.pz - 1:
            for i in range(np):
                data[:, :, bz - np + i] = data[:, :, bz - np - 1]

        return data

    def getVar(self, vname):
        return parcop.parcop.getvar(vname, self.ax, self.ay, self.az)


class parcop_der:
    def __init__(self):
        pass

    def ddx(self, val):
        return parcop.parcop.ddx(val)

    def ddy(self, val):
        return parcop.parcop.ddy(val)

    def ddz(self, val):
        return parcop.parcop.ddz(val)

    def dd4x(self, val):
        return parcop.parcop.dd4x(val)

    def dd4y(self, val):
        return parcop.parcop.dd4y(val)

    def dd4z(self, val):
        return parcop.parcop.dd4z(val)

    def dd8x(self, val):
        return parcop.parcop.dd8x(val)

    def dd8y(self, val):
        return parcop.parcop.dd8y(val)

    def dd8z(self, val):
        return parcop.parcop.dd8z(val)

    def div(self, fx, fy, fz):
        return parcop.parcop.divergence(fx, fy, fz)

    def divT(self, fxx, fxy, fxz, fyx, fyy, fyz, fzx, fzy, fzz):
        return parcop.parcop.divergencetensor(
            fxx, fxy, fxz, fyx, fyy, fyz, fzx, fzy, fzz
        )

    def grad(self, val):
        return parcop.parcop.grads(val)

    def laplacian(self, val):
        return parcop.parcop.plaplacian(val)

    def ring(self, val):
        return parcop.parcop.pring(val)

    def ringV(self, vx, vy, vz):
        return parcop.parcop.pringv(vx, vy, vz)


class parcop_mesh:
    def __init__(self):
        pass

    # def get


class parcop_gfil:
    def __init__(self):
        pass

    def filter(self, val):
        return parcop.parcop.gfilter(val)

    def filterDir(self, val, direction):
        return parcop.parcop.gfilterdir(val, direction)


class parcop_sfil:
    def __init__(self):
        pass

    def filter(self, val):
        return parcop.parcop.sfilter(val)


def factors(n):
    ff = set(
        reduce(
            list.__add__,
            ([i, n // i] for i in range(1, int(pow(n, 0.5) + 1)) if n % i == 0),
        )
    )
    return list(ff)


class pyDecomp:
    def __init__(self, xp, yp, zp, np, nx, ny, nz):

        self.xp = xp
        self.yp = yp
        self.zp = zp
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.np = np

        self.totalProcs = self.xp * self.yp * self.zp

        ax = nx / xp
        ay = ny / yp
        az = nz / zp

        surface = 2 * ax * ay + 2 * ay * az + 2 * az * ax
        volume = ax * ay * az

        self.surf2vol = float(surface) / float(volume)

        # optimal surface to volume ratio
        sopt = float(volume) ** (2.0 / 3.0) * 6.0 / float(volume)

        self.score = 1.0

        # Too many procs
        if self.totalProcs > self.np:
            self.score *= 0.0

        # Too small in any direction
        if (ax > 1) and (ax < 16):
            self.score *= 0.0

        if (ay > 1) and (ay < 16):
            self.score *= 0.0

        if (az > 1) and (az < 16):
            self.score *= 0.0

        # Total procs (under-utilization penalty
        # self.score *= self.totalProcs / self.np
        if self.totalProcs < self.np:
            self.score *= 0.0

        # Surface/Volume
        safrac = 0.1
        self.score -= safrac * (1.0 - sopt / self.surf2vol)


def decomp(nprocs, nx, ny, nz):
    """
    Return the optimal domain decomp for nprocs
    """

    # Get factors in each direction
    xfac = factors(nx)
    yfac = factors(ny)
    zfac = factors(nz)

    # Get all combos
    decomps = []
    for xf in xfac:
        for yf in yfac:
            for zf in zfac:
                xproc = nx / xf
                yproc = ny / yf
                zproc = nz / zf
                decomps.append(pyDecomp(xproc, yproc, zproc, nprocs, nx, ny, nz))

    decomps.sort(key=lambda x: x.score, reverse=True)
    return [decomps[0].xp, decomps[0].yp, decomps[0].zp]


def addCounter(counter1, counter2, datatype):
    for item in counter2:
        if item in counter1:
            counter1[item] += counter2[item]
        else:
            counter1[item] = counter2[item]

    return counter1
