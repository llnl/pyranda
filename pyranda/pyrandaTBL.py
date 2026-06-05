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
from scipy import interpolate
from .pyrandaPackage import pyrandaPackage

from . import parcop


class pyrandaTBL(pyrandaPackage):
    def __init__(self, pysim):

        PackageName = "TBL"
        pyrandaPackage.__init__(self, PackageName, pysim)

        self.pysim = pysim
        self.bcList = {}

        self.BCdata = {}

        # Persistent 2D arrays (mesh size)
        self.A11 = None
        self.A22 = None
        self.A33 = None
        self.A12 = None
        self.Umean = None
        self.Tmean = None
        self.RHO_u = None
        self.RHO_v = None
        self.RHO_w = None
        self.buI = None
        self.bvI = None
        self.bwI = None
        self.buO = None
        self.bvO = None
        self.bwO = None
        self.y_r = None

        # U,V,W
        self.u = "u"
        self.v = "v"
        self.w = "w"

        # Global array
        self.y_G = None

        # Integers
        self.UIx = 0
        self.UIz = 0
        self.VIx = 0
        self.VIz = 0
        self.WIx = 0
        self.WIz = 0

        self.UIy = [0, 0]
        self.VIy = [0, 0]
        self.WIy = [0, 0]

        self.Nbuff = 0

        self.BL_data = None
        self.BL_delim = ","
        # Need boundary layer scales
        self.del_BL = 1.0  # Either BL thickness or y+ delta
        self.U_in = 1.0  # Some velocity scale... u_tau or or u_infty

    def get_sMap(self):
        sMap = {}
        sMap["TBL.inflow("] = "self.packages['TBL'].DFinflow("
        sMap["TBL.setup()"] = "self.packages['TBL'].setup()"
        self.sMap = sMap

    def setup(self):

        TBLfiles = self.BL_data

        # Flow in -x- wall in -y-
        aF = self.pysim.PyMPI.ax
        aW = self.pysim.PyMPI.ay
        aT = self.pysim.PyMPI.az
        nW = self.pysim.PyMPI.ny
        nT = self.pysim.PyMPI.nz

        self.aF = aF
        self.aW = aW
        self.aT = aT
        self.nW = nW
        self.nT = nT

        self.Umean = numpy.zeros((aW, aT))
        self.Tmean = numpy.zeros((aW, aT))
        self.A11 = numpy.zeros((aW, aT))
        self.A22 = numpy.zeros((aW, aT))
        self.A33 = numpy.zeros((aW, aT))
        self.A12 = numpy.zeros((aW, aT))

        self.RHO_u = numpy.zeros((aW, aT))
        self.RHO_v = numpy.zeros((aW, aT))
        self.RHO_w = numpy.zeros((aW, aT))

        self.y_G = numpy.zeros((nW))
        self.y_r = numpy.zeros((aW))

        tmp = self.pysim.PyMPI.emptyScalar()
        if self.pysim.PyMPI.x1proc:  # Flow dir
            tmp = self.pysim.mesh.coords[1].data * 1.0  # Wall dir
            tmp /= tmp.shape[0]
            tmp /= tmp.shape[2]
            tmp1 = numpy.sum(tmp, (0, 2))

            self.y_G = (
                numpy.concatenate(self.pysim.PyMPI.ycom.allgather(tmp1)) / self.del_BL
            )
            self.y_r = tmp1 / self.del_BL  # self.y_G[0:aW]

            UUmean = numpy.zeros((aW, aT))
            VVmean = numpy.zeros((aW, aT))
            WWmean = numpy.zeros((aW, aT))
            UVmean = numpy.zeros((aW, aT))

            # Read in BL data from file (non-dimensional data)
            #   and place on 2D planes Umean/Tmean/UUmean/VVmean/WWmean/UVmean
            fileName = TBLfiles["umean"]
            funcMean = readTBLdata(fileName, delimiter=self.BL_delim)
            yMean = funcMean(self.y_r)
            Uprof = funcMean(self.y_G) * self.U_in
            for ii in range(yMean.shape[0]):
                self.Umean[ii, :] = yMean[ii] * self.U_in

            fileName = TBLfiles["uumean"]
            funcMean = readTBLdata(fileName, delimiter=self.BL_delim)
            yMean = funcMean(self.y_r)
            for ii in range(yMean.shape[0]):
                UUmean[ii, :] = yMean[ii] * self.U_in

            fileName = TBLfiles["vvmean"]
            funcMean = readTBLdata(fileName, delimiter=self.BL_delim)
            yMean = funcMean(self.y_r)
            for ii in range(yMean.shape[0]):
                VVmean[ii, :] = yMean[ii] * self.U_in

            fileName = TBLfiles["wwmean"]
            funcMean = readTBLdata(fileName, delimiter=self.BL_delim)
            yMean = funcMean(self.y_r)
            for ii in range(yMean.shape[0]):
                WWmean[ii, :] = yMean[ii] * self.U_in

            fileName = TBLfiles["uvmean"]
            funcMean = readTBLdata(fileName, delimiter=self.BL_delim)
            yMean = funcMean(self.y_r)
            for ii in range(yMean.shape[0]):
                UVmean[ii, :] = yMean[ii] * self.U_in

            # Recast as componets of the Reynolds stress tensor
            self.A11 = UUmean
            self.A12 = numpy.where(self.A11 == 0.0, 0.0, UVmean**2 / self.A11)
            self.A22 = numpy.sqrt(VVmean**2 - self.A12**2)
            self.A33 = WWmean

            # Calculate the Displacement thickness
            del_star = 0.0
            for i in range(Uprof.shape[0] - 1):
                del_star = del_star + (Uprof[-1] - Uprof[i]) * (
                    self.y_G[i + 1] - self.y_G[i]
                )

            # Non-dimensionalize the local wall normal length scale for filter size selection
            self.y_r = self.y_r / del_star

            self.simtime_old = self.pysim.time
            # self.tauX = float(self.UIx)*del_star / self.U_in

            # Inner coefficients
            self.buI = setFiltCoeffs(self.UIy[0], self.UIz)
            self.bvI = setFiltCoeffs(self.VIy[0], self.VIz)
            self.bwI = setFiltCoeffs(self.WIy[0], self.WIz)

            # Outer coefficients
            self.buO = setFiltCoeffs(self.UIy[1], self.UIz)
            self.bvO = setFiltCoeffs(self.VIy[1], self.VIz)
            self.bwO = setFiltCoeffs(self.WIy[1], self.WIz)

        # END X1proc check

        # Set up a restart directory
        # // TODO //

        # Initialize the temporal averages here.  They will be
        #   over written on restart read. (done)

        # TODO: READ IN OLD DF FILES TO PICK UP THE RHO_u, RHO_v. RHO_w

        # DONE SETUP Digitial filtering

    def DFinflow(self):

        # Get a common seed
        seed = 0
        if self.pysim.PyMPI.master:
            seed = int(numpy.random.rand() * 2**16)

        seed = self.pysim.PyMPI.comm.allreduce(seed, op=MPI.SUM)

        if self.pysim.PyMPI.x1proc:
            # make parcop operator for these
            rands = parcop.parcop.tbl_get_rands(self.nW, self.nT, self.Nbuff, seed)

            iw1 = self.pysim.PyMPI.chunk_3d_lo[1]
            it1 = self.pysim.PyMPI.chunk_3d_lo[2]

            vU = parcop.parcop.tbl_filter(
                self.Nbuff,
                self.buI,
                self.buO,
                rands[0, :, :],
                self.nW,
                self.nT,
                self.aT,
                iw1,
                it1,
                self.y_r,
            )

            vV = parcop.parcop.tbl_filter(
                self.Nbuff,
                self.bvI,
                self.bvO,
                rands[1, :, :],
                self.nW,
                self.nT,
                self.aT,
                iw1,
                it1,
                self.y_r,
            )

            vW = parcop.parcop.tbl_filter(
                self.Nbuff,
                self.bwI,
                self.bwO,
                rands[2, :, :],
                self.nW,
                self.nT,
                self.aT,
                iw1,
                it1,
                self.y_r,
            )

            # Restart IO

            # Get the updated rho_k
            # Time avergaging coefficients and quantities/fluctuations
            # Need to write a restart file with averages stored
            t_nm1 = self.simtime_old * 1.0
            t_n = self.pysim.time * 1.0
            self.simtime_old = t_n * 1.0
            dt_n = t_n - t_nm1
            EXPt = numpy.exp(-numpy.pi * dt_n / self.tauX)
            wgt1 = numpy.sqrt(EXPt)
            wgt2 = numpy.sqrt(1.0 - EXPt)
            self.RHO_u = self.RHO_u * wgt1 + vU * wgt2
            self.RHO_v = self.RHO_v * wgt1 + vV * wgt2
            self.RHO_w = self.RHO_w * wgt1 + vW * wgt2

            # Add the perturbations to mean with given 2 point correlations
            uinlet = self.Umean + self.A11 * self.RHO_u
            vinlet = self.A12 * self.RHO_u + self.A22 * self.RHO_v
            winlet = self.A33 * self.RHO_w

            # Get the temperature perturbation using strong Reynolds Analogy (SRA)
            # MaSq = Mach**2
            # Tinlet = self.Tmean + self.Tmean * ( -gm1*MaSq * (uinlet - self.Umean) / U_in )

            # Add to inlet

            self.pysim.variables[self.u].data[0, :, :] = uinlet
            self.pysim.variables[self.v].data[0, :, :] = vinlet
            self.pysim.variables[self.w].data[0, :, :] = winlet
            # self.pysim.variables[self.T].data[0,:,:] = Tinlet

        # Reconcile EOS outside of this loop/BC call

    def get_rands(self, ny, nz, Nbuff):

        # Get a common seed
        seed = 0
        if self.pysim.PyMPI.master:
            seed = int(numpy.random.rand() * 2**32)
            seed = self.pysim.PyMPI.comm.allreduce(seed, op=MPI.SUM)

        numpy.random.seed(seed)

        return numpy.random.rand(4, ny + Nbuff, nz + Nbuff)


def setFiltCoeffs(N1, N2):

    bij = numpy.zeros((2 * N1 + 1, 2 * N2 + 1))

    den1 = calcDen(N1)
    den2 = calcDen(N2)

    N1d2 = float(N1) / 2.0
    N2d2 = float(N2) / 2.0

    for i in range(1, 2 * N1 + 2):
        for j in range(1, 2 * N2 + 2):
            k1 = abs(i - (N1 + 1))
            k2 = abs(j - (N2 + 1))
            bij[i - 1, j - 1] = (
                numpy.exp(-numpy.pi * float(k1) / N1d2)
                / den1
                * numpy.exp(-numpy.pi * float(k2) / N2d2)
                / den2
            )

    return bij


def calcDen(N):

    den = 0.0
    for i in range(-N, N + 1):  # i=-N,N
        den += (numpy.exp(-numpy.pi * abs(float(i)) / (float(N) / 2.0))) ** 2.0

    return numpy.sqrt(den)


def readTBLdata(fileName, delimiter=None):
    """
    Give a file name, return a numpy.interp object to evalaute to give the value
    """

    ff = numpy.loadtxt(fileName, delimiter=delimiter)
    ii = interpolate.interp1d(ff[:, 0], ff[:, 1], fill_value="extrapolate")

    return ii
