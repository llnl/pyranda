# Copyright (c) 2018, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
#
# Written by: Britton J. Olson, olson45@llnl.gov, modified to use FFTs by Chris Scullard
################################################################################
import numpy
from .pyrandaPackage import pyrandaPackage


class pyrandaPoissonFFT(pyrandaPackage):
    def __init__(self, pysim):
        """
        Solve an elliptical Poisson equation
        \Delta \phi = f(x)
        Second order FD on structured mesh assumed
        """
        PackageName = "Poisson"

        # self.BCtype = BCtype
        self.pysim = pysim
        pyrandaPackage.__init__(self, PackageName, pysim)

        # Setup 1/k^2 array
        pi = numpy.pi
        self.ksqinv = numpy.ones((self.nx, self.ny, int(self.nz / 2 + 1)))
        for k in range(0, int(self.nz / 2)):
            for j in range(0, self.ny):
                for i in range(0, self.nx):
                    kx = 2.0 * pi * i
                    if i > (self.nx / 2):
                        kx = (i - self.nx) * 2.0 * pi
                    ky = 2.0 * pi * j
                    if j > (self.ny / 2):
                        ky = (j - self.ny) * 2.0 * pi
                    kz = 2.0 * pi * k
                    if k > (self.nz / 2):
                        kz = (k - self.nz) * 2.0 * pi
                    ksq = kx * kx + ky * ky + kz * kz
                    if ksq > 0:
                        # u_hat[i,j,k]=-u_hat[i,j,k]/ksq
                        self.ksqinv[i, j, k] = -1.0 / ksq

    def get_sMap(self):
        sMap = {}
        sMap["poisson.solve("] = "self.packages['Poisson'].solve("
        self.sMap = sMap

    def solve(self, rhs):

        # pi = numpy.pi
        # u = numpy.zeros((self.nx,self.ny,self.nz)) #RHS
        # uc = numpy.zeros((self.nx,self.ny,self.nz)) #RHS
        # u_hat = numpy.zeros((self.nx,self.ny,self.nz)) #RHS
        # for k in range(0,self.nz):
        #    for j in range(0,self.ny):
        #        for i in range(0, self.nx):
        #            u[i,j,k] = rhs[i,j,k]

        # u_hat = numpy.fft.rfftn(u) # forward FFT
        # import pdb
        # pdb.set_trace()

        u_hat = numpy.fft.rfftn(rhs)  # forward FFT

        # Divide by -k^2
        # for k in range(0,self.nz/2):
        #    for j in range(0,self.ny):
        #        for i in range(0, self.nx):
        #            kx=2.*pi*i
        #            if i>(self.nx/2):
        #                kx=(i-self.nx)*2.*pi
        #            ky=2.*pi*j
        #            if j>(self.ny/2):
        #                ky=(j-self.ny)*2.*pi
        #            kz=2.*pi*k
        #            if k>(self.nz/2):
        #                kz=(k-self.nz)*2.*pi
        #            ksq=kx*kx+ky*ky+kz*kz
        #            if ksq>0:
        #                u_hat[i,j,k]=-u_hat[i,j,k]/ksq

        u_hat = self.ksqinv * u_hat
        mysol = numpy.fft.irfftn(u_hat, s=(self.nx, self.ny, self.nz))  # inverse FFT

        # sol = self.solver( b )
        # sol = bicgstab(self.A, b )[0]
        # sol = cg(self.A, b )[0]

        # for j in range(0,self.ny):
        #   for i in range(0, self.nx):
        #      mysol[i,j,0] = sol[j + i*self.ny]

        # mysol -= mysol.mean()
        return mysol
