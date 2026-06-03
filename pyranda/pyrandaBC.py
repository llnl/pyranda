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
from .pyrandaPackage import pyrandaPackage


class pyrandaBC(pyrandaPackage):
    def __init__(self, pysim):

        PackageName = "BC"
        pyrandaPackage.__init__(self, PackageName, pysim)

        self.bcList = {}

        self.BCdata = {}

    def get_sMap(self):
        sMap = {}
        sMap["bc.extrap("] = "self.packages['BC'].extrap("
        sMap["bc.const("] = "self.packages['BC'].const("
        sMap["bc.exit("] = "self.packages['BC'].exitbc("
        sMap["bc.slip("] = "self.packages['BC'].slipbc("
        sMap["bc.field("] = "self.packages['BC'].applyBC('field',"
        sMap["bc.farfield("] = "self.packages['BC'].farfieldbc("
        sMap["bc.symm("] = "self.packages['BC'].symm("
        self.sMap = sMap

    def applyBC(self, func, var, direction, val=None):

        pFunc = eval("self." + func)

        if type(var) is not type([]):
            var = [var]

        if type(direction) is not type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                if type(val) is not type(None):
                    pFunc(v, d, val)
                else:
                    pFunc(v, d)

    def extrap(self, var, direction, order=2):

        if type(var) is not type([]):
            var = [var]

        if type(direction) is not type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.extrapolate(v, d, order)

    def extrapolate(self, var, direction, order):
        # Direction switch
        if direction == "x1":
            if self.pyranda.PyMPI.x1proc:
                if order == 2:
                    self.pyranda.variables[var].data[0, :, :] = (
                        2 * self.pyranda.variables[var].data[1, :, :]
                        - self.pyranda.variables[var].data[2, :, :]
                    )
                else:
                    self.pyranda.variables[var].data[0, :, :] = self.pyranda.variables[
                        var
                    ].data[1, :, :]

        if direction == "xn":
            if self.pyranda.PyMPI.xnproc:
                if order == 2:
                    self.pyranda.variables[var].data[-1, :, :] = (
                        2 * self.pyranda.variables[var].data[-2, :, :]
                        - self.pyranda.variables[var].data[-3, :, :]
                    )
                else:
                    self.pyranda.variables[var].data[-1, :, :] = self.pyranda.variables[
                        var
                    ].data[-2, :, :]

        if direction == "y1":
            if self.pyranda.PyMPI.y1proc:
                if order == 2:
                    self.pyranda.variables[var].data[:, 0, :] = (
                        2 * self.pyranda.variables[var].data[:, 1, :]
                        - self.pyranda.variables[var].data[:, 2, :]
                    )
                else:
                    self.pyranda.variables[var].data[:, 0, :] = self.pyranda.variables[
                        var
                    ].data[:, 1, :]

        if direction == "yn":
            if self.pyranda.PyMPI.ynproc:
                if order == 2:
                    self.pyranda.variables[var].data[:, -1, :] = (
                        2 * self.pyranda.variables[var].data[:, -2, :]
                        - self.pyranda.variables[var].data[:, -3, :]
                    )
                else:
                    self.pyranda.variables[var].data[:, -1, :] = self.pyranda.variables[
                        var
                    ].data[:, -2, :]

        if direction == "z1":
            if self.pyranda.PyMPI.z1proc:
                if order == 2:
                    self.pyranda.variables[var].data[:, :, 0] = (
                        2 * self.pyranda.variables[var].data[:, :, 1]
                        - self.pyranda.variables[var].data[:, :, 2]
                    )
                else:
                    self.pyranda.variables[var].data[:, :, 0] = self.pyranda.variables[
                        var
                    ].data[:, :, 1]

        if direction == "zn":
            if self.pyranda.PyMPI.znproc:
                if order == 2:
                    self.pyranda.variables[var].data[:, :, -1] = (
                        2 * self.pyranda.variables[var].data[:, :, -2]
                        - self.pyranda.variables[var].data[:, :, -3]
                    )
                else:
                    self.pyranda.variables[var].data[:, :, -1] = self.pyranda.variables[
                        var
                    ].data[:, :, -2]

    def const(self, var, direction, val):

        if type(var) is not type([]):
            var = [var]

        if type(direction) is not type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.constant(v, d, val)

    def constant(self, var, direction, val):

        # Direction switch
        if direction == "x1":
            if self.pyranda.PyMPI.x1proc:
                self.pyranda.variables[var].data[0, :, :] = val

        if direction == "xn":
            if self.pyranda.PyMPI.xnproc:
                self.pyranda.variables[var].data[-1, :, :] = val

        if direction == "y1":
            if self.pyranda.PyMPI.y1proc:
                self.pyranda.variables[var].data[:, 0, :] = val

        if direction == "yn":
            if self.pyranda.PyMPI.ynproc:
                self.pyranda.variables[var].data[:, -1, :] = val

        if direction == "z1":
            if self.pyranda.PyMPI.z1proc:
                self.pyranda.variables[var].data[:, :, 0] = val

        if direction == "zn":
            if self.pyranda.PyMPI.znproc:
                self.pyranda.variables[var].data[:, :, -1] = val

    def field(self, var, direction, field):

        val = field

        # Direction switch
        if direction == "x1":
            if self.pyranda.PyMPI.x1proc:
                self.pyranda.variables[var].data[0, :, :] = val[0, :, :]

        if direction == "xn":
            if self.pyranda.PyMPI.xnproc:
                self.pyranda.variables[var].data[-1, :, :] = val[-1, :, :]

        if direction == "y1":
            if self.pyranda.PyMPI.y1proc:
                self.pyranda.variables[var].data[:, 0, :] = val[:, 0, :]

        if direction == "yn":
            if self.pyranda.PyMPI.ynproc:
                self.pyranda.variables[var].data[:, -1, :] = val[:, -1, :]

    def slipbc(self, var, direction):

        if type(var) is not type([]):
            var = [var]

        if type(direction) is not type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.slipvelbc(v, d)

    def getMetrics(self):
        """
        Compute the metrics inverse
        """
        dAdx = self.pyranda.getVar("dAx")
        dAdy = self.pyranda.getVar("dAy")
        dAdz = self.pyranda.getVar("dAz")
        dBdx = self.pyranda.getVar("dBx")
        dBdy = self.pyranda.getVar("dBy")
        dBdz = self.pyranda.getVar("dBz")
        dCdx = self.pyranda.getVar("dCx")
        dCdy = self.pyranda.getVar("dCy")
        dCdz = self.pyranda.getVar("dCz")
        detJ = self.pyranda.getVar("dtJ")

        metrics = {}

        metrics["dxdA"] = (-dBdz * dCdy + dBdy * dCdz) * detJ
        metrics["dxdB"] = (dAdz * dCdy - dAdy * dCdz) * detJ
        metrics["dxdC"] = (-dAdz * dBdy + dAdy * dBdz) * detJ

        metrics["dydA"] = (dBdz * dCdx - dBdx * dCdz) * detJ
        metrics["dydB"] = (-dAdz * dCdx + dAdx * dCdz) * detJ
        metrics["dydC"] = (dAdz * dBdx - dAdx * dBdz) * detJ

        metrics["dzdA"] = (-dBdy * dCdx + dBdx * dCdy) * detJ
        metrics["dzdB"] = (dAdy * dCdx - dAdx * dCdy) * detJ
        metrics["dzdC"] = (-dAdy * dBdx + dAdx * dBdy) * detJ

        return metrics

    def getNorms(self, boundary):

        metrics = self.getMetrics()

        if "x" in boundary:
            d1 = "B"
            d2 = "C"
        elif "y" in boundary:
            d1 = "A"
            d2 = "C"
        else:
            d1 = "A"
            d2 = "B"

        a1 = metrics["dxd" + d1]
        a2 = metrics["dyd" + d1]
        a3 = metrics["dzd" + d1]

        b1 = metrics["dxd" + d2]
        b2 = metrics["dyd" + d2]
        b3 = metrics["dzd" + d2]

        # Get normal vector (via cross product)
        n1 = a2 * b3 - a3 * b2
        n2 = -(a1 * b3 - a3 * b1)
        n3 = a1 * b2 - a2 * b1
        mag = numpy.sqrt(n1 * n1 + n2 * n2 + n3 * n3)
        n1 /= mag
        n2 /= mag
        n3 /= mag

        return [n1, n2, n3]

    def slipvelbc(self, velocity, direction):
        """
        Allow tangential velocities
        """

        if direction == "y1":
            if "slipbc-y1" not in self.BCdata.keys():
                [n1, n2, n3] = self.getNorms(direction)
                norms = []
                norms.append(n1[:, 0, :])
                norms.append(n2[:, 0, :])
                norms.append(n3[:, 0, :])
                self.BCdata["slipbc-y1"] = norms
            else:
                norms = self.BCdata["slipbc-y1"]

            # For free slip, always extrapolate first
            for uu in velocity:
                self.extrapolate(uu, "y1", order=2)

            if self.pyranda.PyMPI.y1proc:
                udotn = 0.0
                for uu in range(len(velocity)):
                    u = self.pyranda.variables[velocity[uu]].data[:, 0, :]
                    udotn = udotn + u * norms[uu]

                # magnitude before
                mag0 = 0.0
                for uu in range(len(velocity)):
                    mag0 += self.pyranda.variables[velocity[uu]].data[:, 0, :] ** 2
                mag0 = numpy.sqrt(mag0)  # + 1.0e-12

                for uu in range(len(velocity)):
                    self.pyranda.variables[velocity[uu]].data[:, 0, :] -= (
                        udotn * norms[uu]
                    )

                # magnitude after
                magF = 0.0
                for uu in range(len(velocity)):
                    magF += self.pyranda.variables[velocity[uu]].data[:, 0, :] ** 2
                magF = numpy.sqrt(magF)  # + 1.0e-12

                # Make sure magnitude is NOT larger
                for uu in range(len(velocity)):
                    self.pyranda.variables[velocity[uu]].data[:, 0, :] = numpy.where(
                        magF > mag0,
                        self.pyranda.variables[velocity[uu]].data[:, 0, :]
                        * mag0
                        / magF,
                        self.pyranda.variables[velocity[uu]].data[:, 0, :],
                    )

        if direction == "yn":
            if "slipbc-yn" not in self.BCdata.keys():
                [n1, n2, n3] = self.getNorms(direction)
                norms = []
                norms.append(n1[:, -1, :])
                norms.append(n2[:, -1, :])
                norms.append(n3[:, -1, :])
                self.BCdata["slipbc-yn"] = norms

            else:
                norms = self.BCdata["slipbc-yn"]

            # For free slip, always extrapolate first
            for uu in velocity:
                self.extrapolate(uu, "yn", order=2)

            if self.pyranda.PyMPI.ynproc:
                udotn = 0.0
                for uu in range(len(velocity)):
                    u = self.pyranda.variables[velocity[uu]].data[:, -1, :]
                    udotn = udotn + u * norms[uu]

                # magnitude before
                mag0 = 0.0
                for uu in range(len(velocity)):
                    mag0 += self.pyranda.variables[velocity[uu]].data[:, -1, :] ** 2
                mag0 = numpy.sqrt(mag0)  # + 1.0e-12

                for uu in range(len(velocity)):
                    self.pyranda.variables[velocity[uu]].data[:, -1, :] -= (
                        udotn * norms[uu]
                    )

                # magnitude after
                magF = 0.0
                for uu in range(len(velocity)):
                    magF += self.pyranda.variables[velocity[uu]].data[:, -1, :] ** 2
                magF = numpy.sqrt(magF)  # + 1.0e-12

                # Make sure magnitude is NOT larger
                for uu in range(len(velocity)):
                    self.pyranda.variables[velocity[uu]].data[:, -1, :] = numpy.where(
                        magF > mag0,
                        self.pyranda.variables[velocity[uu]].data[:, -1, :]
                        * mag0
                        / magF,
                        self.pyranda.variables[velocity[uu]].data[:, -1, :],
                    )

        if direction == "x1":
            if "slipbc-x1" not in self.BCdata.keys():
                [n1, n2, n3] = self.getNorms(direction)
                norms = []
                norms.append(n1[0, :, :])
                norms.append(n2[0, :, :])
                norms.append(n3[0, :, :])
                self.BCdata["slipbc-x1"] = norms
            else:
                norms = self.BCdata["slipbc-x1"]

            # For free slip, always extrapolate first
            for uu in velocity:
                self.extrapolate(uu, "x1", order=2)

            if self.pyranda.PyMPI.x1proc:
                udotn = 0.0
                for uu in range(len(velocity)):
                    u = self.pyranda.variables[velocity[uu]].data[0, :, :]
                    udotn = udotn + u * norms[uu]

                # magnitude before
                mag0 = 0.0
                for uu in range(len(velocity)):
                    mag0 += self.pyranda.variables[velocity[uu]].data[0, :, :] ** 2
                mag0 = numpy.sqrt(mag0)  # + 1.0e-12

                for uu in range(len(velocity)):
                    self.pyranda.variables[velocity[uu]].data[0, :, :] -= (
                        udotn * norms[uu]
                    )

                # magnitude after
                magF = 0.0
                for uu in range(len(velocity)):
                    magF += self.pyranda.variables[velocity[uu]].data[0, :, :] ** 2
                magF = numpy.sqrt(magF)  # + 1.0e-12

                # Make sure magnitude is NOT larger
                for uu in range(len(velocity)):
                    self.pyranda.variables[velocity[uu]].data[0, :, :] = numpy.where(
                        magF > mag0,
                        self.pyranda.variables[velocity[uu]].data[0, :, :]
                        * mag0
                        / magF,
                        self.pyranda.variables[velocity[uu]].data[0, :, :],
                    )

        if direction == "xn":
            if "slipbc-xn" not in self.BCdata.keys():
                [n1, n2, n3] = self.getNorms(direction)
                norms = []
                norms.append(n1[-1, :, :])
                norms.append(n2[-1, :, :])
                norms.append(n3[-1, :, :])
                self.BCdata["slipbc-xn"] = norms

            else:
                norms = self.BCdata["slipbc-xn"]

            # For free slip, always extrapolate first
            for uu in velocity:
                self.extrapolate(uu, "xn", order=2)

            if self.pyranda.PyMPI.xnproc:
                udotn = 0.0
                for uu in range(len(velocity)):
                    u = self.pyranda.variables[velocity[uu]].data[-1, :, :]
                    udotn = udotn + u * norms[uu]

                # magnitude before
                mag0 = 0.0
                for uu in range(len(velocity)):
                    mag0 += self.pyranda.variables[velocity[uu]].data[-1, :, :] ** 2
                mag0 = numpy.sqrt(mag0)  # + 1.0e-12

                for uu in range(len(velocity)):
                    self.pyranda.variables[velocity[uu]].data[-1, :, :] -= (
                        udotn * norms[uu]
                    )

                # magnitude after
                magF = 0.0
                for uu in range(len(velocity)):
                    magF += self.pyranda.variables[velocity[uu]].data[-1, :, :] ** 2
                magF = numpy.sqrt(magF)  # + 1.0e-12

                # Make sure magnitude is NOT larger
                for uu in range(len(velocity)):
                    self.pyranda.variables[velocity[uu]].data[-1, :, :] = numpy.where(
                        magF > mag0,
                        self.pyranda.variables[velocity[uu]].data[-1, :, :]
                        * mag0
                        / magF,
                        self.pyranda.variables[velocity[uu]].data[-1, :, :],
                    )

    def exitbc(self, var, direction, norm=False):

        if type(var) is not type([]):
            var = [var]

        if type(direction) is not type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.exit_boundary(v, d, norm)

    def exit_boundary(self, var, direction, norm):

        if direction == "x1":
            if self.pyranda.PyMPI.x1proc:
                u1 = self.pyranda.variables[var].data[0, :, :]
                u2 = self.pyranda.variables[var].data[1, :, :]
                u3 = self.pyranda.variables[var].data[2, :, :]
                self.pyranda.variables[var].data[0, :, :] = self.BENO(u1, u2, u3, norm)
        if direction == "xn":
            if self.pyranda.PyMPI.xnproc:
                u1 = self.pyranda.variables[var].data[-1, :, :]
                u2 = self.pyranda.variables[var].data[-2, :, :]
                u3 = self.pyranda.variables[var].data[-3, :, :]
                self.pyranda.variables[var].data[-1, :, :] = self.BENO(u1, u2, u3, norm)
        if direction == "y1":
            if self.pyranda.PyMPI.y1proc:
                u1 = self.pyranda.variables[var].data[:, 0, :]
                u2 = self.pyranda.variables[var].data[:, 1, :]
                u3 = self.pyranda.variables[var].data[:, 2, :]
                self.pyranda.variables[var].data[:, 0, :] = self.BENO(u1, u2, u3, norm)
        if direction == "yn":
            if self.pyranda.PyMPI.ynproc:
                u1 = self.pyranda.variables[var].data[:, -1, :]
                u2 = self.pyranda.variables[var].data[:, -2, :]
                u3 = self.pyranda.variables[var].data[:, -3, :]
                self.pyranda.variables[var].data[:, -1, :] = self.BENO(u1, u2, u3, norm)

    def BENO(self, u1, u2, u3, norm):
        d1 = u2 - u1
        d2 = (2.0 * u2 - u3) - u1

        beno = numpy.where(numpy.abs(d1) <= numpy.abs(d2), u1 + d1, u1 + d2)
        beno = numpy.where(d1 * d2 <= 0.0, u1, beno)

        if norm:
            beno = numpy.maximum(
                numpy.minimum(beno, numpy.maximum(0.0, u2)), numpy.minimum(0.0, u2)
            )

        return beno

    def getnormal(self, direction):

        [n1, n2, n3] = self.getNorms(direction)

        dataKey = "farfield-%s" % direction
        if dataKey in self.BCdata.keys():
            [xnormal, ynormal, znormal] = self.BCdata[dataKey]

        else:
            if direction == "x1":
                xnormal = n1[0, :, :]
                ynormal = n2[0, :, :]
                znormal = n3[0, :, :]
            if direction == "y1":
                xnormal = n1[:, 0, :]
                ynormal = n2[:, 0, :]
                znormal = n3[:, 0, :]
            if direction == "z1":
                xnormal = n1[:, :, 0]
                ynormal = n2[:, :, 0]
                znormal = n3[:, :, 0]
            if direction == "xn":
                xnormal = n1[-1, :, :]
                ynormal = n2[-1, :, :]
                znormal = n3[-1, :, :]
            if direction == "yn":
                xnormal = n1[:, -1, :]
                ynormal = n2[:, -1, :]
                znormal = n3[:, -1, :]
            if direction == "zn":
                xnormal = n1[:, :, -1]
                ynormal = n2[:, :, -1]
                znormal = n3[:, :, -1]

            self.BCdata[dataKey] = [xnormal, ynormal, znormal]

        return xnormal, ynormal, znormal

    def Reimann(
        self, direction, Rhoref, Uref, Vref, Wref, Pref, gamma, rhoi, Ui, Vi, Wi, Pi
    ):
        """
        Riemann invariants
        index
        b - variables at the boundary index=0,ax,ay,az
        i - variables inside the domain index=1,n-2,...
        o - variables at outiside of the domain (inf)
        """

        # itializing free stream vectors
        rhoo = numpy.zeros((Ui.shape[0], 1)) + Rhoref
        Uo = numpy.zeros((Ui.shape[0], 1)) + Uref
        Vo = numpy.zeros((Ui.shape[0], 1)) + Vref
        Wo = numpy.zeros((Ui.shape[0], 1)) + Wref
        Po = numpy.zeros((Ui.shape[0], 1)) + Pref

        # getting normal vectors
        nx, ny, nz = self.getnormal(direction)

        Vni = Ui * nx + Vi * ny + Wi * nz
        Vno = Uo * nx + Vo * ny + Wo * nz
        SSo = numpy.sqrt(gamma * Po / rhoo)
        SSi = numpy.sqrt(gamma * Pi / rhoi)
        tmp = numpy.sqrt(Ui**2 + Vi**2 + Wi**2) / SSi

        # Compute Riemann invariants for supersonic and subsonic cases
        Rplus = numpy.where(
            tmp >= 1.0, Vno + 2.0 * SSo / (gamma - 1.0), Vni + 2.0 * SSi / (gamma - 1.0)
        )
        Rminus = numpy.where(
            tmp >= 1.0, Vni - 2.0 * SSi / (gamma - 1.0), Vno - 2.0 * SSo / (gamma - 1.0)
        )

        Vnormal = (Rminus + Rplus) * 0.5
        SSb = (Rplus - Rminus) * (gamma - 1.0) * 0.25
        Ub = numpy.where(
            Vnormal >= 0.0, Ui + (Vnormal - Vni) * nx, Uo + (Vnormal - Vno) * nx
        )
        Vb = numpy.where(
            Vnormal >= 0.0, Vi + (Vnormal - Vni) * ny, Vo + (Vnormal - Vno) * ny
        )
        Wb = numpy.where(
            Vnormal >= 0.0, Wi + (Vnormal - Vni) * nz, Wo + (Vnormal - Vno) * nz
        )

        aux = 1.0 / (gamma - 1.0)
        rhob = numpy.where(
            Vnormal >= 0.0,
            ((rhoi**gamma * SSb**2) / (gamma * Pi)) ** aux,
            ((rhoo**gamma * SSb**2) / (gamma * Po)) ** aux,
        )

        Pb = rhob * SSb**2 / gamma
        Etb = Pb / (gamma - 1.0) + 0.5 * rhob * (Ub * Ub + Vb * Vb + Wb * Wb)

        return rhob, Ub, Vb, Wb, Etb, Pb

    def farfieldbc(self, direction):

        if type(direction) is not type([]):
            direction = [direction]

        for d in direction:
            self.farfield(d)

    def farfield(self, direction):

        refdata = self.BCdata["farfield-properties-%s" % direction]
        Rhoref = refdata["rho0"]
        Uref = refdata["u0"]
        Vref = refdata["v0"]
        Wref = refdata["w0"]
        Pref = refdata["p0"]
        gamma = refdata["gamma"]

        rho = refdata["rho"]
        u = refdata["u"]
        v = refdata["v"]
        w = refdata["w"]
        p = refdata["p"]

        # Direction switch
        if direction == "x1":
            if self.pyranda.PyMPI.x1proc:
                rhoi = self.pyranda.variables[rho].data[1, :, :]
                Ui = self.pyranda.variables[u].data[1, :, :]
                Vi = self.pyranda.variables[v].data[1, :, :]
                Wi = self.pyranda.variables[w].data[1, :, :]
                Pi = self.pyranda.variables[p].data[1, :, :]

                [rhob, Ub, Vb, Wb, Etb, Pb] = self.Reimann(
                    direction,
                    Rhoref,
                    Uref,
                    Vref,
                    Wref,
                    Pref,
                    gamma,
                    rhoi,
                    Ui,
                    Vi,
                    Wi,
                    Pi,
                )

                # update variables
                self.pyranda.variables[rho].data[0, :, :] = rhob
                self.pyranda.variables[u].data[0, :, :] = Ub
                self.pyranda.variables[v].data[0, :, :] = Vb
                self.pyranda.variables[w].data[0, :, :] = Wb
                self.pyranda.variables[p].data[0, :, :] = Pb

        if direction == "xn":
            if self.pyranda.PyMPI.xnproc:
                rhoi = self.pyranda.variables[rho].data[-2, :, :]
                Ui = self.pyranda.variables[u].data[-2, :, :]
                Vi = self.pyranda.variables[v].data[-2, :, :]
                Wi = self.pyranda.variables[w].data[-2, :, :]
                Pi = self.pyranda.variables[p].data[-2, :, :]

                [rhob, Ub, Vb, Wb, Etb, Pb] = self.Reimann(
                    direction,
                    Rhoref,
                    Uref,
                    Vref,
                    Wref,
                    Pref,
                    gamma,
                    rhoi,
                    Ui,
                    Vi,
                    Wi,
                    Pi,
                )

                # update variables
                self.pyranda.variables[rho].data[-1, :, :] = rhob
                self.pyranda.variables[u].data[-1, :, :] = Ub
                self.pyranda.variables[v].data[-1, :, :] = Vb
                self.pyranda.variables[w].data[-1, :, :] = Wb
                self.pyranda.variables[p].data[-1, :, :] = Pb

        if direction == "y1":
            if self.pyranda.PyMPI.y1proc:
                rhoi = self.pyranda.variables[rho].data[:, 1, :]
                Ui = self.pyranda.variables[u].data[:, 1, :]
                Vi = self.pyranda.variables[v].data[:, 1, :]
                Wi = self.pyranda.variables[w].data[:, 1, :]
                Pi = self.pyranda.variables[p].data[:, 1, :]

                [rhob, Ub, Vb, Wb, Etb, Pb] = self.Reimann(
                    direction,
                    Rhoref,
                    Uref,
                    Vref,
                    Wref,
                    Pref,
                    gamma,
                    rhoi,
                    Ui,
                    Vi,
                    Wi,
                    Pi,
                )

                # update variables
                self.pyranda.variables[rho].data[:, 0, :] = rhob
                self.pyranda.variables[u].data[:, 0, :] = Ub
                self.pyranda.variables[v].data[:, 0, :] = Vb
                self.pyranda.variables[w].data[:, 0, :] = Wb
                self.pyranda.variables[p].data[:, 0, :] = Pb

        if direction == "yn":
            if self.pyranda.PyMPI.ynproc:
                rhoi = self.pyranda.variables[rho].data[:, -2, :]
                Ui = self.pyranda.variables[u].data[:, -2, :]
                Vi = self.pyranda.variables[v].data[:, -2, :]
                Wi = self.pyranda.variables[w].data[:, -2, :]
                Pi = self.pyranda.variables[p].data[:, -2, :]

                [rhob, Ub, Vb, Wb, Etb, Pb] = self.Reimann(
                    direction,
                    Rhoref,
                    Uref,
                    Vref,
                    Wref,
                    Pref,
                    gamma,
                    rhoi,
                    Ui,
                    Vi,
                    Wi,
                    Pi,
                )

                # update variables
                self.pyranda.variables[rho].data[:, -1, :] = rhob
                self.pyranda.variables[u].data[:, -1, :] = Ub
                self.pyranda.variables[v].data[:, -1, :] = Vb
                self.pyranda.variables[w].data[:, -1, :] = Wb
                self.pyranda.variables[p].data[:, -1, :] = Pb

        if direction == "z1":
            if self.pyranda.PyMPI.z1proc:
                rhoi = self.pyranda.variables[rho].data[:, :, 1]
                Ui = self.pyranda.variables[u].data[:, :, 1]
                Vi = self.pyranda.variables[v].data[:, :, 1]
                Wi = self.pyranda.variables[w].data[:, :, 1]
                Pi = self.pyranda.variables[p].data[:, :, 1]

                [rhob, Ub, Vb, Wb, Etb, Pb] = self.Reimann(
                    direction,
                    Rhoref,
                    Uref,
                    Vref,
                    Wref,
                    Pref,
                    gamma,
                    rhoi,
                    Ui,
                    Vi,
                    Wi,
                    Pi,
                )

                # update variables
                self.pyranda.variables[rho].data[:, :, 0] = rhob
                self.pyranda.variables[u].data[:, :, 0] = Ub
                self.pyranda.variables[v].data[:, :, 0] = Vb
                self.pyranda.variables[w].data[:, :, 0] = Wb
                self.pyranda.variables[p].data[:, :, 0] = Pb

        if direction == "zn":
            if self.pyranda.PyMPI.znproc:
                rhoi = self.pyranda.variables[rho].data[:, :, -2]
                Ui = self.pyranda.variables[u].data[:, :, -2]
                Vi = self.pyranda.variables[v].data[:, :, -2]
                Wi = self.pyranda.variables[w].data[:, :, -2]
                Pi = self.pyranda.variables[p].data[:, :, -2]

                [rhob, Ub, Vb, Wb, Etb, Pb] = self.Reimann(
                    direction,
                    Rhoref,
                    Uref,
                    Vref,
                    Wref,
                    Pref,
                    gamma,
                    rhoi,
                    Ui,
                    Vi,
                    Wi,
                    Pi,
                )

                # update variables
                self.pyranda.variables[rho].data[:, :, -1] = rhob
                self.pyranda.variables[u].data[:, :, -1] = Ub
                self.pyranda.variables[v].data[:, :, -1] = Vb
                self.pyranda.variables[w].data[:, :, -1] = Wb
                self.pyranda.variables[p].data[:, :, -1] = Pb

    def symm(self, var, direction, anti=False, npts=4):

        if type(var) is not type([]):
            var = [var]

        # For anti-symmetric variables
        antiM = 1.0
        if anti:
            antiM = -1.0

        if type(direction) is not type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.symmetric(v, d, antiM, npts)

    def symmetric(self, var, direction, antiMult, npts):
        # Direction switch
        varpt = self.pyranda.variables[var].data

        if direction == "x1":
            if self.pyranda.PyMPI.x1proc:
                varpt[0:npts, :, :] = (
                    varpt[npts : npts * 2, :, :][::-1, :, :] * antiMult
                )
        if direction == "xn":
            if self.pyranda.PyMPI.xnproc:
                varpt[-npts:, :, :] = (
                    varpt[-2 * npts : -npts, :, :][::-1, :, :] * antiMult
                )
        if direction == "y1":
            if self.pyranda.PyMPI.y1proc:
                varpt[:, 0:npts, :] = (
                    varpt[:, npts : npts * 2, :][:, ::-1, :] * antiMult
                )
        if direction == "yn":
            if self.pyranda.PyMPI.ynproc:
                varpt[:, -npts:, :] = (
                    varpt[:, -2 * npts : -npts, :][:, ::-1, :] * antiMult
                )
        if direction == "z1":
            if self.pyranda.PyMPI.z1proc:
                varpt[:, :, 0:npts] = (
                    varpt[:, :, npts : npts * 2][:, :, ::-1] * antiMult
                )
        if direction == "zn":
            if self.pyranda.PyMPI.znproc:
                varpt[:, :, -npts:] = (
                    varpt[:, :, -2 * npts : -npts][:, :, ::-1] * antiMult
                )
