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
import scipy.sparse
from scipy.sparse.linalg import factorized
from .pyrandaPackage import pyrandaPackage


class pyrandaPoisson(pyrandaPackage):
    def __init__(self, pysim, solver_type="scipy"):
        """
        Solve an elliptical Poisson equation
        \Delta \phi = f(x)
        Second order FD on structured mesh assumed
        """
        PackageName = "Poisson"
        self.pysim = pysim
        pyrandaPackage.__init__(self, PackageName, pysim)

        self.BCtype = None
        self.BCdata = {}
        self.BC = {}
        self.BCvar = ""
        self.solver_type = solver_type

        nx = self.nx = pysim.mesh.options["nn"][0]
        ny = self.ny = pysim.mesh.options["nn"][1]
        nz = self.nz = pysim.mesh.options["nn"][2]  # noqa: F841 (used implicitly)

        self.condMesh = None

        # Check petsc solver availability
        if solver_type == "petsc":
            try:
                import petsc4py  # noqa: F401 (ignore attempted imports)
                from petsc4py import PETSc
            except Exception:
                pysim.iprint("Error importing petsc4py:")
                pysim.iprint(" Defaulting to scipy solver.")
                self.solver_type = "scipy"

        ## CUSTOM SCIPY SOLVER.... ONLY 2D serial, non-periodic
        if solver_type == "scipy":
            if self.nz > 1:
                pysim.iprint("Error: Only 2D problems supported")

            dx = (pysim.mesh.options["xn"][0] - pysim.mesh.options["x1"][0]) / nx
            n = self.n = nx * ny
            d = numpy.ones(n)

            d0 = d.copy() * -4.0
            d1_lower = d.copy()[0:-1]
            d1_upper = d1_lower.copy()

            dnx_lower = d.copy()[0:-nx]
            dnx_upper = dnx_lower.copy()

            # Try to set both sides as Nuemann
            d1_upper[nx - 1 :: nx] = 0.0
            d1_upper[::nx] = 2.0

            d1_lower[nx - 1 :: nx] = 0.0
            d1_lower[nx - 2 :: nx] = 2.0

            dnx_upper[0:nx] = 2.0
            dnx_lower[-nx:] = 2.0

            # Avoid singular point
            d0[nx + 1] = 1.0
            d1_upper[nx + 1] = 0.0
            dnx_upper[nx + 1] = 0.0
            d1_lower[nx] = 0.0
            dnx_lower[1] = 0.0

            d0 /= dx * dx
            d1_upper /= dx * dx
            d1_lower /= dx * dx
            dnx_upper /= dx * dx
            dnx_lower /= dx * dx
            A = scipy.sparse.diags(
                [d0, d1_upper, d1_lower, dnx_upper, dnx_lower],
                [0, 1, -1, nx, -nx],
                format="csc",
            )

            self.solver = factorized(A)
            self.A = A

        elif self.solver_type == "petsc":
            px = self.px = pysim.PyMPI.px
            py = self.py = pysim.PyMPI.py
            nx = self.nx
            ny = self.ny
            res = (nx, ny)
            minCoord = pysim.mesh.options["x1"]
            maxCoord = pysim.mesh.options["xn"]
            comm = pysim.PyMPI.comm

            from .pyrandaConduction import ConductionND, MeshVariable

            class pyrandaCond(ConductionND):
                def __init__(self, **kwargs):

                    dim = len(res)
                    extent = numpy.zeros(dim * 2)

                    index = 0
                    for i in range(0, dim):
                        extent[index] = minCoord[i]
                        extent[index + 1] = maxCoord[i]
                        index += 2

                    width = kwargs.pop("stencil_width", 1)

                    dm = PETSc.DMDA().create(
                        dim=dim,
                        sizes=res,
                        stencil_width=width,
                        proc_sizes=[px, py],
                        comm=comm,
                    )
                    dm.setUniformCoordinates(*extent)

                    self.dm = dm
                    self.lgmap = dm.getLGMap()
                    self.lvec = dm.createLocalVector()
                    self.gvec = dm.createGlobalVector()

                    # Setup matrix sizes
                    self.sizes = self.gvec.getSizes(), self.gvec.getSizes()
                    self.dim = dim
                    self.extent = extent

                    # include ghost nodes in local domain
                    # (minI, maxI), (minJ, maxJ), (minK, maxK) = dm.getGhostRanges()
                    ghost_ranges = dm.getGhostRanges()

                    n = numpy.zeros(dim, dtype=PETSc.IntType)
                    nn = 1
                    for i, (gs, ge) in enumerate(ghost_ranges):
                        n[i] = ge - gs
                        nn *= n[i]

                    self.n = n[::-1]
                    self.nn = nn
                    self.npoints = nn

                    # stencil size
                    self.width = width
                    self.stencil_width = 2 * dim * width + 1

                    # create closure array
                    closure = []
                    for w in range(width, 0, -1):
                        closure_array = self._get_closure_array(dim, w, width)
                        closure.extend(closure_array[:-1])
                        closure.append(closure_array[-1])  # centre node at last

                    # create closure object
                    self.closure = self._create_closure_object(closure, width)

                    # local numbering
                    self.nodes = numpy.arange(0, nn, dtype=PETSc.IntType)

                    # set matrix and vector types
                    self.MatType = kwargs.pop(
                        "MatType", "aij"
                    )  # cuda, seqaij, mpiaij, etc.
                    self.VecType = kwargs.pop("VecType", "standard")

                    self._initialise_mesh_variables()
                    self._initialise_boundary_dictionary()
                    self.mat = self._initialise_matrix()
                    self._initialise_COO_vectors(width)
                    self.ksp = self._initialise_ksp(**kwargs)

                    # thermal properties
                    self.diffusivity = MeshVariable("diffusivity", dm)
                    self.heat_sources = MeshVariable("heat_sources", dm)
                    self.temperature = MeshVariable("temperature", dm)

                    # BCvector
                    self.bcarray = MeshVariable("bcarray", dm)

                    # right hand side vector
                    self.rhs = MeshVariable("rhs", dm)

            # Make the conduction object
            # self.condMesh = pyrandaCond(solver="cgs",rtol=1.0e-6,atol=1e-8,
            # precon="lu")

            self.condMesh = pyrandaCond(
                solver="cgs", rtol=1.0e-6, atol=1e-8, precon="lu"
            )

            # We're not using k, just H, which is our RHS
            k = numpy.ones(self.condMesh.nn)
            H = numpy.ones(self.condMesh.nn)
            self.condMesh.update_properties(k, H)
            self.matrix = None

        else:
            pysim.iprint(
                "Error: No valid elliptical solver given: %s" % self.solver_type
            )
            exit()

    def get_sMap(self):
        sMap = {}
        sMap["invDelta("] = "self.packages['Poisson'].solve("
        self.sMap = sMap

    def solve(self, rhs):

        if self.solver_type == "petsc":
            # Some indices to work with
            mx, my = self.condMesh.dm.getSizes()
            (xs, xe), (ys, ye) = self.condMesh.dm.getRanges()

            # Get a working array to move rhs->heat_sources
            b = self.condMesh.dm.getVecArray(self.condMesh.heat_sources._gdata)

            # Set RHS to b (heat_sources)
            b[xs:xe, ys:ye] = -1.0 * rhs[0 : xe - xs, 0 : ye - ys, 0]  # FIX FOR 3D

            # Set RHS from heat_sources
            self.condMesh.rhs[:] = -1.0 * self.condMesh.heat_sources[:]

            # Set the BCs
            # Directly set BC to RHS
            rhs = self.condMesh.dm.getVecArray(self.condMesh.rhs._gdata)
            bcmap = {"x1": "minX", "xn": "maxX", "y1": "minY", "yn": "maxY"}
            for bcstr in bcmap:
                if bcstr in self.BC:
                    val = self.BC[bcstr][0]
                    var = self.BC[bcstr][1]
                    opt = bcmap[bcstr]
                    isFlux = "flux" in val

                    ## FIX for 3D
                    hasBC = False
                    if bcstr == "x1" and xs == 0:
                        # rhs[xs,:] = data[0,:]
                        c_ind = numpy.index_exp[xs, :]
                        m_ind = numpy.index_exp[0, :]
                        hasBC = True
                    if bcstr == "xn" and xe == self.nx:
                        # rhs[xe-1,:] = data[-1,:]
                        c_ind = numpy.index_exp[xe - 1, :]
                        m_ind = numpy.index_exp[-1, :]
                        hasBC = True
                    if bcstr == "y1" and ys == 0:
                        # rhs[:,ys] = data[:,0]
                        c_ind = numpy.index_exp[:, ys]
                        m_ind = numpy.index_exp[:, 0]
                        hasBC = True
                    if bcstr == "yn" and ye == self.ny:
                        # rhs[:,ye-1] = data[:,-1]
                        c_ind = numpy.index_exp[:, ye - 1]
                        m_ind = numpy.index_exp[:, -1]
                        hasBC = True

                    # If any BCs are caught, update them here.
                    if hasBC:
                        val = 0
                        if type(var) is type(""):
                            val = self.pysim.var(var).data[:, :, 0][m_ind]

                        if isFlux:
                            rhs[c_ind] += val
                        else:
                            rhs[c_ind] = val

                    # Must call this to set the dirichlet masks
                    self.condMesh.boundary_condition(opt, 0.0, flux=isFlux)

            if not self.matrix:
                self.matrix = self.condMesh.construct_matrix()
                self.condMesh.ksp.setOperators(self.matrix)

            # Do the solve
            res = self.condMesh.temperature
            self.condMesh.ksp.solve(self.condMesh.rhs._gdata, res._gdata)

            # Unpack the solution onto our mesh
            u = self.condMesh.dm.createNaturalVec()
            self.condMesh.dm.globalToNatural(res._gdata, u)

            # Reshape for pyranda --- FIX FOR 3D
            myshape = (1, int(self.ny / self.py), int(self.nx / self.px))
            mysol = res._gdata[...].reshape(myshape, order="C").T

        elif self.solver_type == "scipy":
            b = numpy.zeros(self.n)  # RHS
            for j in range(0, self.ny):
                for i in range(0, self.nx):
                    b[j + i * self.ny] = rhs[i, j, 0]

            sol = self.solver(b)
            mysol = rhs * 0.0

            for j in range(0, self.ny):
                for i in range(0, self.nx):
                    mysol[i, j, 0] = sol[j + i * self.ny]

        else:
            print("Error: No valid elliptical solver given: %s" % self.solver_type)
            exit()

        return mysol
