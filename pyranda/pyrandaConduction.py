"""
Copyright 2017 Ben Mather

This file is part of Conduction <https://git.dias.ie/itherc/conduction/>

Conduction is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or any later version.

Conduction is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Conduction.  If not, see <http://www.gnu.org/licenses/>.
"""

try:
    range = xrange
except Exception:
    pass

import numpy as np
from petsc4py import PETSc
from mpi4py import MPI

comm = MPI.COMM_WORLD


class MeshVariable(object):
    """
    Mesh variables live on the global mesh
    Every time its data is called a local instance is returned
    """

    def __init__(self, name, dm):
        self._dm = dm
        name = str(name)

        # mesh variable vector
        self._gdata = dm.createGlobalVector()
        self._ldata = dm.createLocalVector()

        self._gdata.setName(name)
        self._ldata.setName(name)

        self.size = self._ldata.getSizes()[0]

    def __delete__(self):
        self._ldata.destroy()
        self._gdata.destroy()

    def __getitem__(self, pos):
        self._dm.globalToLocal(self._gdata, self._ldata)
        return self._ldata[pos]

    def __setitem__(self, pos, value):
        self._ldata[pos] = value
        self._dm.localToGlobal(self._ldata, self._gdata)

    @property
    def array(self):
        self._dm.globalToLocal(self._gdata, self._ldata)
        return self._ldata

    @property
    def data(self):
        pass

    @data.getter
    def data(self):
        self._dm.globalToLocal(self._gdata, self._ldata)
        return self._ldata

    @data.setter
    def data(self, val):
        if type(val) is float:
            self._ldata.set(val)
            self._gdata.set(val)
        else:
            self._ldata.setArray(val)
            self._dm.localToGlobal(self._ldata, self._gdata)

    @data.deleter
    def data(self):
        self._ldata.destroy()
        self._gdata.destroy()

    def getGlobal(self):
        return self._gdata

    def getLocal(self):
        return self._ldata


def sum_duplicates(I, J, V):
    """
    Sum all duplicate entries in the matrix
    """
    order = np.lexsort((J, I))
    I, J, V = I[order], J[order], V[order]
    unique_mask = (I[1:] != I[:-1]) | (J[1:] != J[:-1])
    unique_mask = np.append(True, unique_mask)
    (unique_inds,) = np.nonzero(unique_mask)
    return I[unique_mask], J[unique_mask], np.add.reduceat(V, unique_inds)


class ConductionND(object):
    """
    Implicit N-dimensional solver for the steady-state heat equation
    over a structured grid using PETSc data structures.

    Parameters
    ----------
     minCoord : tuple, minimum Cartesian coordinates at edge of domain
     maxCoord : tuple, maximum Cartesian coordinates at edge of domain
     res      : tuple, resolution in each dimension
     kwargs   : dict, keyword arguments to pass to KSP method and preconditioner
        see PETSc documentaion for KSPType and PCType options...
        http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html
        http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html
    """

    def __init__(self, minCoord, maxCoord, res, **kwargs):

        dim = len(res)
        extent = np.zeros(dim * 2)

        index = 0
        for i in range(0, dim):
            extent[index] = minCoord[i]
            extent[index + 1] = maxCoord[i]
            index += 2

        width = kwargs.pop("stencil_width", 1)

        dm = PETSc.DMDA().create(dim=dim, sizes=res, stencil_width=width, comm=comm)
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

        n = np.zeros(dim, dtype=PETSc.IntType)
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
        self.nodes = np.arange(0, nn, dtype=PETSc.IntType)

        # set matrix and vector types
        self.MatType = kwargs.pop("MatType", "aij")  # cuda, seqaij, mpiaij, etc.
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

        # right hand side vector
        self.rhs = MeshVariable("rhs", dm)

    def __delete__(self):

        del self.rhs, self.diffusivity, self.heat_sources, self.temperature
        self.mat.destroy()
        self.dm.destroy()
        self.lvec.destroy()
        self.gvec.destroy()
        self.lgmap.destroy()

    def _initialise_ksp(self, matrix=None, atol=1e-10, rtol=1e-50, **kwargs):
        """
        Initialise linear solver object
        """
        if matrix is None:
            matrix = self.mat

        solver = kwargs.pop("solver", "gmres")
        precon = kwargs.pop("pc", None)

        ksp = PETSc.KSP().create(comm)
        ksp.setType(solver)
        ksp.setOperators(matrix)
        ksp.setTolerances(atol, rtol)
        if precon is not None:
            pc = ksp.getPC()
            pc.setType(precon)
        ksp.setFromOptions()
        return ksp

    def _initialise_COO_vectors(self, pad=1):

        nn = self.nn
        n = self.n

        self.index = np.pad(self.nodes.reshape(n), pad, "constant", constant_values=-1)

        self.rows = np.empty((self.stencil_width, nn), dtype=PETSc.IntType)
        self.cols = np.empty((self.stencil_width, nn), dtype=PETSc.IntType)
        self.vals = np.empty((self.stencil_width, nn))

    def _initialise_mesh_variables(self):

        dim = self.dim
        bbox = self.dm.getBoundingBox()

        extent = np.zeros(dim * 2)

        index = 0
        for bs, be in bbox:
            extent[index] = bs
            extent[index + 1] = be
            index += 2

        self.extent = extent

        # local coordinates
        self.coords = self.dm.getCoordinatesLocal().array.reshape(-1, dim)

        grid_coords = [None] * dim
        for i in range(0, dim):
            grid_coords[i] = np.unique(self.coords[:, i])

        self.grid_coords = grid_coords

    def _initialise_boundary_dictionary(self):

        coords = self.coords
        grid_coords = self.grid_coords
        dim = self.dim

        minCoords = coords.min(axis=0)
        maxCoords = coords.max(axis=0)

        bbox = self.dm.getBoundingBox()
        sizes = self.dm.getSizes()

        # Setup boundary dictionary
        bc = dict()

        wall = [("minX", "maxX"), ("minY", "maxY"), ("minZ", "maxZ")]

        for i in range(0, dim):
            w0, w1 = wall[i]
            c0, c1 = bbox[i]
            m0, m1 = self.coords[:, i] == c0, self.coords[:, i] == c1
            d0 = d1 = (c1 - c0) / (sizes[i] - 1)

            bc[w0] = {"val": 0.0, "delta": d0, "flux": True, "mask": m0}
            bc[w1] = {"val": 0.0, "delta": d1, "flux": True, "mask": m1}

        self.bc = bc
        self.dirichlet_mask = np.zeros(self.nn, dtype=bool)

    def _initialise_matrix(self, nnz=None):
        """
        There should be no mallocs but we turn off the error just to be sure.
        If there is it will be from users adjusting the BCs.

        Could push zeros into the matrix to allocate all potential entries
        but that would lengthen the build stage.
        """
        if nnz is None:
            nnz = (self.stencil_width, self.dim * 2)

        mat = PETSc.Mat().create(comm=comm)
        mat.setType(self.MatType)
        mat.setSizes(self.sizes)
        mat.setLGMap(self.lgmap)
        mat.setPreallocationNNZ(nnz)
        mat.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, 0)
        mat.setFromOptions()

        return mat

    def _initialise_vector(self, sizes):

        vec = PETSc.Vec().create(comm=comm)
        vec.setType(self.VecType)
        vec.setSizes(self.sizes[0])
        vec.setLGMap(self.lgmap)
        vec.setFromOptions()

        return vec

    def _create_closure_object(self, closure, pad=1):

        nc = len(closure)
        n = self.n
        p2 = 2 * pad
        obj = [[0] * self.dim for i in range(nc)]

        for i in range(0, nc):
            # construct slicing object
            for j in range(0, self.dim):
                start, end = closure[i - j]
                obj[i][j] = slice(start, n[j] + end + p2)

        return obj

    def _get_closure_array(self, dim, width=1, pad=1):
        w, p = width, pad
        if w > p:
            raise ValueError("width exceeds padding")

        if dim == 1:
            closure = [(p - w, -p - w), (p + w, -(p - w)), (p, -p)]
        elif dim == 2:
            closure = [(p - w, -p - w), (p, -p), (p + w, -(p - w)), (p, -p), (p, -p)]
        elif dim == 3:
            closure = [
                (p - w, -p - w),
                (p, -p),
                (p, -p),
                (p + w, -(p - w)),
                (p, -p),
                (p, -p),
                (p, -p),
            ]
        else:
            raise ValueError("{} is an invalid number of dimensions".format(dim))
        return closure

    def refine(self, fn, axis):
        """
        Pass a function to apply to the x,y,z coordinates on the mesh.
        The domain will be redefined accordingly.

        Notes
        -----
         We do it this way to make sure the domain is balanced across
         processors. Adding new nodes would imbalance the matrix.
        """
        v = self.dm.getCoordinatesLocal()
        coords = v.array.reshape(-1, self.dim)

        coords[:, axis] = fn(coords[:, axis])

        if not np.isfinite(coords).all():
            raise ValueError("This function has created NaNs or Inf numbers")

        v.setArray(coords.ravel())

        self.dm.setCoordinatesLocal(v)

        self._initialise_mesh_variables()
        self._initialise_boundary_dictionary()
        self.mat = self._initialise_matrix()

    def create_meshVariable(self, name):
        return MeshVariable(name, self.dm)

    def update_properties(self, diffusivity, heat_sources):
        """
        Update diffusivity and heat sources
        """

        self.diffusivity[:] = diffusivity
        self.heat_sources[:] = heat_sources

    def boundary_condition(self, wall, val, flux=True):
        """
        Set the boundary conditions on each wall of the domain.
        By default each wall is a Neumann (flux) condition.
        If flux=True, positive val indicates a flux vector towards the centre
        of the domain.

        val can be a vector with the same number of elements as the wall
        """
        wall = str(wall)

        if wall in self.bc:
            self.bc[wall]["val"] = np.array(val, copy=True)
            self.bc[wall]["flux"] = bool(flux)
            d = self.bc[wall]

            mask = d["mask"]

            if flux:
                self.dirichlet_mask[mask] = False
                self.bc[wall]["val"] /= -d["delta"]
            else:
                self.dirichlet_mask[mask] = True

        else:
            raise ValueError("Wall should be one of {}".format(self.bc.keys()))

    def find_neighbours(self, width=1):
        """
        Find node neighbours for each point in the domain

        Returns a point cloud of neighbours shape(n,nneighbours)
        -1 indicates a null value
        """

        nodes = self.nodes
        dim = self.dim
        n = self.n

        # setup new stencil
        stencil_width = 2 * self.dim * width + 1
        neighbours = np.empty((stencil_width, self.nn), dtype=PETSc.IntType)
        index = np.pad(nodes.reshape(n), width, "constant", constant_values=-1)

        closure = []
        for w in range(width, 0, -1):
            closure_array = self._get_closure_array(dim, w, width)
            closure.extend(closure_array[:-1])
        closure.append(closure_array[-1])  # centre node at last

        # create closure object
        closure = self._create_closure_object(closure, width)

        for i in range(0, stencil_width):
            obj = closure[i]
            neighbours[i] = index[obj].ravel()

        return neighbours.T

    def construct_matrix(self, in_place=True, derivative=False):
        """
        Construct the coefficient matrix
        i.e. matrix A in Ax = b

        We vectorise the 7-point stencil for fast matrix insertion.
        An extra border of dummy values around the domain allows for automatic
        Neumann (flux) boundary creation.
        These are stomped on if there are any Dirichlet conditions.

        """

        if in_place:
            mat = self.mat
        else:
            mat = self._initialise_matrix()

        nodes = self.nodes
        nn = self.nn
        n = self.n
        dim = self.dim

        index = self.index

        rows = self.rows
        cols = self.cols
        vals = self.vals

        dirichlet_mask = self.dirichlet_mask

        u = self.diffusivity[:].reshape(n)
        k = np.pad(u, self.width, "constant", constant_values=0)

        for i in range(0, self.stencil_width):
            obj = tuple(self.closure[i])

            rows[i] = nodes
            cols[i] = index[obj].ravel()

            distance = np.linalg.norm(self.coords[cols[i]] - self.coords, axis=1)
            distance[distance == 0] = 1e-12  # protect against dividing by zero
            delta = 1.0 / (2.0 * distance**2)

            vals[i] = delta * (k[obj] + u).ravel()

        # Dirichlet boundary conditions (duplicates are summed)
        cols[:, dirichlet_mask] = nodes[dirichlet_mask]
        vals[:, dirichlet_mask] = 0.0

        # zero off-grid coordinates
        vals[cols < 0] = 0.0

        # centre point
        vals[-1] = 0.0
        if derivative:
            vals[-1][dirichlet_mask] = 0.0
        else:
            vals[-1][dirichlet_mask] = -1.0

        row = rows.ravel()
        col = cols.ravel()
        val = vals.ravel()

        # mask off-grid entries and sum duplicates
        mask = col >= 0
        row, col, val = sum_duplicates(row[mask], col[mask], val[mask])

        # indptr, col, val = coo_tocsr(row, col, val)
        nnz = np.bincount(row)
        indptr = np.insert(np.cumsum(nnz), 0, 0)

        mat.assemblyBegin()
        mat.setValuesLocalCSR(indptr.astype(PETSc.IntType), col, val)
        mat.assemblyEnd()

        # set diagonal vector
        diag = mat.getRowSum()
        diag.scale(-1.0)
        mat.setDiagonal(diag)

        return mat

    def construct_rhs(self, in_place=True):
        """
        Construct the right-hand-side vector
        i.e. vector b in Ax = b

        Boundary conditions are grabbed from the dictionary and
        summed to the rhs.
        Be careful of duplicate entries on the corners!!
        """
        if in_place:
            rhs = self.rhs
        else:
            rhs = MeshVariable("rhs", self.dm)

        vec = -1.0 * self.heat_sources[:]

        for wall in self.bc:
            val = self.bc[wall]["val"]
            flux = self.bc[wall]["flux"]
            mask = self.bc[wall]["mask"]
            if flux:
                vec[mask] += val
            else:
                vec[mask] = val

        rhs[:] = vec
        return rhs

    def solve(self, matrix=None, rhs=None):
        """
        Construct the matrix A and vector b in Ax = b
        and solve for x

        GMRES method is default
        """
        if matrix is None:
            matrix = self.construct_matrix()
        if rhs is None:
            rhs = self.construct_rhs()
        res = self.temperature

        ksp = self.ksp
        ksp.setOperators(matrix)
        ksp.solve(rhs._gdata, res._gdata)
        # We should hand this back to local vectors
        return res[:]

    def sync(self, vector):
        """
        Synchronise a vector field across all processors
        """
        self.lvec.setArray(vector)
        self.dm.localToGlobal(self.lvec, self.gvec)
        self.dm.globalToLocal(self.gvec, self.lvec)
        return self.lvec.array.copy()

    def gradient(self, vector, **kwargs):
        """
        Computes the gradient of a vector field

        Parameters
        ----------
         vector : array shape(n,) size of the mesh

        Returns
        -------
         grad   : array shape(dim,n) gradient in each direction
        """
        return np.gradient(vector.reshape(self.n), *self.grid_coords[::-1], **kwargs)

    def heatflux(self):
        """
        Compute the heat flux based on stored vector fields
        - temperature
        - diffusivity

        Returns
        -------
         Q : array shape(dim,n) heat flux in each direction
        """
        T = self.temperature[:]
        k = self.diffusivity[:] * -1
        divT = np.array(self.gradient(T))
        q = []
        for i in range(0, self.dim):
            div = k * divT[i].ravel()
            q.append(div)

        return np.array(q)

    def isosurface(
        self, vector, isoval, axis=0, interp="nearest", return_indices=False
    ):
        """
        Calculate an isosurface along a given axis
        (So far this is only working for axis=0 and in serial)

        Parameters
        ----------
         vector : array, the same size as the mesh (n,)
         isoval : float, isosurface value
         axis   : int, axis to generate the isosurface
         interp : str, method can be either
            'nearest' - nearest neighbour interpolation
            'linear'  - linear interpolation
         return_indices : bool (default=False)
            return the coordinate index

        Returns
        -------
         z_interp : isosurface the same size as the specified axis
         indices  : int, same size as axis (if return_indices is True)
        """
        Vcube = vector.reshape(self.n)
        Zcube = self.coords[:, ::-1][:, axis].reshape(self.n)
        sort_idx = ((Vcube - isoval) ** 2).argsort(axis=axis)
        i0 = sort_idx[0]
        # z0 = Zcube.take(i0)

        obj = []
        for d in range(0, self.dim):
            obj.append(slice(0, self.n[d]))
        obj.pop(axis)

        idx = list(np.mgrid[obj])
        idx.insert(axis, i0)
        z0 = Zcube[idx]

        z_interp = z0
        if interp == "linear":
            v0 = Vcube[idx]

            # identify next nearest node
            i1 = sort_idx[1]
            idx[axis] = i1
            z1 = Zcube[idx]
            v1 = Vcube[idx]

            vmin = np.minimum(v0, v1)
            vmax = np.maximum(v0, v1)
            ratio = np.vstack([np.ones_like(vmax) * isoval, vmin, vmax])
            ratio -= ratio.min(axis=0)
            ratio /= ratio.max(axis=0)
            z_interp = ratio[0] * z1 + (1.0 - ratio[0]) * z0

        if return_indices:
            Bcube = np.zeros_like(Vcube, dtype=bool)
            idx[axis] = i0
            Bcube[idx] = True
            indices = np.where(Bcube.ravel())[0]
            return z_interp, indices
        else:
            return z_interp

    def save_mesh_to_hdf5(self, filename):
        """
        Save important mesh information to an HDF5 file
        including bounding box and resolution.

        These are saved under the topology group

        Parameters
        ----------
         filename : file to put .h5 file
        """

        import h5py

        filename = str(filename)
        if not filename.endswith(".h5"):
            filename += ".h5"

        ViewHDF5 = PETSc.Viewer()
        ViewHDF5.createHDF5(filename, mode="w")
        ViewHDF5.view(obj=self.dm)
        ViewHDF5.destroy()

        if comm.rank == 0:
            # Every processor is writing the same thing
            f = h5py.File(filename, "r+")
            f.create_group("topology")
            topo = f["topology"]

            # create attributes
            extent = self.extent.reshape(self.dim, -1)
            minCoord = extent[:, 0]
            maxCoord = extent[:, 1]
            shape = self.dm.getSizes()

            topo.attrs.create("minCoord", minCoord[::-1])
            topo.attrs.create("maxCoord", maxCoord[::-1])
            topo.attrs.create("shape", np.array(shape)[::-1])

            f.close()

    def save_field_to_hdf5(self, filename, *args, **kwargs):
        """
        Saves data on the mesh to an HDF5 file
         e.g. height, rainfall, sea level, etc.

        Pass these as arguments or keyword arguments for
        their names to be saved to the hdf5 file
        """
        import os.path

        filename = str(filename)
        if not filename.endswith(".h5"):
            filename += ".h5"

        # write mesh if it doesn't exist
        # if not os.path.isfile(file):
        #     self.save_mesh_to_hdf5(file)

        kwdict = kwargs
        for i, arg in enumerate(args):
            key = "arr_{}".format(i)
            if key in kwdict.keys():
                raise ValueError(
                    "Cannot use un-named variables\
                                  and keyword: {}".format(key)
                )
            kwdict[key] = arg

        vec = self.gvec.duplicate()

        # change mode to append if file already exists
        # set mode to "a" after first write
        if os.path.isfile(filename):
            mode = "a"
        else:
            mode = "w"

        for key in kwdict:
            val = kwdict[key]
            try:
                vec.setArray(val)
            except Exception:
                self.lvec.setArray(val)
                self.dm.localToGlobal(self.lvec, vec)

            vec.setName(key)

            ViewHDF5 = PETSc.Viewer()
            ViewHDF5.createHDF5(filename, mode=mode)
            ViewHDF5.view(obj=vec)
            ViewHDF5.destroy()
            mode = "a"

        vec.destroy()

    def save_vector_to_hdf5(self, filename, *args, **kwargs):
        """
        Saves vector on the mesh to an HDF5 file
         e.g. heat flux field.

        Pass these as arguments or keyword arguments for
        their names to be saved to the hdf5 file

        Each argument with x,y,z direction tuple
         e.g. Q=(Qx, Qy, Qz)
        """
        import os.path

        filename = str(filename)
        if not filename.endswith(".h5"):
            filename += ".h5"

        kwdict = kwargs
        for i, arg in enumerate(args):
            key = "arr_{}".format(i)
            if key in kwdict.keys():
                raise ValueError(
                    "Cannot use un-named variables\
                                  and keyword: {}".format(key)
                )
            kwdict[key] = arg

        # change mode to append if file already exists
        # set mode to "a" after first write
        if os.path.isfile(filename):
            mode = "a"
        else:
            mode = "w"

        # This is a flattened dim x n global vector
        gvec = self.dm.getCoordinates().duplicate()

        for key in kwdict:
            val = np.array(kwdict[key]).T.ravel()

            # vx, vy, vz = kwdict[key]
            # val = np.column_stack([vx, vy, vz]).ravel()

            gvec.assemblyBegin()
            gvec.setValuesLocal(np.arange(val.size, dtype=PETSc.IntType), val)
            gvec.assemblyEnd()
            gvec.setName(key)

            ViewHDF5 = PETSc.Viewer()
            ViewHDF5.createHDF5(filename, mode=mode)
            ViewHDF5.view(obj=gvec)
            ViewHDF5.destroy()
            mode = "a"

        gvec.destroy()
