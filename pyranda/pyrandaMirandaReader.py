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
import sys, os, glob
from os.path import join
from .pyranda import pyrandaSim

# Gaurd against hdf5 dependency
try:
    import h5py
except Exception:
    print("h5py is required for mirandaReader to work")
    pass

# Static paths for sidre/blueprint import
path_proc = "/sidre/groups/marbl/groups/miranda/groups/blueprint/views/"
path_coord = "/sidre/groups/marbl/groups/miranda/groups/blueprint/groups/coordsets/groups/coords/groups/dims/views/"
path_grid = "/sidre/groups/marbl/groups/miranda/groups/blueprint/groups/coordsets/groups/coords/groups/spacing/views/"
path_field = "/sidre/external/marbl/miranda/blueprint/fields/"
path_mats = "/sidre/external/marbl/miranda/blueprint/matsets/mats/volume_fractions/"

# Hard code in miranda native variable names/list
miranda_viz_names = {}
miranda_viz_names["u"] = 0
miranda_viz_names["v"] = 1
miranda_viz_names["w"] = 2
miranda_viz_names["density"] = 3
miranda_viz_names["energy"] = 4
miranda_viz_names["pressure"] = 5
miranda_viz_names["e_pressure"] = 6
miranda_viz_names["sound_speed"] = 7
miranda_viz_names["shear_viscosity"] = 8
miranda_viz_names["bulk_viscosity"] = 9
miranda_viz_names["thermal_conductivity"] = 10

# if ns > 1
miranda_viz_names["diffusivity"] = 11

# if there are no other physics
miranda_viz_names["mat1"] = 12
miranda_viz_names["mat2"] = 13


class mirandaReader:
    def __init__(self, rootdir, var_list, periodic=[False, False, False]):

        self.rootdir = rootdir  # Directory containing restart dirs
        self.var_list = var_list  # List of variables to be read in
        self.native_viz_files = []  # List of available viz files to be read in
        self.var_index_dict = {}

        # Rootfile is directory where the marbl_* files are located
        cycle = 0
        self.data_dir = join(rootdir, "marbl_" + str(cycle).zfill(7))

        # Check if file exists
        if not os.path.isdir(self.data_dir):
            print("Error: Cant read restart file %s" % self.data_dir)
            return None

        # Read in domain meta-data from master proc
        proc_num = 0
        filename = join(
            self.data_dir,
            "marbl_" + str(cycle).zfill(7) + "_" + str(proc_num).zfill(7) + ".hdf5",
        )
        header = "datagroup"  # + str(proc_num).zfill(7)

        # Load hdf5 file
        f = h5py.File(filename, "r")

        # Processor arrangement
        px = int(f[header + path_proc + "px/value"][0])
        py = int(f[header + path_proc + "py/value"][0])
        pz = int(f[header + path_proc + "pz/value"][0])

        # Grid spacing
        try:
            dx = f[header + path_grid + "dx/value"][0]
            self.ax = ax = int(f[header + path_coord + "i/value"][0]) - 1
        except Exception:
            dx = 1.0
            ax = 1

        try:
            dy = f[header + path_grid + "dy/value"][0]
            self.ay = ay = int(f[header + path_coord + "j/value"][0]) - 1
        except Exception:
            dy = 1.0
            ay = 1

        try:
            dz = f[header + path_grid + "dz/value"][0]
            self.az = az = int(f[header + path_coord + "k/value"][0]) - 1
        except Exception:
            dz = 1.0
            az = 1

        # Correction for x-z 2D domains
        if (ax > 1) and (ay > 1) and (pz > 1) and (py == 1):
            py = pz
            pz = 1

        self.procs = [px, py, pz]

        # Need to get a mesh object here
        x1, xn = 0, dx * ax * px
        y1, yn = 0, dy * ay * py
        z1, zn = 0, dz * az * pz
        nx, ny, nz = ax * px, ay * py, az * pz
        mesh_options = {}
        mesh_options["x1"] = [x1, y1, z1]
        mesh_options["xn"] = [xn, yn, zn]
        mesh_options["nn"] = [nx, ny, nz]
        mesh_options["periodic"] = periodic

        # Make pyranda instance
        self.pysim = pyrandaSim(rootdir + "_pyranda", mesh_options)

        # Add variables to pyranda-state
        for v in self.var_list:
            self.pysim.addVar(v)
        self.pysim.allocate()

        # Make a procmap for the source data/MPI
        self.procMap = {}
        ipp = 0
        for ipz in range(pz):
            for ipy in range(py):
                for ipx in range(px):
                    gx1 = ax * ipx + 1 - 1
                    gy1 = ay * ipy + 1 - 1
                    gz1 = az * ipz + 1 - 1

                    gxn = ax * (ipx + 1) - 1
                    gyn = ay * (ipy + 1) - 1
                    gzn = az * (ipz + 1) - 1

                    self.procMap["%s-g1" % ipp] = numpy.array([gx1, gy1, gz1])
                    self.procMap["%s-gn" % ipp] = numpy.array([gxn, gyn, gzn])

                    ipp += 1

    def update(self, cycle, fileType=0):
        """
        Import Miranda data from specified cycle
        """

        # Update data directory to current cycle
        if fileType == 0:  # marbl retart files
            self.data_dir = join(self.rootdir, "marbl_" + str(cycle).zfill(7))
        elif fileType == 1:  # native miranda vis files
            self.data_dir = join(self.rootdir, "vis" + str(cycle).zfill(7))
        elif fileType == 2:  # native miranda vis files
            self.data_dir = join(self.rootdir, "visit_dump." + str(cycle).zfill(5))

        # Check if file exists
        if not os.path.isdir(self.data_dir):
            self.pysim.iprint("Error: Cant read restart file %s" % self.data_dir)
            return None

        readChunkMiranda(
            self.pysim,
            self.procs,
            self.procMap,
            self.data_dir,
            cycle,
            self.var_list,
            self.var_index_dict,
            fileType=fileType,
        )

        return self.pysim

    def getNativeVizFiles(self):

        # Use glob to find list of native viz files and sort
        self.native_viz_files = glob.glob(join(self.rootdir, "vis0*"))
        self.native_viz_files.sort()

    def setNativeVizDict(self, rosetta):
        varDict = {}
        cnt = 0
        for ii in self.var_list:
            varDict[ii] = miranda_viz_names[rosetta[ii]]
        self.var_index_dict = varDict

    def getMirandaSamraiVariables(self):

        # Read in the first vis file to get variable list
        data_dir = join(self.rootdir, "visit_dump." + str(0).zfill(5))
        filename = join(data_dir, "processor_cluster." + str(0).zfill(5) + ".samrai")

        # Load hdf5 file
        data = h5py.File(filename, "r")

        # Use the 0-proc for variables
        sproc = str(0).zfill(5)
        header = "processor." + sproc + "/level.00000/patch." + sproc + "/"
        visit_vars = [dd for dd in data[header].keys()]
        return visit_vars


def readChunkMiranda(
    pysim, procs, procMap, dump, cycle, var_list, var_index_dict, fileType=0
):
    """
    Same as readData but only reads in global range of data given by
    irange.
    """

    Rx = [0] * 2
    Ry = [0] * 2
    Rz = [0] * 2

    # This procs extents
    Rx[0] = pysim.PyMPI.chunk_3d_lo[0]  # irange[0]
    Rx[1] = pysim.PyMPI.chunk_3d_hi[0] + 1  # irange[1]
    Ry[0] = pysim.PyMPI.chunk_3d_lo[1]  # irange[2]
    Ry[1] = pysim.PyMPI.chunk_3d_hi[1] + 1  # irange[3]
    Rz[0] = pysim.PyMPI.chunk_3d_lo[2]  # irange[4]
    Rz[1] = pysim.PyMPI.chunk_3d_hi[2] + 1  # irange[5]

    # Restart deomain info
    nprocs = procs[0] * procs[1] * procs[2]
    ax = pysim.nx / procs[0]
    ay = pysim.ny / procs[1]
    az = pysim.nz / procs[2]
    nshape = (ax, ay, az, len(pysim.variables))

    for iproc in range(nprocs):
        g1 = procMap["%s-g1" % iproc]
        gn = procMap["%s-gn" % iproc] + 1

        # Shift left point if node data
        iff = 0
        jff = 0
        kff = 0
        c1 = Rx[1] in range(g1[0], gn[0])
        c2 = Rx[0] in range(g1[0], gn[0])
        c3 = (g1[0] and gn[0]) in range(Rx[0], Rx[1] + 1)
        CX = c1 or c2 or c3

        c1 = Ry[1] in range(g1[1], gn[1])
        c2 = Ry[0] in range(g1[1], gn[1])
        c3 = (g1[1] and gn[1]) in range(Ry[0], Ry[1] + 1)
        CY = c1 or c2 or c3

        c1 = Rz[1] in range(g1[2], gn[2])
        c2 = Rz[0] in range(g1[2], gn[2])
        c3 = (g1[2] and gn[2]) in range(Rz[0], Rz[1] + 1)
        CZ = c1 or c2 or c3

        if CX and CY and CZ:
            Li1 = numpy.max((0, Rx[0] - g1[0])) + iff
            Lif = numpy.min((Rx[1], gn[0])) - g1[0] + iff
            Ki1 = numpy.max((Rx[0], g1[0])) - Rx[0]
            Kif = Ki1 + (Lif - Li1)

            Lj1 = numpy.max((0, Ry[0] - g1[1])) + jff
            Ljf = numpy.min((Ry[1], gn[1])) - g1[1] + jff
            Kj1 = numpy.max((Ry[0], g1[1])) - Ry[0]
            Kjf = Kj1 + (Ljf - Lj1)

            Lk1 = numpy.max((0, Rz[0] - g1[2])) + kff
            Lkf = numpy.min((Rz[1], gn[2])) - g1[2] + kff
            Kk1 = numpy.max((Rz[0], g1[2])) - Rz[0]
            Kkf = Kk1 + (Lkf - Lk1)

            # Read proc file
            proc_num = iproc
            if fileType == 0:  # MARBL restart files
                filename = join(
                    dump,
                    "marbl_"
                    + str(cycle).zfill(7)
                    + "_"
                    + str(proc_num).zfill(7)
                    + ".hdf5",
                )
                # Load hdf5 file
                f = h5py.File(filename, "r")

                for var in var_list:
                    DATA = import_field(f, iproc, var)
                    pysim.variables[var].data[Ki1:Kif, Kj1:Kjf, Kk1:Kkf] = DATA[
                        Li1:Lif, Lj1:Ljf, Lk1:Lkf
                    ]

            elif fileType == 1:  # Miranda native graphics files
                filename = join(dump, "p" + str(proc_num).zfill(8))

                # Load in the binary file
                fd = open(filename, "rb")
                f = numpy.fromfile(file=fd, dtype=numpy.single)
                fd.close()

                ibx = 0
                iby = 0
                ibz = 0

                shape = (int(ax + ibx), int(ay + iby), int(az + ibz))
                stride = numpy.product(shape) + 2
                for var in var_list:
                    ivar = var_index_dict[var]
                    istart = int(ivar * stride + 1)
                    iend = int((ivar + 1) * stride - 1)
                    DATA = f[istart:iend].reshape(shape, order="F")
                    pysim.variables[var].data[Ki1:Kif, Kj1:Kjf, Kk1:Kkf] = DATA[
                        Li1:Lif, Lj1:Ljf, Lk1:Lkf
                    ]

            elif fileType == 2:  # Samrai, single level only data
                filename = join(
                    dump, "processor_cluster." + str(proc_num).zfill(5) + ".samrai"
                )

                # Load hdf5 file
                f = h5py.File(filename, "r")
                file_grid_shape = (int(ax), int(ay), int(az))

                for var in var_list:
                    DATA = import_field_samrai(f, iproc, var, file_grid_shape)
                    pysim.variables[var].data[Ki1:Kif, Kj1:Kjf, Kk1:Kkf] = DATA[
                        Li1:Lif, Lj1:Ljf, Lk1:Lkf
                    ]


def convert(str_array):  # Convert numpy array of binary encoded chars to string
    s = ""
    for ele in str_array:
        s += ele.decode("UTF-8")
    return s


def import_field(data, proc_num, var):

    header = "datagroup"  # + str(proc_num).zfill(7)
    try:
        nx = int(data[header + path_coord + "i/value"][0]) - 1
    except Exception:
        nx = 1
    try:
        ny = int(data[header + path_coord + "j/value"][0]) - 1
    except Exception:
        ny = 1
    try:
        nz = int(data[header + path_coord + "k/value"][0]) - 1
    except Exception:
        nz = 1

    n = (nx, ny, nz)

    if var[0] == "Y":
        tmp = data[header + path_mats + str(var[1:]).zfill(3)]  # Temp data holder
        arr = numpy.zeros(
            tmp.shape, dtype="float64"
        )  # Initialize array using dataset shape
        tmp.read_direct(arr)  # Directly read in hdf5 dataset to numpy array
        # print('Finished importing ' + var + str(spec_num).zfill(3) + '!')

    else:
        tmp = data[header + path_field + var + "/values"]  # Temp data holder
        arr = numpy.zeros(
            tmp.shape, dtype="float64"
        )  # Initialize array using dataset shape
        tmp.read_direct(arr)  # Directly read in hdf5 dataset to numpy array
        # print('Finished importing ' + var + '!')

    return numpy.reshape(arr, n, order="F")


def import_field_samrai(data, proc_num, var, n):

    sproc = str(proc_num).zfill(5)
    header = "processor." + sproc + "/level.00000/patch." + sproc + "/"

    tmp = data[header + var]  # Temp data holder
    arr = numpy.zeros(
        tmp.shape, dtype="float64"
    )  # Initialize array using dataset shape
    tmp.read_direct(arr)  # Directly read in hdf5 dataset to numpy array

    ax = n[0]
    ay = n[1]
    az = n[2]

    xslc = slice(0, None)
    yslc = slice(0, None)
    zslc = slice(0, None)

    if ax > 1:
        ax += 2
        xslc = slice(1, -1)
    if ay > 1:
        ay += 2
        yslc = slice(1, -1)
    if az > 1:
        az += 2
        zslc = slice(1, -1)

    n_buffer = (ax, ay, az)

    amrdata = numpy.reshape(arr, n_buffer, order="F")

    return amrdata[xslc, yslc, zslc]
