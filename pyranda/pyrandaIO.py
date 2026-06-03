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
import os
import numpy
import struct


class pyrandaIO:
    """
    Class to handle IO of pyranda objects and files
    """

    def __init__(self, rootname, pympi):

        self.rootname = rootname
        self.PyMPI = pympi
        self.ioformat = "BINARY"  # "BINARY" or "ASCII" only
        self.variables = []  # List of variables to be included

        if self.PyMPI.master == 1:
            try:
                os.mkdir(rootname)
            except Exception:
                pass
        self.PyMPI.comm.barrier()  # Wait for directory to be made

    def makeDump(self, data, dumpName):

        dumpFile = dumpName
        if self.PyMPI.master == 1:
            try:
                os.mkdir(os.path.join(self.rootname, dumpFile))
            except Exception:
                pass

        rank = self.PyMPI.comm.rank
        file_name = os.path.join(self.rootname, dumpFile, "p" + str(rank).zfill(8))

        fd = open(file_name, "wb")
        # fwrite(fd, data.size, data)
        data.tofile(fd)
        fd.close()

    def makeDumpVTK(self, mesh, variables, varList, dumpFile, cycle, time):

        if self.ioformat == "BINARY":
            self.makeDumpVTKbin(mesh, variables, varList, dumpFile, cycle, time)
        elif self.ioformat == "ASCII":
            self.makeDumpVTKascii(mesh, variables, varList, dumpFile, cycle, time)

    def makeDumpVTKascii(self, mesh, variables, varList, dumpFile, cycle, time):

        ghost = True

        if ghost:
            fx = self.PyMPI.ghost(mesh.coords[0].data[:, :, :])
            fy = self.PyMPI.ghost(mesh.coords[1].data[:, :, :])
            fz = self.PyMPI.ghost(mesh.coords[2].data[:, :, :])
        else:
            fx = mesh.coords[0].data[:, :, :]
            fy = mesh.coords[1].data[:, :, :]
            fz = mesh.coords[2].data[:, :, :]

        ax = fx.shape[0]
        ay = fx.shape[1]
        az = fx.shape[2]

        fid = open(dumpFile + ".vtk", "w")

        fid.write("# vtk DataFile Version 3.0 \n")
        fid.write("vtk output \n")
        fid.write("ASCII \n")
        fid.write("DATASET STRUCTURED_GRID \n")

        # For ASCII vtk, this gives proper time stamps
        fid.write("FIELD FieldData 2 \n")
        fid.write("CYCLE 1 1 int \n")
        fid.write("%d \n" % cycle)
        fid.write("TIME 1 1 double \n")
        fid.write("%s \n" % time)

        fid.write("DIMENSIONS  %s %s %s  \n" % (ax, ay, az))
        fid.write("POINTS %s float  \n" % (ax * ay * az))

        for k in range(az):
            for j in range(ay):
                for i in range(ax):
                    fid.write("%s " % fx[i, j, k])
                    fid.write("%s " % fy[i, j, k])
                    fid.write("%s " % fz[i, j, k])

        fid.write("\nPOINT_DATA %s  " % (ax * ay * az))

        for var in varList:
            fid.write("\nSCALARS %s float \n" % var)
            fid.write("LOOKUP_TABLE default \n")
            if ghost:
                gdata = self.PyMPI.ghost(variables[var].data)
            else:
                gdata = variables[var].data
            fid = open(dumpFile + ".vtk", "a")
            for k in range(az):
                for j in range(ay):
                    for i in range(ax):
                        fid.write("%s " % gdata[i, j, k])

        fid.close()

    def makeDumpVTKbin(
        self, mesh, variables, varList, dumpFile, cycle, time, make2d=True
    ):

        ghost = True

        if ghost:
            fx = self.PyMPI.ghost(mesh.coords[0].data[:, :, :])
            fy = self.PyMPI.ghost(mesh.coords[1].data[:, :, :])
            fz = self.PyMPI.ghost(mesh.coords[2].data[:, :, :])
        else:
            fx = mesh.coords[0].data[:, :, :]
            fy = mesh.coords[1].data[:, :, :]
            fz = mesh.coords[2].data[:, :, :]

        ax = fx.shape[0]
        ay = fx.shape[1]
        az = fx.shape[2]

        # Make2d option - for simulations with only 2 dimensional data, force them to be in first/second dimension
        #  Only valid when az > 1 and either ax==1 or ay==1
        reshape = None
        if make2d and (az > 1) and ((ax == 1) or (ay == 1)):
            if ax == 1:
                ax = ay
                ay = az
                az = 1
                fx[:, :, :] = fy
                fy[:, :, :] = fz
                fz[:, :, :] = 0.0
            elif ay == 1:
                ay = az
                az = 1
                fy[:, :, :] = fz
                fz[:, :, :] = 0.0

            # 2D dataset always looks like this
            reshape = (ax, ay, 1)

            # Reshape mesh
            fx = fx.reshape(reshape, order="F")
            fy = fy.reshape(reshape, order="F")
            fz = fz.reshape(reshape, order="F")

        # For nz == 1, check for z value an make it zero for viz
        if az == 1 and (fz.max() == fz.min()):
            fz[:, :, :] = 0.0

        fid = open(dumpFile + ".vtk", "w")

        fid.write("# vtk DataFile Version 3.0 \n")
        fid.write("vtk output \n")
        fid.write("BINARY \n")
        fid.write("DATASET STRUCTURED_GRID \n")
        fid.write("DIMENSIONS  %s %s %s  \n" % (ax, ay, az))

        # Add time/cycle stamps to vtk files
        fid.write("FIELD FieldData 2 \n")
        # Time meta
        fid.write("TIME 1 1 float \n")
        fid.close()
        # Time in binary
        fid = open(dumpFile + ".vtk", "ab")
        fid.write(struct.pack(">f", time))
        fid.close()
        # Cycle meta
        fid = open(dumpFile + ".vtk", "a")
        fid.write("CYCLE 1 1 int \n")
        fid.close()
        # Cycle in binary
        fid = open(dumpFile + ".vtk", "ab")
        fid.write(struct.pack(">i", cycle))
        fid.close()

        fid = open(dumpFile + ".vtk", "a")
        fid.write("POINTS %s float  \n" % (ax * ay * az))

        fid.close()
        fid = open(dumpFile + ".vtk", "ab")

        for k in range(az):
            for j in range(ay):
                for i in range(ax):
                    fid.write(struct.pack(">f", fx[i, j, k]))
                    fid.write(struct.pack(">f", fy[i, j, k]))
                    fid.write(struct.pack(">f", fz[i, j, k]))

        fid.close()

        fid = open(dumpFile + ".vtk", "a")
        fid.write("\nPOINT_DATA %s  " % (ax * ay * az))
        fid.close()

        for var in varList:
            fid = open(dumpFile + ".vtk", "a")
            fid.write("\nSCALARS %s float \n" % var)
            fid.write("LOOKUP_TABLE default \n")
            fid.close()
            if ghost:
                gdata = self.PyMPI.ghost(variables[var].data)
            else:
                gdata = variables[var].data
            if reshape:
                gdata = gdata.reshape(reshape, order="F")
            fid = open(dumpFile + ".vtk", "ab")
            for k in range(az):
                for j in range(ay):
                    for i in range(ax):
                        fid.write(struct.pack(">f", gdata[i, j, k]))
            fid.close()

    def makeDumpTec(self, mesh, variables, varList, dumpFile):

        fx = mesh.coords[0].data[:, :, :].flatten()
        fy = mesh.coords[1].data[:, :, :].flatten()
        fz = mesh.coords[2].data[:, :, :].flatten()

        fid = open(dumpFile + ".tec", "w")

        ax = self.PyMPI.chunk_3d_size[0]
        ay = self.PyMPI.chunk_3d_size[1]
        az = self.PyMPI.chunk_3d_size[2]

        fid.write('VARIABLES = "X", "Y", "Z" ')
        for var in varList:
            fid.write(', "%s" ' % var)
        fid.write("\n")

        fid.write("ZONE  I=%s, J=%s, K=%s, DATAPACKING=BLOCK  \n" % (ax, ay, az))

        # Mesh
        numpy.savetxt(fid, fx, fmt="%f")
        numpy.savetxt(fid, fy, fmt="%f")
        numpy.savetxt(fid, fz, fmt="%f")

        # Variables
        for var in varList:
            numpy.savetxt(fid, variables[var].data.flatten(), fmt="%f")

        fid.close()

    def makeMIR(self):

        # Write meta file for viz
        form = """
VERSION 2.0
zonal: yes                # Zonal or Nodal?
curvilinear: yes          # Read grid files?
gridfiles: grid/p%08d     # grid files
datafiles: vis%07d/p%08d  # data files
fileorder: XYZ # processor order
domainsize:       #AX#    #AY#    #AZ# # overall nodal dimensions
blocksize:        #AX#    #AY#    #AZ# # nodal dimensions per miranda processor
interiorsize:     #AX#    #AY#    #AZ# # nodal dimensions per interior processor
bndrysize:        #AX#    #AY#    #AZ# # nodal dimensions per boundary processor
origin:    #X1#  #Y1#  #Z1#            # xmin, ymin, zmin for rectilinear
spacing:   #DX#  #DY#  #DZ#            # dx, dy, dz for rectilinear
variables:   #NVARS#  # number of variables
  #VARS#
timesteps:    #NVIZ#  # number of times to plot
  #CYCTIME#
"""
        # form = form.replace('#AX#',str(self.PyMPI.chunk_3d_size[0]))
        # form = form.replace('#AY#',str(self.PyMPI.chunk_3d_size[1]))
        # form = form.replace('#AZ#',str(self.PyMPI.chunk_3d_size[2]))

        # form = form.replace('#X1#',str(self.mesh.options['x1'][0]))
        # form = form.replace('#Y1#',str(self.mesh.options['x1'][1]))
        # form = form.replace('#Z1#',str(self.mesh.options['x1'][2]))

        # form = form.replace('#DX#',str(self.dx))
        # form = form.replace('#DY#',str(self.dy))
        # form = form.replace('#DZ#',str(self.dz))

        # form = form.replace('#NVARS#',str(wlen))
        # svars = ''
        # for ivar in wVars:
        #    svars += "  %s 1 %f %f \n" % (ivar,
        #                                  self.PyMPI.min3D( self.variables[ivar].data),
        #                                  self.PyMPI.max3D( self.variables[ivar].data) )
        # form = form.replace("#VARS#",svars)

        # Time history
        # form = form.replace("#NVIZ#",str(len( self.vizDumpHistory) ))
        # sviz = ''
        # for iv in self.vizDumpHistory:
        #    sviz += "  %s  %f \n" % ( str(iv[0]).zfill(7), iv[1] )
        # form = form.replace("#CYCTIME#",sviz)

        # if self.PyMPI.master:
        #    mfid = open( os.path.join( self.PyIO.rootname, 'pyranda.mir'),'w+')
        #    mfid.write(form)
        #    mfid.close()
