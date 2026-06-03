from __future__ import print_function
import sys
import time
import numpy
import matplotlib.pyplot as plt

from pyranda import pyrandaSim, pyrandaBC

# Try to get args
try:
    Npts = int(sys.argv[1])
except Exception:
    Npts = 16

try:
    test = bool(sys.argv[2])
except Exception:
    test = False

## Define a mesh
L = numpy.pi * 2.0
Lp = L * (Npts - 1.0) / Npts


mesh_options = {}
mesh_options["coordsys"] = 0
mesh_options["periodic"] = numpy.array([False, True, True])
mesh_options["dim"] = 1
mesh_options["x1"] = [0.0, 0.0, 0.0]
mesh_options["xn"] = [Lp, Lp, Lp]
mesh_options["nn"] = [Npts, 1, 1]

# Initialize a simulation object on a mesh
ss = pyrandaSim("heat_equation", mesh_options)

ss.addPackage(pyrandaBC(ss))


# Define the equations of motion
eom = """
# Heat equation
ddt(:phi:)  =  :c: * lap(:phi:)
# Boundary condition
bc.const(['phi'],['x1'],2.0)
bc.const(['phi'],['xn'],1.0)
"""
ss.EOM(eom)

# Initialize variables
ic = """
xnn = meshx[-1,0,0]
:phi: = 1.0 + 1.0*(xnn - meshx)/xnn
:c:   = 1.0
"""
ss.setIC(ic)

x = ss.mesh.coords[0].data
xx = ss.PyMPI.zbar(x)

# Time step size
dt_max = L / ss.mesh.nn[0] * 0.005
tt = dt_max * 500


# Main time loop for physics
dt = dt_max
cnt = 1
time = 0.0
viz = True
while tt > time:
    time = ss.rk4(time, dt)
    dt = min(dt_max, (tt - time))

    if not test:
        ss.iprint("%s -- %s" % (cnt, time))

    # Plot animation of advection
    xnn = xx[-1, 0]
    anl = 1.0 + 1.0 * (xnn - xx) / xnn
    v = ss.PyMPI.zbar(ss.variables["phi"].data)

    error = numpy.abs(anl - v)

    cnt += 1
    if viz:
        if (ss.PyMPI.master and (cnt % 50 == 0)) and (not test):
            plt.figure(1)
            plt.clf()

            plt.plot(xx[:, 0], v[:, 0], "k.-")
            plt.plot(xx[:, 0], anl[:, 0], "b.-")
            plt.plot(xx[:, 0], error[:, 0], "b.-")
            plt.pause(0.001)

print(numpy.sum(error[:, 0]))
