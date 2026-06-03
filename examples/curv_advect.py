import numpy
import matplotlib.pyplot as plt
from matplotlib import cm
from pyranda import pyrandaSim


## Define a mesh
L = numpy.pi * 2.0
Npts = 200
Lp = L * (Npts - 1.0) / Npts

dx = L / float(Npts)


def myMesh(i, j, k):
    x = (i - 1) * dx + 0.5 * numpy.sin(j / float(Npts) * 2.0 * numpy.pi)
    y = (j - 1) * dx - 0.5 * numpy.sin(i / float(Npts) * 2.0 * numpy.pi)
    z = 0.0
    return x, y, z


imesh = {}

imesh["x1"] = [0.0, 0.0, 0.0]
imesh["xn"] = [Lp, Lp, 1]
imesh["nn"] = [Npts, Npts, 1]
imesh["periodic"] = numpy.array([True, True, True])
imesh["coordsys"] = 3
imesh["dim"] = 3
imesh["function"] = myMesh


ss = pyrandaSim("curved_advection", imesh)


eom = """
ddt(:phi:) = -div( :phi: * :u: , :phi:*:v: )
"""

ss.EOM(eom)

ic = """
r   = sqrt( (meshx-pi)**2 + (meshy-pi)**2  )
:u: = 1.0
:v: = 0.0
:phi: = 1.0 + 0.1 * exp( -(r/(pi/4.0))**2 )
"""

# Set the initial conditions
ss.setIC(ic)


# Mesh for viz on master
x = ss.mesh.coords[0].data
y = ss.mesh.coords[1].data
xx = ss.PyMPI.zbar(x)
yy = ss.PyMPI.zbar(y)

# Time step size
v = 1.0
dt_max = v / ss.mesh.nn[0] * L * 0.590
tt = L / v * 1.0


# Main time loop for physics
dt = dt_max
cnt = 1
time = 0.0
viz = True
viz_freq = 25
while tt > time:
    # raw_input('Pause...')

    time = ss.rk4(time, dt)
    dt = min(dt_max, (tt - time))

    ss.iprint("%s -- %s" % (cnt, time))

    phi = ss.PyMPI.zbar(ss.variables["phi"].data)
    if (ss.PyMPI.master) and (cnt % viz_freq == 1):  # or True:
        plt.figure(1)
        plt.clf()
        plt.contourf(xx, yy, phi, 64, cmap=cm.jet)
        plt.pause(0.001)

    cnt += 1


phi = ss.PyMPI.zbar(ss.variables["phi"].data)
if ss.PyMPI.master:
    plt.figure(1)
    plt.clf()
    plt.contourf(xx, yy, phi, 64, cmap=cm.jet)
    plt.pause(0.001)
