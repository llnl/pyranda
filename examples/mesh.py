import numpy
import matplotlib.pyplot as plt
from pyranda import pyrandaSim


Npts = 64
L = numpy.pi * 2.0
Lp = L * (Npts - 1.0) / Npts

mesh_options = {}
mesh_options["type"] = "cartesian"
mesh_options["periodic"] = numpy.array([True, True, True])
mesh_options["dim"] = 3
mesh_options["x1"] = [0.0, 0.0, 0.0]
mesh_options["xn"] = [Lp, Lp, Lp]
mesh_options["nn"] = [Npts, Npts, 1]


ss = pyrandaSim("advection", mesh_options)

ss.addVar("phi", kind="conserved")
ss.addVar("Cx", kind="constant")
ss.addVar("Cy", kind="constant")

eq = "ddt(:phi:) =  - :Cx: * ddx( :phi: ) - :Cy: * ddy(:phi:)"
ss.addEqu(eq)

eq = ":phi:      =  fbar( :phi: )"
ss.addEqu(eq)

ss.allocate()


x = ss.mesh.coords[0].data
y = ss.mesh.coords[1].data
z = ss.mesh.coords[2].data
rad = numpy.sqrt((x - numpy.pi) ** 2 + (y - numpy.pi) ** 2)  # + (z-numpy.pi)**2  )

ss.variables["phi"].data = numpy.exp(-((rad) ** 2) / (0.8**2))
# ss.variables['phi'].data = numpy.sin( x )

theta = 45.0 * numpy.pi / 180.0

v = 1.0
ss.variables["Cx"].data = v * numpy.cos(theta)
ss.variables["Cy"].data = v * numpy.sin(theta)

d = L / numpy.cos(theta)

# plt.ion()
time = 0.0
viz = False


tt = d / v

dt_max = v / ss.mesh.nn[0] * 0.75


phi1 = ss.PyMPI.zbar(ss.variables["phi"].data)
xx = ss.PyMPI.zbar(x)
yy = ss.PyMPI.zbar(y)
if ss.PyMPI.master:
    plt.figure(1)
    plt.contour(xx, yy, phi1, 32)
    plt.figure(2)
    plt.plot(xx[:, Npts / 2], phi1[:, Npts / 2], "k--")


dt = dt_max
cnt = 1
while tt > time:
    time = ss.rk4(time, dt)
    dt = min(dt_max, (tt - time))

    # if time >= tt:
    #    break

    ss.iprint("%s -- %s" % (cnt, time))
    cnt += 1
    if viz:
        plt.figure(1)
        plt.clf()
        plt.plot(ss.variables["phi"].data[:, 16, 0])
        plt.pause(0.001)
        input("Press Enter")


phi = ss.PyMPI.zbar(ss.variables["phi"].data)
if ss.PyMPI.master:
    plt.figure(1)
    plt.contour(xx, yy, phi, 32)
    plt.figure(2)
    plt.plot(xx[:, Npts / 2], phi[:, Npts / 2], "b-")

    plt.figure(3)
    plt.plot(xx[:, Npts / 2], numpy.abs(phi[:, Npts / 2] - phi1[:, Npts / 2]), "b-")

    plt.show()
