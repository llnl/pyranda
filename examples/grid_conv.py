import numpy
import matplotlib.pyplot as plt
from pyranda import pyrandaSim


def gridAdv(npts, cfl=1.0):

    # Define the domain/mesh
    xf = 1.0 * (npts - 1) / npts
    domain = "xdom = (0.0 , xf , npts , periodic=True)"
    domain = domain.replace("xf", str(xf))
    domain = domain.replace("npts", str(npts))

    # Initialize a simulation object on a mesh
    pysim = pyrandaSim("advection", domain, silent=True)

    # Define the equations of motion
    pysim.EOM(" ddt(:phi:)  =  -:c: * ddx(:phi:) ")

    # Initialize variables
    ic = """
:phi: = 1.0 + 0.1 * exp( -(abs(meshx-.5)/.1 )**2 )
:phi0: = :phi:
:c:   = 1.0
"""
    pysim.setIC(ic)

    # Integrate in time
    dt = 1.0 / (npts - 1) / 1.0 * cfl
    time = 0.0
    tfinal = 1.0
    while time < tfinal:
        time = pysim.rk4(time, dt)
        dt = min(dt, (tfinal - time))

    phi = pysim.variables["phi"].data
    phi0 = pysim.variables["phi0"].data

    error = numpy.sum(numpy.abs(phi - phi0))
    return error


npts = [16, 32, 64, 128, 256]  # ,512]
error = []

cfl = 0.25

for npt in npts:
    error.append(gridAdv(npt, cfl))

plt.loglog(npts, error, "k-o")
# 10th order fiducial
plt.loglog(npts, 1.0e10 * (1.0 / numpy.array(npts)) ** 10, "k--")

# 10th order fiducial
fourth = 1.0e-2 * (1.0 / numpy.array(npts)) ** 4
fact = fourth[-1] / error[-1] * 10.0
plt.loglog(npts, fourth / fact, "k--")
plt.xlabel("Npts")
plt.ylabel("L2-Error")
plt.title("CFL = %s" % cfl)
plt.show()
