from __future__ import print_function
import sys
import numpy
import matplotlib.pyplot as plt
from pyranda import pyrandaSim


# Try to get args
try:
    Npts = int(sys.argv[1])
except Exception:
    Npts = 10

try:
    test = bool(int(sys.argv[2]))
except Exception:
    test = False


# Define the domain/mesh
domain = "xdom = (0.0 , 1.0 , 1 )"

# Initialize a simulation object on a mesh
pysim = pyrandaSim("integralTests", domain)

# Define the equations of motion
eom = """
ddt(:a:)  =  :a:
ddt(:b:)  =  cos( 6.0*simtime )
"""
pysim.EOM(eom)

# Initialize variables
a0 = 1.0
b0 = 0.0
ic = """
:a: = %s
:b: = %s
""" % (a0, b0)
pysim.setIC(ic)

# Integrate in time
dt = 1.0 / float(Npts)
time = 0.0
tt = [time]
aa = [a0]
bb = [b0]
while time < 1.0:
    time = pysim.rk4(time, dt)
    tt.append(time)
    aa.append(pysim.variables["a"].data[0, 0, 0])
    bb.append(pysim.variables["b"].data[0, 0, 0])

# Plot the initial/final solution
tt = numpy.array(tt)
aa = numpy.array(aa)
a_sol = 1.0 - a0 + numpy.exp(tt)

bb = numpy.array(bb)
b_sol = b0 + 1.0 / 6.0 * numpy.sin(6.0 * tt)


if not test:
    plt.figure()
    plt.plot(tt, aa, "b-")
    plt.plot(tt, a_sol, "k--")

    plt.figure()
    plt.plot(tt, bb, "b-")
    plt.plot(tt, b_sol, "k--")


if test:
    a_error = numpy.sum(numpy.abs(aa - a_sol))
    b_error = numpy.sum(numpy.abs(bb - b_sol))
    print(a_error)


plt.show()
