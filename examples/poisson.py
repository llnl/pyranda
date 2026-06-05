import sys
import numpy
import matplotlib.pyplot as plt

from pyranda import pyrandaSim

# Try to get args
try:
    Npts = int(sys.argv[1])
except Exception:
    Npts = 100

try:
    test = bool(int(sys.argv[2]))
except Exception:
    test = False

## Define a mesh
L = numpy.pi * 2.0


# Define the domain/mesh
imesh = """
Lp = %s
Npts = %d
xdom = (0.0, Lp,  Npts, periodic=False)
ydom = (0.0, Lp,  Npts, periodic=False)
""" % (L, Npts)

# Initialize a simulation object on a mesh
ss = pyrandaSim("Poisson", imesh)

# Define the equations of motion
eom = """
:wx: = meshx * ( 2.0 * pi - meshx )
:wy: = meshy * ( 2.0 * pi - meshy )
:phi: =  cos( meshx ) * :wy: + cos( meshy ) * :wx:
Delta(:ff:) = :phi:
:gg: = lap(:ff:)
"""
ss.EOM(eom)

# Perform the solve at update
ss.updateVars()


if not test:
    ss.plot.figure(1)
    ss.plot.clf()
    ss.plot.contourf("phi", 32)
    ss.plot.figure(2)
    ss.plot.clf()
    ss.plot.contourf("ff", 32)
    plt.pause(0.001)

    ss.plot.figure(3)
    ss.plot.clf()
    ss.plot.contourf("gg", 32)
    plt.pause(0.001)
