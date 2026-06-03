import os
import shutil
import sys
import tempfile
import time

import numpy

if "MPLCONFIGDIR" not in os.environ:
    mplconfigdir = os.path.join(tempfile.gettempdir(), "mpl-pyranda")
    os.makedirs(mplconfigdir, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = mplconfigdir

if not any(os.environ.get(name) for name in ("DISPLAY", "WAYLAND_DISPLAY")):
    os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt

headless = "agg" in plt.get_backend().lower()

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
Lp = L * (Npts - 1.0) / Npts


# Define the domain/mesh
imesh = """
Lp = %s
Npts = %d
xdom = (0.0, Lp,  Npts, periodic=True)
""" % (Lp, Npts)

# Initialize a simulation object on a mesh
ss = pyrandaSim("advection", imesh)

# Define the equations of motion
ss.EOM(" ddt(:phi:)  =  -:c: * ddx(:phi:) ")


# Initialize variables
ic = """
r   = sqrt( (meshx-pi)**2  )
:phi: = 1.0 + 0.1 * exp( -(r/(pi/4.0))**2 )
:phi2: = 1.0*:phi:
:c:   = 1.0
"""
ss.setIC(ic)

# ss.variables["u"].data += 1.0

x = ss.mesh.coords[0].data
xx = ss.PyMPI.zbar(x)

# Time step size
v = 1.0
dt_max = v / ss.mesh.nn[0] * L * 0.90
tt = L / v * 1.0


# Main time loop for physics
dt = dt_max
cnt = 1
time = 0.0
viz = not headless
while tt > time:
    # raw_input('Pause...')

    time = ss.rk4(time, dt)
    dt = min(dt_max, (tt - time))

    if not test:
        ss.iprint("%s -- %s" % (cnt, time))

    # Plot animation of advection
    cnt += 1
    if viz:
        v = ss.PyMPI.zbar(ss.variables["phi"].data)
        if (ss.PyMPI.master and (cnt % 5 == 0)) and (not test):
            plt.figure(1)
            plt.clf()
            plt.plot(xx[:, 0], v[:, 0], "k.-")
            plt.pause(0.001)


phi = ss.variables["phi"].data
phi2 = ss.variables["phi2"].data
error = numpy.sum((phi - phi2) ** 2)
ss.iprint(error)


if test == 0 and not headless:
    # Latex/PDF generation
    ss.setupLatex()
    ss.latex.tMap[":phi:"] = r"\phi"
    ss.latex.tMap[":phi2:"] = r"\phi_2"

    intro = ss.latex.addSection("Introduction")
    intro.body = """
    This brief document shows the ability to integrate the simulation with formal 
    documentation.  The advection equations is solved and plotted below.
    """

    equations = ss.latex.addSection("Equations of Motion")
    equations.body = """
    Pyranda solves the following equations using a 10th order compact finite difference
    stencil and a 5 stage RK4 temporal integration scheme.
    """
    equations.body += ss.latex.renderEqu()

    initial = ss.latex.addSection("Initial Conditions", 1)
    initial.body = """
    The problem is initialized given the following formulae:
    """
    initial.body += ss.latex.renderEqu(kind="IC")

    results = ss.latex.addSection("Results")
    results.addFigure("Profiles", size=0.45)
    results.body = "Here in the figure, we see the final profile of the advection"

    plt.figure()
    plt.plot(xx[:, 0], ((phi - phi2) ** 2)[:, 0, 0])
    results.addFigure("Error", size=0.45)

    texcache = os.path.join(tempfile.gettempdir(), "texlive-pyranda")
    os.makedirs(texcache, exist_ok=True)
    os.environ.setdefault("TEXMFVAR", texcache)
    os.environ.setdefault("TEXMFCONFIG", texcache)

    ss.latex.makeDoc()
    if any(shutil.which(cmd) for cmd in ("xpdf", "open")) and not headless:
        ss.latex.showPDF()
