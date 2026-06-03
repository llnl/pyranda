import sys
import time
import numpy
import matplotlib.pyplot as plt
from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep, pyrandaPoisson


# Try to get args for testing
try:
    Npts = int(sys.argv[1])
except Exception:
    Npts = 40

try:
    test = bool(int(sys.argv[2]))
except Exception:
    test = False

try:
    testName = sys.argv[3]
except Exception:
    testName = None

try:
    solver = sys.argv[4]
except Exception:
    solver = "scipy"


## Define a mesh
L = 1.0
Re = 100.0


mesh_options = {}
mesh_options["coordsys"] = 0
mesh_options["periodic"] = numpy.array([False, False, True])
mesh_options["x1"] = [0.0, 0.0, 0.0]
mesh_options["xn"] = [L, L, L]
mesh_options["nn"] = [Npts, Npts, 1]


# Initialize a simulation object on a mesh
pysim = pyrandaSim("lidDrivenCavity", mesh_options)
pysim.addPackage(pyrandaBC(pysim))
pysim.addPackage(pyrandaTimestep(pysim))

# pysim.linearSolver = solver
pyPoisson = pyrandaPoisson(pysim, solver)
pyPoisson.BC = {}
pyPoisson.BC["x1"] = ["val", 0]
pyPoisson.BC["xn"] = ["val", 0]
pyPoisson.BC["y1"] = ["val", 0]
pyPoisson.BC["yn"] = ["val", 0]
pysim.addPackage(pyPoisson)

eom = """
# Primary Equations of motion here
:us: = :u: +  :dt: * ( - :u: * ddx( :u: ) - :v: * ddy( :u: ) + lap(:u:)*:nu:  )
:vs: = :v: +  :dt: * ( - :u: * ddx( :v: ) - :v: * ddy( :v: ) + lap(:v:)*:nu:  )
:ps: = invDelta( 1.0/:dt: * ( ddx(:us:) + ddy(:vs:) ) )
:u: = :us: - :dt: * ddx(:ps:)
:v: = :vs: - :dt: * ddy(:ps:)
:u: = fbar(:u:)
:v: = fbar(:v:)
# Boundary conditions
bc.const(['u','v'],['x1','xn','y1'],0.0)
bc.const(['u'],['yn'],1.0)
bc.const(['v'],['yn'],0.0)
:umag: = sqrt( :u:*:u: + :v:*:v: )
:dt: = dt.courant(:u:,:v:,0.0,0.0)*.9
:dt: = numpy.minimum(:dt:,0.1 * dt.diff(:nu:,1.0))
"""

pysim.EOM(eom)

ics = """
:u:  *= 0.0
:v:  *= 0.0
:nu: += 1.0/REY
:dt: = .0001
:bcvar: *= 0.0
""".replace("REY", str(Re))
pysim.setIC(ics)


dt = 0.001
time = 0.0
tt = 5.0
cyc = 1

if test:
    tt = 2.0

while tt > time:
    # Simple euler - No RK4 (steady state solution)
    pysim.updateVars()

    if cyc % 50 == 0 and not test:
        pysim.plot.figure(1)
        pysim.plot.clf()
        pysim.plot.contourf("umag", 32)

        X = pysim.mesh.coords[0].data[::2, ::2, 0]
        Y = pysim.mesh.coords[1].data[::2, ::2, 0]
        u = pysim.var("u").data[::2, ::2, 0]
        v = pysim.var("v").data[::2, ::2, 0]
        plt.quiver(X, Y, u, v)
        plt.pause(0.01)

    dt = pysim.var("dt").data
    time += dt
    cyc += 1

    print("Cycle: %s ---  Time: %f ---  dt:  %f " % (cyc, time, dt))

# Plot final solution
ioff = int(Npts / 2)
uF = pysim.var("u")[ioff, :, 0]
yF = pysim.mesh.coords[1][ioff, :, 0]

if not test:
    pysim.plot.figure(2)
    pysim.plot.plot("u", slice2d="i=%s;k=0" % ioff)
    plt.plot(yF, uF, "ko")

# Save file for regression testing
if test:
    fname = testName + ".dat"
    numpy.savetxt(fname, (yF, uF))
    print(fname)
