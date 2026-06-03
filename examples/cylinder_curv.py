from __future__ import print_function
import sys
import time
import numpy
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep, pyrandaIBM, pyrandaProbes


# Try to get args for testing
try:
    Npts = int(sys.argv[1])
except Exception:
    Npts = 64

try:
    test = bool(int(sys.argv[2]))
except Exception:
    test = False

try:
    testName = sys.argv[3]
except Exception:
    testName = None


## Define a mesh
# Npts = 32
L = numpy.pi * 2.0
dim = 2
gamma = 1.4

problem = "cylinder_curvilinear"

Lp = L * (Npts - 1.0) / Npts

from meshTest import zoomMesh_solve

dxf = 4 * Lp / float(Npts) * 0.3
xS = zoomMesh_solve(Npts, -2.0 * Lp, 2.0 * Lp, -2.0, 2.0, 1.0, dxf)


def zoomMesh(i, j, k):
    x = xS[i]
    y = xS[j]
    z = 0.0
    return x, y, z


mesh_options = {}
mesh_options["coordsys"] = 3
mesh_options["function"] = zoomMesh
mesh_options["periodic"] = numpy.array([False, False, True])
mesh_options["periodicGrid"] = False
mesh_options["dim"] = 3
mesh_options["x1"] = [-2 * Lp, -2 * Lp, 0.0]
mesh_options["xn"] = [2 * Lp, 2 * Lp, Lp]
mesh_options["nn"] = [Npts, 1, 1]
if dim == 2:
    mesh_options["nn"] = [Npts, Npts, 1]


# Initialize a simulation object on a mesh
ss = pyrandaSim(problem, mesh_options)
ss.addPackage(pyrandaBC(ss))
ss.addPackage(pyrandaIBM(ss))
ss.addPackage(pyrandaTimestep(ss))
del mesh_options


rho0 = 1.0
p0 = 1.0
gamma = 1.4
mach = 2.0
s0 = numpy.sqrt(p0 / rho0 * gamma)
u0 = s0 * mach
e0 = p0 / (gamma - 1.0) + rho0 * 0.5 * u0 * u0


# Define the equations of motion
eom = """
# Primary Equations of motion here
ddt(:rho:)  =  -div(:rho:*:u:,  :rho:*:v:)
ddt(:rhou:) =  -div(:rhou:*:u: + :p: - :tau:, :rhou:*:v:)
ddt(:rhov:) =  -div(:rhov:*:u:, :rhov:*:v: + :p: - :tau:)
ddt(:Et:)   =  -div( (:Et: + :p: - :tau:)*:u: - :tx:*:kappa:, (:Et: + :p: - :tau:)*:v: - :ty:*:kappa: )
# Level set equation
#ddt(:phi:)  =  - :gx: * :u1: - :gy: * :v1: 
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v:) ) * ( :gamma: - 1.0 )
:T:         = :p: / (:rho: * :R: )
# Artificial bulk viscosity (old school way)
:div:       =  div(:u:,:v:)
:beta:      =  gbar( ring(:div:) * :rho: ) * 7.0e-3
:tau:       =  :beta: * :div: 
[:tx:,:ty:,:tz:] = grad(:T:)
:kappa:     = gbar( ring(:T:)* :rho:*:cv:/(:T: * :dt: ) ) * 1.0e-3
# Apply constant BCs
[:u:,:v:,:w:] = ibmV( [:u:,:v:,:w:], :phi:, [:gx:,:gy:,:gz:], [:u1:,:u2:,0.0] )
:rho: = ibmS( :rho: , :phi:, [:gx:,:gy:,:gz:] )
:p:   = ibmS( :p:   , :phi:, [:gx:,:gy:,:gz:] )
bc.extrap(['rho','p','u'],['xn'])
#bc.extrap(['u','rho','p'],['y1','yn'])
bc.const(['u'],['x1','y1','yn'],u0)
bc.const(['v'],['x1','xn','y1','yn'],0.0)
bc.const(['rho'],['x1','y1','yn'],rho0)
bc.const(['p'],['x1','y1','yn'],p0)
:Et:  = :p: / ( :gamma: - 1.0 )  + .5*:rho:*(:u:*:u: + :v:*:v:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
:dtB: = 0.2* dt.diff(:beta:,:rho:)
:dt: = numpy.minimum(:dt:,:dtB:)
:umag: = sqrt( :u:*:u: + :v:*:v: )
"""
eom = eom.replace("u0", str(u0)).replace("p0", str(p0)).replace("rho0", str(rho0))


# Add the EOM to the solver
ss.EOM(eom)


# Initialize variables
ic = """
:gamma: = 1.4
:R: = 1.0
:cp: = :R: / (1.0 - 1.0/:gamma: )
:cv: = :cp: - :R:
#rad = sqrt( (meshx-pi)**2  +  (meshy-pi)**2 ) 
rad = sqrt( meshx**2  +  meshy**2 ) 
:phi: = rad - pi/4.0
:rho: = 1.0 + 3d()
:p:  =  1.0 + 3d() #exp( -(meshx-1.5)**2/.25**2)*.1
:u: = where( :phi:>0.5, mach * sqrt( :p: / :rho: * :gamma:) , 0.0 )
#:u: = mach * sqrt( :p: / :rho: * :gamma:)
:u: = gbar( gbar( :u: ) )
:v: = 0.0 + 3d()
:Et: = :p:/( :gamma: - 1.0 ) + .5*:rho:*(:u:*:u: + :v:*:v:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
[:gx:,:gy:,:gz:] = grad( :phi: )
:gx: = gbar( :gx: )
:gy: = gbar( :gy: )
"""
ic = ic.replace("mach", str(mach))

# Set the initial conditions
ss.setIC(ic)

# Length scale for art. viscosity
# Initialize variables
x = ss.mesh.coords[0].data
y = ss.mesh.coords[1].data
z = ss.mesh.coords[2].data
dx = x[1, 0, 0] - x[0, 0, 0]
# ss.variables['dx2'].data += dx**2


# Write a time loop
time = 0.0
viz = True

# Approx a max dt and stopping time
tt = 3.0  #

# Mesh for viz on master
xx = ss.PyMPI.zbar(x)
yy = ss.PyMPI.zbar(y)

# Start time loop
cnt = 1
viz_freq = 100
pvar = "umag"

# import pdb
# pdb.set_trace()
CFL = 1.0
dt = ss.variables["dt"].data * CFL * 0.01


v = ss.PyMPI.zbar(ss.variables[pvar].data)
phi = ss.PyMPI.zbar(ss.variables["phi"].data)
if ss.PyMPI.master and (not test):
    plt.figure(2)
    plt.clf()
    plt.contourf(xx, yy, v, 64, cmap=cm.jet)
    plt.contour(xx, yy, phi, [0.0])
    plt.plot(xx, yy, "k-", lw=0.5, alpha=0.5)
    plt.plot(xx.T, yy.T, "k-", lw=0.5, alpha=0.5)
    plt.title(pvar)
    plt.pause(0.001)


# Make a probe set for diagnostics
x = [5.0] * 20
y = numpy.linspace(-10, 10, 20)
probes = pyrandaProbes(ss, x=x, y=y, z=None)


while tt > time:
    # Update the EOM and get next dt
    time = ss.rk4(time, dt)
    dt = min(ss.variables["dt"].data * CFL, dt * 1.1)
    dt = min(dt, (tt - time))

    # Print some output
    ss.iprint("%s -- %s --- %f" % (cnt, time, dt))
    cnt += 1
    if viz and (not test):
        v = ss.PyMPI.zbar(ss.variables[pvar].data)
        phi = ss.PyMPI.zbar(ss.variables["phi"].data)
        if cnt % viz_freq == 1:
            ss.write()

        if (ss.PyMPI.master) and (cnt % viz_freq == 1):  # or True:
            plt.figure(1)
            plt.clf()
            plt.plot(probes.get("p"))

            plt.figure(2)
            plt.clf()
            plt.contourf(xx, yy, v, 64, cmap=cm.jet)
            plt.contour(xx, yy, phi, [0.0])
            # plt.plot(xx, yy, 'k-', lw=0.5, alpha=0.5)
            # plt.plot(xx.T, yy.T, 'k-', lw=0.5, alpha=0.5)
            plt.title(pvar)
            plt.axis("equal")
            plt.pause(0.001)


# ss.writeRestart()#ivars=rvars)
# from pyranda import pyrandaRestart as pr
# tt = pr("cylinder_curvilinear")


# Curve test.  Write file and print its name at the end
if test:
    v = ss.PyMPI.zbar(ss.variables[pvar].data)
    ny = ss.PyMPI.ny
    v1d = v[:, int(ny / 2)]
    x1d = xx[:, int(ny / 2)]
    fname = (testName or problem) + ".dat"
    numpy.savetxt(fname, (x1d, v1d))
    print(fname)
