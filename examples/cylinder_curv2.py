from __future__ import print_function
import sys
import time
import numpy
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep, pyrandaIBM


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

problem = "cylinder_omesh"

Lp = L * (Npts - 1.0) / Npts

Ri = 1.0
Rf = 10.0

NX = Npts
NY = Npts


def cylMesh(i, j, k):
    theta = float(i) / float(NX) * 2.0 * numpy.pi
    r = float(j) / float(NY - 1) * (Rf - Ri) + Ri
    x = r * numpy.cos(theta)
    y = r * numpy.sin(theta)
    z = 0.0
    return x, y, z


mesh_options = {}
mesh_options["coordsys"] = 3
mesh_options["function"] = cylMesh
mesh_options["periodic"] = numpy.array([True, False, False])
mesh_options["periodicGrid"] = False
mesh_options["dim"] = 3
mesh_options["x1"] = [-2 * Lp, -2 * Lp, 0.0]
mesh_options["xn"] = [2 * Lp, 2 * Lp, Lp]
mesh_options["nn"] = [NX, NY, 1]


# Initialize a simulation object on a mesh
ss = pyrandaSim(problem, mesh_options)
ss.addPackage(pyrandaBC(ss))
ss.addPackage(pyrandaIBM(ss))
ss.addPackage(pyrandaTimestep(ss))


rho0 = 1.0
p0 = 1.0
gamma = 1.4
mach = 1.5
s0 = numpy.sqrt(p0 / rho0 * gamma)
u0 = s0 * mach
e0 = p0 / (gamma - 1.0) + rho0 * 0.5 * u0 * u0


# Define the equations of motion
eom = """
# Primary Equations of motion here
ddt(:rho:)  =  -div(:rho:*:u:,  :rho:*:v:)
ddt(:rhou:) =  -div(:rhou:*:u: + :p: - :tau:, :rhou:*:v:)
ddt(:rhov:) =  -div(:rhov:*:u:, :rhov:*:v: + :p: - :tau:)
ddt(:Et:)   =  -div( (:Et: + :p: - :tau:)*:u: , (:Et: + :p: - :tau:)*:v:  )
# kappa terms - :tx:*:kappa:, - :ty:*:kappa:
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v:) ) * ( :gamma: - 1.0 )
#:T:         = :p: / (:rho: * :R: )
# Artificial bulk viscosity (old school way)
:div:       =  div(:u:,:v:)
:beta:      =  gbar( ring(:div:) * :rho: ) * 7.0e-2
:tau:       = :beta:*:div:
#[:tx:,:ty:,:tz:] = grad(:T:)
#:kappa:     = gbar( ring(:T:)* :rho:*:cv:/(:T: * :dt: ) ) * 1.0e-3
bc.extrap(['u','v','rho','p'],['yn'])
bc.const(['u'],['yn'],u0)
bc.const(['v'],['yn'],0.0)
bc.extrap(['rho','p'],['y1'])
bc.slip([ ['u','v']  ],['y1'])
:Et:  = :p: / ( :gamma: - 1.0 )  + .5*:rho:*(:u:*:u: + :v:*:v:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
:dtB: = dt.diff(:beta:,:rho:)
#:dtK: = dt.diff(:kappa:,:rho:*:R:)
#:dt: = numpy.minimum(:dt:,:dtB:)
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
#rad = sqrt( (meshx-numpy.pi)**2  +  (meshy-numpy.pi)**2 ) 
rad = sqrt( meshx**2  +  meshy**2 ) 
:rho: = 1.0 + 3d()
:p:  =  1.0 + 3d() #exp( -(meshx-1.5)**2/.25**2)*.1
:u: = mach * sqrt( :p: / :rho: * :gamma:)
:v: = 0.0 + 3d()
:Et: = :p:/( :gamma: - 1.0 ) + .5*:rho:*(:u:*:u: + :v:*:v:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
"""
ic = ic.replace("mach", str(mach))


# Set the initial conditions
ss.setIC(ic)

# Write a time loop
time = 0.0
viz = True

# Approx a max dt and stopping time
tt = 3.0  #

# Start time loop
cnt = 1
viz_freq = 100
pvar = "umag"

CFL = 0.8
dt = ss.variables["dt"].data * CFL * 0.1


if not test:
    ss.plot.figure(2)
    ss.plot.clf()
    ss.plot.contourf(pvar, 64, cmap=cm.jet)
    ss.plot.contour("p", 8, colors="black")
    # ss.plot.showGrid()
    ss.plot.title(pvar)


while tt > time:
    # Update the EOM and get next dt
    time = ss.rk4(time, dt)

    # Time step controls
    dt_msg = "cfl"
    dt = ss.variables["dt"].data * CFL
    dtt = min(dt, 1.1 * dt)
    if dtt != dt:
        dt_msg = "ramp"
    dt = dtt
    dtt = min(0.2 * ss.variables["dtB"].data, dt)
    if dtt != dt:
        dt_msg = "bulk"
    dt = dtt
    # dtt = min( 0.2 * ss.variables['dtK'].data, dt)
    # if ( dtt != dt ):
    #    dt_msg = 'kappa'
    # dt = dtt
    dtt = min(dt, (tt - time))
    if dtt != dt:
        dt_msg = "vis-dump"
    dt = dtt

    # Print some output
    ss.iprint(
        "Cycle:%s -- Time:%.4e --- dt:%.4e --- Limit:%s" % (cnt, time, dt, dt_msg)
    )
    cnt += 1
    if viz and (not test):
        if cnt % viz_freq == 1:
            ss.write(["beta", "u", "v", "rho", "umag"])

            ss.plot.figure(2)
            ss.plot.clf()
            ss.plot.contourf(pvar, 64, cmap=cm.jet)
            ss.plot.contour("p", 8, colors="black")
            # ss.plot.showGrid()
            ss.plot.title(pvar)


# Curve test.  Write file and print its name at the end
if test:
    # Initialize variables
    x = ss.mesh.coords[0].data
    y = ss.mesh.coords[1].data

    # Mesh for viz on master
    xx = ss.PyMPI.zbar(x)
    yy = ss.PyMPI.zbar(y)

    v = ss.PyMPI.zbar(ss.variables[pvar].data)
    ny = ss.PyMPI.ny
    v1d = v[:, int(ny / 2)]
    x1d = xx[:, int(ny / 2)]
    fname = testName + ".dat"
    numpy.savetxt(fname, (x1d, v1d))
    print(fname)
