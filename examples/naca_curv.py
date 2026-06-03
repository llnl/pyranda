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
    Npts = 100

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

problem = "naca_test2"

Lp = L * (Npts - 1.0) / Npts


# Make and load an NACA mesh
from omesh import naca_omesh

NACA = "2412"
nx = 300
ny = 100
[xnaca, ynaca] = naca_omesh(
    NACA, nx, ny, PLOT=False, dr0=0.0025, fact=1.05, iS=5, te=0.08, teSig=12, dratio=10
)


def nacaMesh(i, j, k):
    x = xnaca[i, j]
    y = ynaca[i, j]
    z = 0.0
    return x, y, z


mesh_options = {}
mesh_options["coordsys"] = 3
mesh_options["function"] = nacaMesh
mesh_options["periodic"] = numpy.array([True, False, True])
mesh_options["dim"] = 2
mesh_options["x1"] = [-2 * Lp, -2 * Lp, 0.0]
mesh_options["xn"] = [2 * Lp, 2 * Lp, Lp]
mesh_options["nn"] = [nx, ny, 1]


# Initialize a simulation object on a mesh
ss = pyrandaSim(problem, mesh_options)
ss.addPackage(pyrandaBC(ss))
ss.addPackage(pyrandaIBM(ss))
ss.addPackage(pyrandaTimestep(ss))


rho0 = 1.0
p0 = 1.0
gamma = 1.4
mach = 0.7
s0 = numpy.sqrt(p0 / rho0 * gamma)
u0 = s0 * mach
e0 = p0 / (gamma - 1.0) + rho0 * 0.5 * u0 * u0


# Define the equations of motion
eom = """
# Primary Equations of motion here
ddt(:rho:)  =  -div(:rho:*:u:,  :rho:*:v:)
ddt(:rhou:) =  -div(:rhou:*:u: - :tauxx:, :rhou:*:v: - :tauxy:)
ddt(:rhov:) =  -div(:rhov:*:u: - :tauxy:, :rhov:*:v: - :tauyy:)
ddt(:Et:)   =  -div( (:Et: - :tauxx:)*:u: -:tauxy:*:v:- :tx:*:kappa:, (:Et: - :tauyy:)*:v: - :tauxy:*:u: - :ty:*:kappa:)
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
:beta:      =  gbar( ring(:div:) * :rho: ) * 7.0e-2
# Artificial thermal conductivity
[:tx:,:ty:,:tz:] = grad(:T:)
:kappa:     = gbar( ring(:T:)* :rho:*:cv:/(:T: * :dt: ) ) * 1.0e-3
# Artificial shear
:S:    = sqrt( :ux:*:ux: + :vy:*:vy: + .5*((:uy:+:vx:)**2 ) )
:mu:   =  gbar( ringV(:u:,:v:,:w:) * :rho: ) * 1.0e-4 #+ .0005
#:mu:   =  gbar( ring(:S:  ) * :rho: ) * 1.0e-4
# Tau
[:ux:,:uy:,:tz:] = grad(:u:)
[:vx:,:vy:,:tz:] = grad(:v:)
:taudia:    =  (:beta:-2./3.*:mu:) *:div: - :p:
:tauxx:     =  2.0*:mu:*:ux:   + :taudia:
:tauyy:     =  2.0*:mu:*:vy:   + :taudia:
:tauxy:     = :mu:*(:uy:+:vx:)
# Apply constant BCs
bc.extrap(['u','v','rho','p'],['yn'])
bc.const(['u'],['yn'],u0)
bc.extrap(['rho','p'],['y1'])
bc.slip([ ['u','v']  ],['y1'])
# Update the conserved quantities
:Et:  = :p: / ( :gamma: - 1.0 )  + .5*:rho:*(:u:*:u: + :v:*:v:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dtC: = dt.courant(:u:,:v:,:w:,:cs:)
:dtB: = 2.5* dt.diff(:beta:,:rho:)
:dt: = numpy.minimum(:dtC:,:dtB:)
:dtM: = 0.2* dt.diff(:mu:,:rho:)
:dt: = numpy.minimum(:dt:,:dtM:)
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
:rho: = 1.0 + 3d()
:p:  =  1.0 + 3d() 
:u: = mach * sqrt( :p: / :rho: * :gamma:)
bc.slip([ ['u','v']  ],['y1'])
:u: = gbar(gbar(gbar(gbar(gbar(gbar( gbar( :u: ) ) )))))
:v: = gbar(gbar(gbar(gbar(gbar(gbar( gbar( :v: ) ) )))))
:Et: = :p:/( :gamma: - 1.0 ) + .5*:rho:*(:u:*:u: + :v:*:v:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
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
tt = 18.0  #

# Mesh for viz on master
xx = ss.PyMPI.zbar(x)
yy = ss.PyMPI.zbar(y)

# Start time loop
cnt = 1
viz_freq = 500
pvar = "umag"


CFL = 1.0
dt = ss.variables["dt"].data * CFL


v = ss.PyMPI.zbar(ss.variables[pvar].data)
p = ss.PyMPI.zbar(ss.variables["p"].data)

if not test:
    ss.plot.figure(2)
    ss.plot.clf()
    ss.plot.contourf(pvar, 64, cmap=cm.jet)
    # ss.plot.plot(xx, yy, 'k-', lw=0.5, alpha=0.5)
    ss.plot.showGrid()


wvars = ["p", "rho", "u", "v", "beta", "mu", "umag"]

ss.write(wvars)

# import pdb
# pdb.set_trace()

ramp = 0.01
dt *= ramp
while tt > time:
    # Update the EOM and get next dt

    time = ss.rk4(time, dt)
    dt = ss.variables["dt"].data * CFL
    dt = min(dt, (tt - time))

    dt *= ramp
    ramp = min(1.01 * ramp, 1.0)
    # ss.iprint(ramp)

    # airfoil extrap
    tramp1 = 0.0  # .03
    tramp2 = 0.0  # .35
    if time < tramp1:
        ss.parse(":u: = gbar(:u:)")
        ss.parse(":v: = gbar(:v:)")
        ss.parse(":rho: = gbar(:rho:)")
        ss.parse(":rhou: = :rho:*:u:")
        ss.parse(":rhov: = :rho:*:v:")

    if time > tramp1 and time < tramp2:
        wgt = 1.0 - (time - tramp1) / (tramp2 - tramp1)  # Linear ramp
        ss.parse(":u:   = :u:*(1.0-%s)   + %s*gbar(:u:)" % (wgt, wgt))
        ss.parse(":v:   = :v:*(1.0-%s)   + %s*gbar(:v:)" % (wgt, wgt))
        ss.parse(":rho: = :rho:*(1.0-%s) + %s*gbar(:rho:)" % (wgt, wgt))
        ss.parse(":rhou: = :rho:*:u:")
        ss.parse(":rhov: = :rho:*:v:")

    # Print some output
    stab = "None"
    if ss.variables["dt"].data == ss.variables["dtB"].data:
        stab = "bulk"
    if ss.variables["dt"].data == ss.variables["dtM"].data:
        stab = "shear"
    if ss.variables["dt"].data == ss.variables["dtC"].data:
        stab = "CFL"

    ss.iprint("%s -- %f --- %f --- %s" % (cnt, time, dt, stab))
    cnt += 1
    if viz and (not test):
        v = ss.PyMPI.zbar(ss.variables[pvar].data)
        if cnt % viz_freq == 1:  # or True:
            ss.plot.figure(2)
            ss.plot.clf()
            ss.plot.contourf(pvar, 64, cmap=cm.jet)
            ss.plot.showGrid()

            ss.write(wvars)
