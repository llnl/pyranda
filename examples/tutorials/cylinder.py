import sys
import numpy
from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep, pyrandaIBM

sys.path.append("../")
from meshTest import zoomMesh_solve


## Define a mesh
Npts = 64
L = numpy.pi * 2.0
dim = 2
gamma = 1.4

problem = "cylinder"

Lp = L * (Npts - 1.0) / Npts

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
mesh_options["gridPeriodic"] = numpy.array([False, False, False])
mesh_options["dim"] = 3
mesh_options["x1"] = [-2 * Lp, -2 * Lp, 0.0]
mesh_options["xn"] = [2 * Lp, 2 * Lp, Lp]
mesh_options["nn"] = [Npts, Npts, 1]


# Initialize a simulation object on a mesh
ss = pyrandaSim(problem, mesh_options)
ss.addPackage(pyrandaBC(ss))
ss.addPackage(pyrandaIBM(ss))
ss.addPackage(pyrandaTimestep(ss))


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
# Apply BCs and internal IBM
[:u:,:v:,:w:] = ibmV( [:u:,:v:,0.0], :phi:, [:gx:,:gy:,:gz:], [:u1:,:u2:,0.0] )
:rho: = ibmS( :rho: , :phi:, [:gx:,:gy:,:gz:] )
:p:   = ibmS( :p:   , :phi:, [:gx:,:gy:,:gz:] )
bc.extrap(['rho','p','u'],['xn'])
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


# Write a time loop
time = 0.0
viz = True

# Approx a max dt and stopping time
tt = 3.0  #

# Start time loop
cnt = 1
viz_freq = 100
pvar = "umag"

CFL = 1.0
dt = ss.var("dt").data * CFL * 0.01
wvars = ["p", "rho", "u", "v", "phi"]
ss.write(wvars)
while tt > time:
    # Update the EOM and get next dt
    time = ss.rk4(time, dt)
    dt = min(ss.variables["dt"].data * CFL, dt * 1.1)
    dt = min(dt, (tt - time))

    # Print some output
    ss.iprint("%s -- %s --- %f" % (cnt, time, dt))
    cnt += 1
    if cnt % viz_freq == 1:
        ss.write(wvars)


ss.writeRestart()
