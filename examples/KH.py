# KH.py - Kelvin-Helmholtz paper of McNally et. al.
from pyranda import pyrandaSim, pyrandaTimestep
import sys
import numpy

problem = "KelvinHelmholtz"

# Some global problem definitions
gamma = 5.0 / 3.0  # Ratio of specific heats
p0 = 2.5  # Constant initial pressure value
u1 = 0.5  # Forward velocity
u2 = -0.5  # Backward velocity
rho1 = 1.0  # Density of light
rho2 = 2.0  # Density of heavy
L = 0.025  # Length scale of transition
Vmean = 0.0  # Convective offset
tFinal = 1.5  # Final time
R = 1.0  # Gas constant
test = False  # Run in test mode?

cp = R / (1.0 - 1.0 / gamma)
cv = cp - R

# Grid size
Npts = 128

dt_max = 1.0  # e-3

try:
    Npts = int(sys.argv[1])
except Exception:
    pass

try:
    Vmean = float(sys.argv[2]) * 1.0 / tFinal
except Exception:
    pass

try:
    problem += sys.argv[3]
except Exception:
    pass
try:
    test = bool(int(sys.argv[4]))
except Exception:
    pass


# Make a mesh and a pyrandaSim object
mesh = """
xdom = (0.0, len,  Npts, periodic=True)
ydom = (0.0, len,  Npts, periodic=True)
""".replace("Npts", str(Npts)).replace("len", str(1.0 * (Npts - 1.0) / Npts))

ss = pyrandaSim(problem, mesh)

# Add "Timestep" package for dt.* functions
ss.addPackage(pyrandaTimestep(ss))

# Define the 2D Euler equations with AFLES terms

eom = """
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)            - ddy(:rho:*:v:)               
ddt(:rhou:) =  -ddx(:rhou:*:u: - :tauxx:) - ddy(:rhou:*:v: - :tauxy:)    
ddt(:rhov:) =  -ddx(:rhov:*:u: - :tauxy:) - ddy(:rhov:*:v: - :tauyy:)    
ddt(:Et:)   =  -ddx( (:Et: - :tauxx:)*:u: - :tauxy:*:v:  - :tx:*:kappa: )  - ddy( (:Et: - :tauyy:)*:v: -:tauxy:*:u:- :ty:*:kappa: ) 
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v:) ) * ( gamma - 1.0 )
# Artificial bulk viscosity 
:ux:        =  ddx(:u:)
:vy:        =  ddy(:v:)
:div:       =  :ux: + :vy: 
# Remaining cross derivatives
:uy:        =  ddy(:u:)
:vx:        =  ddx(:v:)
:enst:      = sqrt( (:uy:-:vx:)**2 )
:tke:       = :rho:*(:u:*:u: + :v:*:v: )
:S:         = sqrt( :ux:*:ux: + :vy:*:vy:  + .5*((:uy:+:vx:)**2  ) )
:mu:        =  gbar( abs(ring(:S:  )) ) * :rho: * 1.0e-4
:beta:      =  gbar( abs(ring(:div:)) * :rho: )  * 7.0e-3
:taudia:    =  (:beta:-2./3.*:mu:) *:div: - :p:
:tauxx:     =  2.0*:mu:*:ux:   + :taudia:
:tauyy:     =  2.0*:mu:*:vy:   + :taudia:
:tauxy:     = :mu:*(:uy:+:vx:)
:T:         = :p: / (:rho: * R0 )
[:tx:,:ty:,:tz:] = grad(:T:)
:kappa:     = gbar( ring(:T:)* :rho:*cv/(:T: * :dt: ) ) * 1.0e-3
:cs:  = sqrt( :p: / :rho: * gamma )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)*1.0
:dt: = numpy.minimum(:dt:,0.2 * dt.diff(:beta:,:rho:))
:dt: = numpy.minimum(:dt:,0.2 * dt.diff(:mu:,:rho:))
"""

# Add the EOM to the solver
ss.EOM(eom, {"gamma": gamma, "R0": R, "cv": cv})


# Define the initial conditions here
ic = """
Um = (u1-u2)/2.0
rhoM = (rho1-rho2)/2.0
:u: =                     u1-Um*exp( -(meshy-.75)/L)
:u: = where( meshy < .75, u2+Um*exp( -(.75-meshy)/L) , :u: )
:u: = where( meshy < .50, u2+Um*exp( (-meshy+.25)/L) , :u: )
:u: = where( meshy < .25, u1-Um*exp(  (meshy-.25)/L) , :u: )
:v: = 0.01*sin( 4.*pi*meshx ) + Vmean
:rho: =                     rho1-rhoM*exp( -(meshy-.75)/L)
:rho: = where( meshy < .75, rho2+rhoM*exp( -(.75-meshy)/L) , :rho: )
:rho: = where( meshy < .50, rho2+rhoM*exp( (-meshy+.25)/L) , :rho: )
:rho: = where( meshy < .25, rho1-rhoM*exp(  (meshy-.25)/L) , :rho: )
:p: += p0
# Form conserved variables
:Et: = :p:/( gamma - 1.0 ) + .5*:rho:*(:u:*:u: + :v:*:v:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:cs:  = sqrt( :p: / :rho: * gamma )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)*.1
"""
icDict = {
    "gamma": gamma,
    "u1": u1,
    "u2": u2,
    "p0": p0,
    "rho1": rho1,
    "rho2": rho2,
    "L": L,
    "Vmean": Vmean,
}

ss.setIC(ic, icDict)


# Main time stepping loop
viz_freq = 0.1

# Variables to write to viz.
wvars = ["rho", "u", "v", "p", "beta", "kappa", "mu", "enst"]


time = 0.0
dt = ss.var("dt").data
viz_dump = viz_freq

if not test:
    ss.write(wvars)

while tFinal > time:
    time = ss.rk4(time, dt)
    dt = min(ss.variables["dt"].data, 1.1 * dt)
    dt = min(dt_max, dt)
    dt = min(dt, (tFinal - time))

    # Simulation heart-beat
    ss.iprint("Cycle: %5d --- Time: %10.4e --- deltat: %10.4e" % (ss.cycle, time, dt))

    # Constant time
    if (time > viz_dump) and (not test):
        ss.write(wvars)
        viz_dump += viz_freq

        # ss.plot.figure(1)
        # ss.plot.clf()
        # ss.plot.contourf('rho',64 )
if not test:
    ss.write(wvars)
    ss.plot.figure(1)
    ss.plot.clf()
    ss.plot.contourf("rho", 64)

# Curve test.  Write file and print its name at the end
if test:
    ny = ss.PyMPI.ny
    x = ss.mesh.coords[0]
    rho = ss.variables["rho"]
    islice = "j=%s;k=0" % (int(3 * ny / 4))
    x1d = ss.plot.getLine(x.data, islice)
    v1d = ss.plot.getLine(rho.data, islice)
    fname = problem + ".dat"
    numpy.savetxt(fname, (x1d, v1d))
    print(fname)
