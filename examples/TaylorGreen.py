from __future__ import print_function
import sys
import numpy
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaTimestep


# Try to get args for testing
try:
    Npts = int(sys.argv[1])
except Exception:
    Npts = 32

try:
    test = bool(int(sys.argv[2]))
except Exception:
    test = False

problem = "TGvortex"

## Define a mesh
# Npts = 32
imesh = (
    """
xdom = (0.0, 2*pi*FF,  Npts, periodic=True)
ydom = (0.0, 2*pi*FF,  Npts, periodic=True)
zdom = (0.0, 2*pi*FF,  Npts, periodic=True)
""".replace("Npts", str(Npts))
    .replace("pi", str(numpy.pi))
    .replace("FF", str(float(Npts - 1) / Npts))
)


# Initialize a simulation object on a mesh
ss = pyrandaSim(problem, imesh)
ss.addPackage(pyrandaTimestep(ss))


# Define the equations of motion
eom = """
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)            - ddy(:rho:*:v:)               - ddz(:rho:*:w:)
ddt(:rhou:) =  -ddx(:rhou:*:u: - :tauxx:) - ddy(:rhou:*:v: - :tauxy:)    - ddz(:rhou:*:w: - :tauxz:)
ddt(:rhov:) =  -ddx(:rhov:*:u: - :tauxy:) - ddy(:rhov:*:v: - :tauyy:)    - ddz(:rhov:*:w: - :tauyz:)
ddt(:rhow:) =  -ddx(:rhow:*:u: - :tauxz:) - ddy(:rhow:*:v: - :tauyz:)    - ddz(:rhow:*:w: - :tauzz:)
ddt(:Et:)   =  -ddx( (:Et: - :tauxx:)*:u: - :tauxy:*:v: - :tauxz:*:w: )  - ddy( (:Et: - :tauyy:)*:v: -:tauxy:*:u: - :tauyz:*:w:) - ddz( (:Et: - :tauzz:)*:w: - :tauxz:*:u: - :tauyz:*:v: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:rhow:      =  fbar( :rhow: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:w:         =  :rhow: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v: + :w:*:w:) ) * ( :gamma: - 1.0 )
# Artificial bulk viscosity 
:ux:        =  ddx(:u:)
:vy:        =  ddy(:v:)
:wz:        =  ddz(:w:)
:div:       =  :ux: + :vy: + :wz:
# Remaining cross derivatives
:uy:        =  ddy(:u:)
:uz:        =  ddz(:u:)
:vx:        =  ddx(:v:)
:vz:        =  ddz(:v:)
:wy:        =  ddy(:w:)
:wx:        =  ddx(:w:)
:enst:      = sqrt( (:uy:-:vx:)**2 + (:uz: - :wx:)**2 + (:vz:-:wy:)**2 )
:tke:       = :rho:*(:u:*:u: + :v:*:v: + :w:*:w:)
:S:         = sqrt( :ux:*:ux: + :vy:*:vy: + :wz:*:wz: + .5*((:uy:+:vx:)**2 + (:uz: + :wx:)**2 + (:vz:+:wy:)**2) )
:mu:        =  gbar( abs(ring(:S:  )) ) * :rho: * 1.0e-4
:beta:      =  gbar( abs(ring(:div:)) * :rho: )  * 7.0e-3
:taudia:    =  (:beta:-2./3.*:mu:) *:div: - :p:
:tauxx:     =  2.0*:mu:*:ux:   + :taudia:
:tauyy:     =  2.0*:mu:*:vy:   + :taudia:
:tauzz:     =  2.0*:mu:*:wz:   + :taudia:
:tauxy:     = :mu:*(:uy:+:vx:)
:tauxz:     = :mu:*(:uz:+:wx:)
:tauyz:     = :mu:*(:vz:+:wz:)
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)*1.0
:dt: = numpy.minimum(:dt:,0.2 * dt.diff(:beta:,:rho:))
:dt: = numpy.minimum(:dt:,0.2 * dt.diff(:mu:,:rho:))
"""

# Add the EOM to the solver
ss.EOM(eom)


# Initialize variables
ic = """
:gamma: = 1.4
u0 = 1.0
p0 = 100.0
rho0 = 1.0
L = 1.0
:u: =  u0*sin(meshx/L)*cos(meshy/L)*cos(meshz/L)
:v: = -u0*cos(meshx/L)*sin(meshy/L)*cos(meshz/L)
:w: = 0.0*:u:
:p:  = p0 + rho0/16.0*( ( cos(2.*meshx/L) + cos(2.*meshy/L) ) * ( cos(2.*meshz/L) + 2.0 ) - 2.0 )
:rho: = rho0 + 0.0*:u:
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:rhow: = :rho:*:w:
:Et:  = :p: / (:gamma:-1.0) + 0.5*:rho:*(:u:*:u: + :v:*:v: + :w:*:w:)
:cs:  = sqrt( :p: / :rho: * :gamma: )
:tke: = :rho:*(:u:*:u: + :v:*:v: + :w:*:w:)
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
"""

# Set the initial conditions
ss.setIC(ic)

# Length scale for art. viscosity
# Initialize variables
x = ss.mesh.coords[0].data
y = ss.mesh.coords[1].data
z = ss.mesh.coords[2].data

# Write a time loop
time = 0.0
viz = True

# Start time loop
CFL = 0.5
dt = ss.variables["dt"].data * CFL

# Viz
cnt = 1
viz_freq = 20
pvar = "u"

tke0 = ss.var("tke").sum()
enst0 = ss.var("enst").sum()
TKE = []
ENST = []
TIME = []

tstop = 25.0
if test:
    tstop = 0.1

while time < tstop:
    # Update the EOM and get next dt
    time = ss.rk4(time, dt)
    dt = ss.variables["dt"].data * CFL

    # Print some output
    tke = ss.var("tke").sum() / tke0
    enst = ss.var("enst").sum() / enst0

    TIME.append(time)
    ENST.append(enst)
    TKE.append(tke)

    ss.iprint("%s -- %s --- TKE: %s " % (cnt, time, tke))
    cnt += 1
    if viz:
        # Write viz data and plot
        if cnt % viz_freq == 0:
            ss.write()

            ss.plot.figure(2)
            ss.plot.clf()
            ss.plot.contourf(pvar, 64, slice3d="k=16", cmap=cm.jet)
            ss.plot.title(pvar + ",Time=%f" % time)


data = ss.PyMPI.subsum3xz(ss.mesh.coords[1].data) / (Npts * Npts)
if ss.PyMPI.master:
    print(data)

if test:
    print(enst)
else:
    if ss.PyMPI.master:
        plt.figure()
        plt.plot(TIME, TKE, "k--")
        plt.plot(TIME, ENST, "b--")
        plt.show()
