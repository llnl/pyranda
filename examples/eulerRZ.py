from __future__ import print_function
import re
import sys
import time
import numpy
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC
from pyranda.pyrandaMesh import defaultMeshOptions


# Try to get args for testing
try:
    Npts = int(sys.argv[1])
except Exception:
    Npts = 128

try:
    test = bool(int(sys.argv[2]))
except Exception:
    test = False

try:
    testName = sys.argv[3]
except Exception:
    testName = None


## Define a mesh
L = numpy.pi * 2.0
dx2 = L / (Npts - 1) / 2.0
gamma = 1.4
problem = "sod2dRZ"
mesh_options = defaultMeshOptions()
mesh_options["coordsys"] = 1
mesh_options["periodic"] = numpy.array([False, False, True])
# mesh_options['dim'] = 2
mesh_options["x1"] = [dx2, -0.1, 0.0]
mesh_options["xn"] = [L, 0.1, L]
mesh_options["nn"] = [Npts, 1, Npts]
mesh_options["symmetric"][0][0] = True

# Initialize a simulation object on a mesh
ss = pyrandaSim("advection", mesh_options)
ss.addPackage(pyrandaBC(ss))

# Define the equations of motion
eom = """
# Primary Equations of motion here
ddt(:rho:)  =  -:mass:
ddt(:rhou:) =  -:xmom:
ddt(:rhow:) =  -:zmom:
ddt(:Et:)   =  -div( (:Et: + :p: - :tau:)*:u:,:v:, (:Et: + :p: - :tau:)*:w: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhow:      =  fbar( :rhow: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:w:         =  :rhow: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :w:*:w:) ) * ( :gamma: - 1.0 )
# Artificial bulk viscosity (old school way)
:div:       =  div(:u:,:v:,:w:)
:beta:      =  gbar(abs(ring(:div:))) * :rho: * 7.0e-2
:tau:       =  :beta:*:div:
:v:         = :v: * 0.0
# Compute divT of momentum
:mass: = div(:rho:*:u:               , :v:, :rho:*:w:)
[:fxx:,:fxy:,:fxz:] = [:rhou:*:u: + :p: - :tau:, :v:, :rhou:*:w:]
[:fyx:,:fyy:,:fyz:] = [      :v:               , :p:-:tau: , :v:       ]
[:fzx:,:fzy:,:fzz:] = [:rhow:*:u:,               :v:, :rhow:*:w: + :p: - :tau:]
[:xmom:,:ymom:,:zmom:] = divT(:fxx:,:fxy:,:fxz:,:fyx:,:fyy:,:fyz:,:fzx:,:fzy:,:fzz:)
"""

# eom += """# Apply constant BCs
# bc.extrap(['rho','Et'],['x1','xn','z1','zn'],order=1)
# bc.const(['u','w'],['x1','xn','z1','zn'],0.0)
# """

print(eom)

# Add the EOM to the solver
ss.EOM(eom)


# Initialize variables
ic = "rad = sqrt( (meshx-pi)**2  +  (meshz-pi)**2 ) "


# SOD shock tube in 1d and 2d
ic += """
:gamma: = 1.4
wgt = .5*(1.0-tanh( (rad-pi/2.0)/.1) )   # [0,1]
:Et:  = 1.0/(:gamma:-1.0) * (wgt*.9 + .1)
#:Et:  = gbar( where( rad < pi/2.0, 1.0/(:gamma:-1.0) , .1 /(:gamma:-1.0) ) )
#:rho: = gbar( where( rad < pi/2.0, 1.0    , .125 ) )
:rho: = 1.0*wgt + (1.0-wgt)*.125
"""

# Set the initial conditions
ss.setIC(ic)

# Write a time loop
time = 0.0
viz = True

# Approx a max dt and stopping time
v = 1.0
dt_max = v / ss.mesh.nn[0] * 0.75
tt = 2.0

if test:
    tt = 0.5

# Start time loop
dt = dt_max
cnt = 1
viz_freq = 10


pvar = "rho"

while tt > time:
    # Update the EOM and get next dt
    time = ss.rk4(time, dt)
    dt = min(dt_max, (tt - time))

    # Print some output
    ss.iprint("%s -- %s" % (cnt, time))
    cnt += 1
    if viz and (not test):
        if (cnt % viz_freq == 0) and True:
            ss.plot.figure(1)
            ss.plot.clf()
            ss.plot.plot(pvar, slice2d="j=0;k=50", style="k.-")
            ss.plot.figure(2)
            ss.plot.clf()
            ss.plot.contourf(pvar, 64)


if test:
    half = int(Npts / 2.0)
    x = ss.mesh.coords[0][:, 0, 0]
    p = ss.var("p")[half, 0, :]
    fname = testName + ".dat"
    numpy.savetxt(fname, (x, p))
    print(fname)
