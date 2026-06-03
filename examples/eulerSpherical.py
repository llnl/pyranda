from __future__ import print_function
import re
import sys
import time
import numpy
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep
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
spike_angle = numpy.pi / 16
problem = "2dSpherical"
mesh_options = defaultMeshOptions()
mesh_options["coordsys"] = 2
mesh_options["periodic"] = numpy.array([False, False, True])
# mesh_options['dim'] = 2
mesh_options["x1"] = [L / 10, numpy.pi / 2 - 0.01, 0 * numpy.pi - spike_angle]
mesh_options["xn"] = [L, numpy.pi / 2 + 0.01, 0 * numpy.pi + spike_angle]
mesh_options["nn"] = [Npts * 4, 1, Npts]
# mesh_options['symmetric'][0][0] = True

# Initialize a simulation object on a mesh
ss = pyrandaSim(problem, mesh_options)
ss.addPackage(pyrandaBC(ss))
ss.addPackage(pyrandaTimestep(ss))

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
[:fxx:,:fxy:,:fxz:] = [:rhou:*:u: + :p: - :tau:, :v:       , :rhou:*:w:]
[:fyx:,:fyy:,:fyz:] = [      :v:               , :p:-:tau: , :v:       ]
[:fzx:,:fzy:,:fzz:] = [:rhow:*:u:,               :v:       , :rhow:*:w: + :p: - :tau:]
[:xmom:,:ymom:,:zmom:] = divT(:fxx:,:fxy:,:fxz:,:fyx:,:fyy:,:fyz:,:fzx:,:fzy:,:fzz:)
# Boundary conditions
bc.const(['u'],['x1','xn'],0.0)
# Compute some max time steps
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
"""

# eom += """# Apply constant BCs
# bc.extrap(['rho','Et'],['x1','xn','z1','zn'],order=1)
# bc.const(['u','w'],['x1','xn','z1','zn'],0.0)
# """

print(eom)

# Add the EOM to the solver
ss.EOM(eom)


# Initialize variables
ic = "rad = sqrt( (meshx-pi)**2  +  (meshy-0.0)**2 ) "


# SOD shock tube in 1d and 2d
ic += """
:gamma: = 1.4
wgt = .5*(1.0-tanh( (rad-pi/16.0)/.1) )   # [0,1]
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
dt_max = v / ss.mesh.nn[0] * 0.05
tt = 2.0

if test:
    tt = 0.1


# Start time loop
dt = dt_max
cnt = 1
viz_freq = 50
viz_per = tt / 20.0
viz_time = 0.0

pvar = "rho"


while tt > time:
    # Update the EOM and get next dt
    time = ss.rk4(time, dt)
    dt = min(dt_max, ss.var("dt").data)
    dt = min(dt, (tt - time))

    # Print some output
    ss.iprint("%s -- %s" % (cnt, time))
    cnt += 1
    if viz and (not test):
        # if (cnt%viz_freq == 0) and True:
        if time >= viz_time:
            viz_time += viz_per
            ss.write()
            ss.plot.figure(1)
            ss.plot.clf()
            # ss.plot.plot(pvar,slice2d="j=0;k=50",style='k.-')
            # ss.plot.figure(2)
            # ss.plot.clf()
            # ss.plot.contourf(pvar,64)

            y1d = ss.var(pvar)[:, 0, :].mean(axis=1)
            x1d = ss.mesh.coordsNative[0][:, 0, int(Npts / 2)]
            if ss.PyMPI.master:
                plt.plot(x1d, y1d, "k-o")

            ss.plot.figure(2)
            ss.plot.clf()
            x2d = ss.var("meshx")[:, 0, :]
            y2d = ss.var("meshy")[:, 0, :]
            v2d = ss.var(pvar)[:, 0, :]
            if ss.PyMPI.master:
                plt.contourf(x2d, y2d, v2d, 64)
                plt.axis("equal")
                plt.pause(0.1)


if test:
    p = ss.var(pvar)[:, 0, :].mean(axis=1)
    x = ss.mesh.coordsNative[0][:, 0, int(Npts / 2)]
    fname = testName + ".dat"
    numpy.savetxt(fname, (x, p))
    print(fname)
