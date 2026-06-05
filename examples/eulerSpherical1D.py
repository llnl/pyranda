from __future__ import print_function
import sys
import numpy
import matplotlib.pyplot as plt

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
problem = "sod2dSpherical"
# mesh_options = defaultMeshOptions()
# mesh_options['coordsys'] = 2
# mesh_options['periodic'] = numpy.array([False, False, True])
# mesh_options['dim'] = 2
# mesh_options['x1'] = [ L/5 , numpy.pi/2-.01  ,  0*numpy.pi - spike_angle ]
# mesh_options['xn'] = [ L   , numpy.pi/2+.01  ,  0*numpy.pi + spike_angle ]
# mesh_options['nn'] = [ Npts*4,   1  ,  Npts]
# mesh_options['symmetric'][0][0] = True

# Initialize a simulation object on a mesh
# ss = pyrandaSim(problem,mesh_options)
# ss.addPackage( pyrandaBC(ss) )
# ss.addPackage( pyrandaTimestep(ss) )

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
# bc.const(['u'],['x1','xn'],0.0)
# Compute some max time steps
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
:dt: = numpy.minimum( :dt:, 0.1 * dt.diff(:beta:,:rho:) )
"""

# eom += """# Apply constant BCs
# bc.extrap(['rho','Et'],['x1','xn','z1','zn'],order=1)
# bc.const(['u','w'],['x1','xn','z1','zn'],0.0)
# """

# print(eom)

# Add the EOM to the solver
# ss.EOM(eom)


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
# ss.setIC(ic)

# Write a time loop
time = 0.0
viz = True

# Approx a max dt and stopping time
v = 1.0
tt = 2.0

if test:
    tt = 0.1

pvar = "u"

# import cProfile
# cProfile.run("time = ss.rk4(time,dt)")
# dx = numpy.diff( ss1D.var('meshx')[0:2,0,0] )[0]
# x1 = float( ss1D.var('meshx')[0,0,0] )

# L1 = x1 + dx*7.0
# Mpts = int(L1/dx)
# x0 = (L1/dx - Mpts )*dx
# x0 =

Mpts = Npts * 4
L1 = L
x0 = L / Mpts


### MAKE A R\in (0,Rmin) domain in 1D.
mesh_options1D = defaultMeshOptions()
mesh_options1D["coordsys"] = 2
mesh_options1D["periodic"] = [False, False, True]
# mesh_options1D['symmetric'] = [[True,False],[False,False],[False,False]]
mesh_options1D["x1"] = [x0, numpy.pi / 2 - 0.01, 0 * numpy.pi - spike_angle]
mesh_options1D["xn"] = [L1, numpy.pi / 2 + 0.01, 0 * numpy.pi + spike_angle]
mesh_options1D["nn"] = [Mpts, 1, 1]


# Initialize a simulation object on a mesh
ss1D = pyrandaSim(problem + "_1D", mesh_options1D)
ss1D.addPackage(pyrandaBC(ss1D))
ss1D.addPackage(pyrandaTimestep(ss1D))

eom1D = eom + ":rhow: = :rhow: * 0.0\n"
eom1D += "bc.const(['u'],['x1'],0.0) \n"
eom1D += "bc.extrap(['rho','p'],['x1'],order=1) \n"
eom1D += ":rhou: = :rho:*:u: \n"
eom1D += ":Et: = :p:/( :gamma: - 1.0 ) + .5*:rho:*(:u:*:u: + :w:*:w:) \n"

ss1D.EOM(eom1D)
ss1D.setIC(ic)


# Start time loop
dt_max = v / ss1D.mesh.nn[0] * 0.05
dt = dt_max
cnt = 1
viz_freq = 50
viz_per = tt / 50.0
viz_time = 0.0


while tt > time:
    # Update the EOM and get next dt
    time = ss1D.rk4(time, dt)
    dt = min(dt_max, ss1D.var("dt").data)
    dt = min(dt, (tt - time))

    # Print some output
    ss1D.iprint("%s -- %s" % (cnt, time))
    cnt += 1
    if viz and (not test):
        # if (cnt%viz_freq == 0) and True:
        if time >= viz_time:
            viz_time += viz_per

            fig = plt.figure(1)
            fig.clf()
            ax = fig.subplots(2, 2)

            pvar = "rho"
            x1d = ss1D.plot.getGrid1d("j=0;k=0")
            v1d = ss1D.plot.getLine(ss1D.var(pvar).data, "j=0;k=0")
            ax[0][0].plot(x1d, v1d, "b-x")
            # ax[0][0].title = pvar

            pvar = "p"
            x1d = ss1D.plot.getGrid1d("j=0;k=0")
            v1d = ss1D.plot.getLine(ss1D.var(pvar).data, "j=0;k=0")
            ax[1][0].plot(x1d, v1d, "b-x")
            # ax[1][0].title = pvar

            pvar = "u"
            x1d = ss1D.plot.getGrid1d("j=0;k=0")
            v1d = ss1D.plot.getLine(ss1D.var(pvar).data, "j=0;k=0")
            ax[0][1].plot(x1d, v1d, "b-x")
            # ax[0][1].title = pvar

            pvar = "Et"
            x1d = ss1D.plot.getGrid1d("j=0;k=0")
            v1d = ss1D.plot.getLine(ss1D.var(pvar).data, "j=0;k=0")
            ax[1][1].plot(x1d, v1d, "b-x")
            # ax[1][1].title = pvar

            plt.pause(0.1)
            # ss1D.plot.plot(pvar,style='b-x')

if test:
    pvar = "p"
    p = ss1D.var(pvar)[:, 0, 0]
    x = ss1D.mesh.coordsNative[0][:, 0, 0]
    fname = testName + ".dat"
    numpy.savetxt(fname, (x, p))
    print(fname)
