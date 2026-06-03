import numpy
import matplotlib.pyplot as plt

from pyranda import pyrandaSim, pyrandaBC


## Define a mesh
L = numpy.pi * 2.0
Npts = 200
Lp = L * (Npts - 1.0) / Npts

imesh = """
xdom = (0.0, Lp, Npts)
""".replace("Lp", str(Lp)).replace("Npts", str(Npts))

# Initialize a simulation object on a mesh
ss = pyrandaSim("sod", imesh)
ss.addPackage(pyrandaBC(ss))


# Define the equations of motion
eom = """
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)
ddt(:rhou:) =  -ddx(:rhou:*:u: + :p: - :tau:)
ddt(:Et:)   =  -ddx( (:Et: + :p: - :tau:)*:u: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u:) ) * ( :gamma: - 1.0 )
# Artificial bulk viscosity (old school way)
:div:       =  ddx(:u:) 
:beta:      =  gbar( ring(:div:) * :rho:) * 7.0e-2
:tau:       =  :beta:*:div:
# Apply constant BCs
bc.extrap(['rho','Et'],['x1'])
bc.const(['u'],['x1','xn'],0.0)
"""

# Add the EOM to the solver
ss.EOM(eom)


# Initial conditions SOD shock tube in 1d
ic = """
:gamma: = 1.4
:Et:  = gbar( where( meshx < pi, 1.0/(:gamma:-1.0) , .1 /(:gamma:-1.0) ) )
:rho: = gbar( where( meshx < pi, 1.0    , .125 ) )
"""

# Set the initial conditions
ss.setIC(ic)


# Write a time loop
time = 0.0

# Approx a max dt and stopping time
v = 1.0
dt_max = v / ss.mesh.nn[0] * 0.75
tt = L / v * 0.25  # dt_max

# Start time loop
dt = dt_max
cnt = 1
viz_freq = 25
pvar = "rho"
viz = True

while tt > time:
    # Update the EOM and get next dt
    time = ss.rk4(time, dt)
    dt = min(dt_max, (tt - time))

    # Print some output
    ss.iprint("%s -- %s" % (cnt, time))
    cnt += 1
    if viz:
        if cnt % viz_freq == 0:
            ss.plot.figure(1)
            plt.clf()
            ss.plot.plot(pvar, "b.-")


ss.writeGrid()
ss.write()
