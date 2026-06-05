import numpy

from pyranda import pyrandaSim, pyrandaBC, pyrandaLEOS, pyrandaTimestep

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
ss.addPackage(pyrandaTimestep(ss))


# LEOS package and materials
eos = pyrandaLEOS(ss)
matID = 70  # 70 - Nitrogen
eos.addMaterial(matID)
ss.addPackage(eos)

rho0 = eos.materials[matID]["rho0"]
T0 = 300.0

##  Look up pressure at these rho/T
P0 = eos.materials[matID]["Pt"].eval(rho0, T0, N=1)[0]
E0 = eos.materials[matID]["Et"].eval(rho0, T0, N=1)[0]

parms = {}
parms["rho0"] = rho0
parms["E0"] = E0
parms["P0"] = P0
parms["T0"] = T0
parms["matID"] = matID

# Define the equations of motion
eom = """
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)
ddt(:rhou:) =  -ddx(:rhou:*:u: + :p: - :tau:)
ddt(:Et:)   =  -ddx( (:Et: + :p: - :tau:)*:u: - :tx:*:kappa:)
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:ke:        = .5*:rho:*:u:*:u:
:ie:        = (:Et: - :ke:)/:rho:
#:p:        =  ( :Et: - .5*:rho:*(:u:*:u:) ) * ( :gamma: - 1.0 )
#:T:         = LEOS.temperature(matID,:rho:,:ie:)
:T:         = LEOS.inv.energy(matID,:rho:,:ie:)
:p:         = LEOS.pressure(matID,:rho:,:T:)
# Artificial bulk viscosity (old school way)
:div:       =  ddx(:u:) 
:beta:      =  gbar( ring(:div:) * :rho:) * 7.0e-1
:tau:       =  :beta:*:div:
# Artificial thermal conductivity
[:tx:,:ty:,:tz:] = grad(:ie:)
:kappa:     = 1.0e-3 * :rho: * gbar( ring(:ie:) / (abs(:ie:)+1.0e-6 )) / :dt:
# Apply constant BCs
bc.extrap(['rho','Et'],['x1'])
bc.const(['u'],['x1','xn'],0.0)
# Time step
:cs: = LEOS.soundSpeed( matID, :rho:, :T: ) 
:dt: = dt.courant(:u:,:v:,:w:,:cs:)*0.1
:dt: = numpy.minimum(:dt:,0.01 * dt.diff(:beta: ,:rho:))
#:dt: = numpy.minimum(:dt:,0.1 * dt.diff(:kappa:,:rho:))
"""

# Add the EOM to the solver
ss.EOM(eom, eomDict=parms)


# Initial conditions SOD shock tube in 1d
ic = """
:gamma: = 1.4
:p:   = gbar( gbar( where( meshx < pi, P0  *10.    , P0   ) ))
:rho: = gbar( gbar( where( meshx < pi, rho0*2.0    , rho0 ) ))

:T:   = LEOS.inv.pressure(matID,:rho:,:p:)
:ie:  = LEOS.energy(matID,:rho:,:T:)
:Et:  = :rho: * :ie:
"""

# Set the initial conditions
ss.var("dt").data = 1.0e-10
ss.setIC(ic, icDict=parms)


# Write a time loop
time = 0.0

# Approx a max dt and stopping time
dt_max = ss.var("dt").data * 1.0

tt = 5.0e-6  # 1000 * dt_max

# Start time loop
dt = dt_max / 10.0
cnt = 1
viz_freq = 100
pvar = "rho"
viz = True

cplots = ["rho", "p", "u", "T", "ie"]


def doPlots():
    ifig = 1  # Plot figure counter
    for cc in cplots:  # Loops over plots
        ss.plot.figure(ifig)  #   - make a new figure
        ss.plot.clf()  #   - clear the figure
        ss.plot.plot(cc, "k-o")  #   - plot a filled contour
        ifig += 1


doPlots()

while tt > time:
    # Update the EOM and get next dt
    time = ss.rk4(time, dt)
    dt = min(dt_max, (tt - time))

    # Print some output
    ss.iprint("%s -- %s" % (cnt, time))
    cnt += 1
    if viz:
        if cnt % viz_freq == 0:
            doPlots()


doPlots()
ss.writeGrid()
ss.write()
