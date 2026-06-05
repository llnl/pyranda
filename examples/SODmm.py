import numpy
import matplotlib.pyplot as plt
import pdb
import matplotlib
from pyranda import pyrandaSim, pyrandaBC, pyrandaIBM


def move_figure(f, x, y):
    """Move figure's upper left corner to pixel (x, y)"""
    backend = matplotlib.get_backend()
    if backend == "TkAgg":
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == "WXAgg":
        f.canvas.manager.window.SetPosition((x, y))
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        f.canvas.manager.window.move(x, y)
    # plt.show()


def plotFix(plots):
    px_max = 1300
    py_max = 1000
    """
    Arrange the plots in reasonable array
    """
    xw = 640
    yw = 550
    x1 = y1 = 0
    for p in plots:
        # print(x1)
        # print(y1)
        move_figure(p, x1, y1)
        x1 += xw
        if x1 > px_max:
            x1 = 0
            y1 += yw
            if y1 > py_max:
                y1 = 0


## Define a mesh
L = numpy.pi * 2.0
Npts = 500
Lp = L * (Npts - 1.0) / Npts

imesh = """
xdom = (0.0, Lp, Npts)
""".replace("Lp", str(Lp)).replace("Npts", str(Npts))

# Initialize a simulation object on a mesh
ss = pyrandaSim("sodMM", imesh)
ss.addPackage(pyrandaBC(ss))
ss.addPackage(pyrandaIBM(ss))

deltaPhi = 2.0 * L / float(Npts)
# Define the equations of motion
eom = """
# Primary Equations of motion here
ddt(:rhoA:)  =  -ddx(:rhoA:*:uA:)
ddt(:rhoB:)  =  -ddx(:rhoB:*:uB:)
ddt(:rhouA:) =  -ddx(:rhouA:*:uA: + :pA: - :tauA:) - :FA:
ddt(:rhouB:) =  -ddx(:rhouB:*:uB: + :pB: - :tauB:) - :FB:
ddt(:EtA:)   =  -ddx( (:EtA: + :pA: - :tauA:)*:uA: ) - :FA:*:uA:
ddt(:EtB:)   =  -ddx( (:EtB: + :pB: - :tauB:)*:uB: ) - :FB:*:uB:
ddt(:phi:)   = -:gx:*:uphi: - .1 * sign(:phi:)*(:mgp:-1.0) * deltaPhi / (deltat+1.0e-10)
#ddt(:uphi:)  =  -ddx( :uphi:*:uphi: - :tauphi:) - :Fphi:/(:rhoA: + :rhoB:)
# Conservative filter of the EoM
:rhoA:       =  fbar( :rhoA:  )
:rhoB:       =  fbar( :rhoB:  )
:rhouA:      =  fbar( :rhouA: )
:rhouB:      =  fbar( :rhouB: )
:EtA:        =  fbar( :EtA:   )
:EtB:        =  fbar( :EtB:   )
#:phi:        =  .1 * gbar( :phi: ) + :phi: * .9
# Phi gradient
[:gx:,:gy:,:gz:] = grad( :phi: )
# Update the primatives and enforce the EOS
:uA:         =  :rhouA: / :rhoA:
:uB:         =  :rhouB: / :rhoB:
:pA:         =  ( :EtA: - .5*:rhoA:*(:uA:*:uA:) ) * ( :gamma: - 1.0 )
:pB:         =  ( :EtB: - .5*:rhoB:*(:uB:*:uB:) ) * ( :gamma: - 1.0 )
:wgt:        = (1.0-tanh(:phi:/deltaPhi)**2)
:dphi:       = :wgt:*.5/ deltaPhi
:Fphi:       = 1.0*:dphi:*:gx:*(:pA:-:pB:)
:mgp:        =  sqrt( :gx:**2 + :gy:**2 )
:rhophi: = where(:phi: > 0.0, :rhoA:, :rhoB:)
:rhophi:         = ibmS( :rhophi: , :phi:, [:gx:,:gy:,:gz:] )
:rhophi:         = ibmS( :rhophi: , -:phi:, [-:gx:,-:gy:,-:gz:] )
:FA: = :Fphi: * :rhophi: / :rhoA:
:FB: = :Fphi: * :rhophi: / :rhoB:
# Close level set velocity
#:uphi: = (:rhouA: + :rhouB:) / (:rhoA: + :rhoB:)
#:uphi: = where(abs(:phi:) < deltaPhi, (:rhouA: + :rhouB:) / (:rhoA: + :rhoB:), :uphi: )
#:rhouphi: = where(:phi: > 0.0, :rhouA:, :rhouB:)
#:rhouphi:         = ibmS( :rhouphi: , :phi:, [:gx:,:gy:,:gz:] )
#:rhouphi:         = ibmS( :rhouphi: , -:phi:, [-:gx:,-:gy:,-:gz:] )
#:uphi: = :rhouphi: / :rhophi:
:uphi: = where(:phi: > 0.0, :uA:, :uB:)
:uphi:         = ibmS( :uphi: , :phi:, [:gx:,:gy:,:gz:] )
:uphi:         = ibmS( :uphi: , -:phi:, [-:gx:,-:gy:,-:gz:] )
#:uphi:  = gbar( gbar( gbar( :uphi: ) ) )
#:uphi:  = gbar( gbar( gbar( :uphi: ) ) )
#:uphi:  = gbar( gbar( gbar( :uphi: ) ) )
#:uphi:  = gbar( gbar( gbar( :uphi: ) ) )
#:uphi: = :wgt:*:uphi: #+ gbar(:uphi:)*(1.0-:wgt:**2)
#:uphi: = max(:uphi:) + :uphi:*0.0
#:uphi: = where( abs(:phi:) < deltaPhi, :uphi:, gbar(:uphi:) )
#:uphi: = where( abs(:phi:) < deltaPhi, :uphi:, gbar(:uphi:) )
#:uphi: = where( abs(:phi:) < deltaPhi, :uphi:, gbar(:uphi:) )
#:wgt: = .5*(1.0 - tanh( (:phi: - 8.0*deltaPhi)/(deltaPhi) ))
#:wgt2: = .5*(1.0 - tanh( (-:phi: - 8.0*deltaPhi)/(deltaPhi) ))
#:wgt3: = numpy.minimum( :wgt:, :wgt2: )
#:uphi: = gbar( :uphi:*:wgt3: )
# Apply constant BCs
bc.extrap(['rhoA','EtA'],['x1'])
bc.extrap(['rhoB','EtB'],['xn'])
bc.const(['uA'],['x1'],0.0)
#bc.const(['uA'],['xn'],0.0)
bc.const(['uB'],['xn'],0.0)
# IBM
:rhoA:         = ibmS( :rhoA: , :phi:, [:gx:,:gy:,:gz:] )
[:uA:,:v:,:w:] = ibmV( [:uA:,:v:,0.0], :phi:, [:gx:,:gy:,:gz:], [:uphi:,0.0,0.0] )
:pA:          = ibmS( :pA: , :phi:, [:gx:,:gy:,:gz:] )
:EtA:         = :pA: / ( :gamma: - 1.0 ) + .5*:rhoA:*(:uA:*:uA:)
:rhouA:        = :rhoA:*:uA:
:rhoB:         = ibmS( :rhoB: , -:phi:, [-:gx:,-:gy:,-:gz:] )
[:uB:,:v:,:w:] = ibmV( [:uB:,:v:,0.0], -:phi:, [-:gx:,-:gy:,-:gz:], [:uphi:,0.0,0.0] )
:pB:          = ibmS( :pB: , -:phi:, [-:gx:,-:gy:,-:gz:] )
:EtB:         = :pB: / ( :gamma: - 1.0 ) + .5*:rhoB:*(:uB:*:uB:)
:rhouB:        = :rhoB:*:uB:
# Artificial bulk viscosity (old school way)
:divA:       =  ddx(:uA:)
:divB:       =  ddx(:uB:)
#:divphi:     =  ddx(:uphi:)
:betaA:      =  gbar( ring(:divA:) * :rhoA:) * 7.0e-2
:betaB:      =  gbar( ring(:divB:) * :rhoB:) * 7.0e-2
#:betaphi:      =  gbar( ring(:divphi:)) * 7.0e-2
:tauA:       =  :betaA:*:divA:
:tauB:       =  :betaB:*:divB:
#:tauphi:       =  :betaphi:*:divphi:
""".replace("deltaPhi", str(deltaPhi))

# Add the EOM to the solver
ss.EOM(eom)

# Initial conditions SOD shock tube in 1d
ic = """
:gamma: = 1.4 + 3d()
pert    = .0*exp( -abs(meshx-1.0)**2/(.2**2))
:pA:    =  gbar(where(meshx < pi, 10.0, 1.0)) * ( 1.0 + pert)
:rhoA:  = (1.0 + 3d()) * ( 1.0 + pert*(:gamma:-1) )
:uA:    = pert
:rhouA: = :rhoA: * :uA:
:EtA:   = :pA: / ( :gamma: - 1.0 ) + .5*:rhoA:*(:uA:*:uA:)
:pB: = gbar(where(meshx < pi, 10.0, 1.0)) * ( 1.0 + pert)
:EtB:   =  :pB:/(:gamma:-1.0)
:rhoB:  = 0.125 + 3d()
:phi:   = pi - meshx
"""

# Set the initial conditions
ss.setIC(ic)


# Write a time loop
time = 0.0

# Approx a max dt and stopping time
v = 1.0
dt_max = v / ss.mesh.nn[0] * 0.05
tt = L / v * 2.025  # dt_max

# Mesh for viz on master
x = ss.mesh.coords[0].data
xx = ss.PyMPI.zbar(x)

# Start time loop
dt = dt_max
cnt = 1
viz_freq = 30
pvar = "rhoA"
viz = True

pdb.set_trace()


def updatePlots():
    rhoA = ss.PyMPI.zbar(ss.variables["rhoA"].data)
    rhoB = ss.PyMPI.zbar(ss.variables["rhoB"].data)
    pA = ss.PyMPI.zbar(ss.variables["pA"].data)
    pB = ss.PyMPI.zbar(ss.variables["pB"].data)
    uA = ss.PyMPI.zbar(ss.variables["uA"].data)
    uB = ss.PyMPI.zbar(ss.variables["uB"].data)
    phi = ss.PyMPI.zbar(ss.variables["phi"].data)
    uphi = ss.PyMPI.zbar(ss.variables["uphi"].data)

    imin = numpy.argmin(numpy.abs(phi), 0)

    f1 = plt.figure(1)
    plt.clf()
    plt.plot(xx[:, 0], rhoA[:, 0], "k.-")
    plt.plot(xx[:, 0], rhoB[:, 0], "b.-")
    plt.plot([xx[imin, 0], xx[imin, 0]], [rhoB[imin, 0], rhoA[imin, 0]], "r-")
    plt.title("rho")
    #
    f2 = plt.figure(2)
    plt.clf()
    plt.plot(xx[:, 0], pA[:, 0], "k.-")
    plt.plot(xx[:, 0], pB[:, 0], "b.-")
    plt.title("p")
    #
    f3 = plt.figure(3)
    plt.clf()
    plt.plot(xx[:, 0], uA[:, 0], "k.-")
    plt.plot(xx[:, 0], uB[:, 0], "b.-")
    plt.title("u")
    #
    f4 = plt.figure(4)
    plt.clf()
    plt.plot(xx[:, 0], uphi[:, 0], "k.-")
    plt.plot(xx[imin, 0], uphi[imin, 0], "ro")
    plt.plot(xx[:, 0], xx[:, 0] * 0.0, "r--")
    plt.title("U-phi")
    #
    f5 = plt.figure(5)
    plt.clf()
    plt.plot(xx[:, 0], phi[:, 0], "k.-")
    plt.plot(xx[imin, 0], phi[imin, 0], "ro")
    plt.plot(xx[:, 0], xx[:, 0] * 0.0, "r--")
    plt.title("phi")
    #
    # f, ax = plt.subplots()
    # move_figure(f1, 500, 500)
    plotFix([f1, f2, f3, f4, f5])
    plt.pause(0.001)


while tt > time:
    # Update the EOM and get next dt
    time = ss.rk4(time, dt)
    dt = min(dt_max, (tt - time))

    # Print some output
    ss.iprint("%s -- %s" % (cnt, time))
    cnt += 1
    if viz:
        if (ss.PyMPI.master and (cnt % viz_freq == 0)) and True:
            poo = raw_input("fff")  # noqa: F821 (marked as undefined)
            updatePlots()


ss.writeGrid()
ss.write()
