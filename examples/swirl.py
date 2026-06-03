import sys
import time
import numpy
import matplotlib.pyplot as plt

from pyranda import pyrandaSim, pyrandaIBM


# Try to get args
try:
    Npts = int(sys.argv[1])
except Exception:
    Npts = 128

# import pdb
# pdb.set_trace()
try:
    test = bool(int(sys.argv[2]))
except Exception:
    test = False

try:
    diffusive = bool(int(sys.argv[3]))
except Exception:
    diffusive = False


# print sys.argv[2]

## Define a mesh
L = numpy.pi * 2.0
Lp = L * (Npts - 1.0) / Npts
twoD = True

mesh_options = {}
mesh_options["coordsys"] = 0
mesh_options["periodic"] = numpy.array([twoD, twoD, True])
mesh_options["dim"] = 3
mesh_options["x1"] = [0.0, 0.0, 0.0]
mesh_options["xn"] = [Lp, Lp, Lp]
# mesh_options['x1'] = [ -Lp/2.0 , -Lp/2.0  ,  0.0 ]
# mesh_options['xn'] = [ Lp/2.0   , Lp/2.0   ,  Lp ]

mesh_options["nn"] = [Npts, 1, 1]
if twoD:
    mesh_options["nn"] = [Npts, Npts, 1]


# Initialize a simulation object on a mesh
ss = pyrandaSim("swirl", mesh_options)
ss.addPackage(pyrandaIBM(ss))


# Define the equations of motion
eom = """
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u: - :tau:) - ddy(:rho:*:v: - :tau:) 
ddt(:rhoA:)  =  -ddx(:rhoA:*:uA: - :tauA:) -ddy(:rhoA:*:vA: - :tauA:) 
ddt(:phi:)  =  - :gx: * :u1: - :gy: * :v1:  - .01*sign(:phi:)*(:mgp:-1.0)
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhoA:      =  fbar( :rhoA:  )
#:phi:       =  fbar( :phi: )
# Update the primatives and enforce the EOS
:dr:        = ddx(:rho:)
:drA:       = ddx(:rhoA:)
:tau:       = 5.0e-7*gbar(abs(lap(lap(:rho:)))) * :dx4: * :u: * :dr: * :dx:
:tauA:      = 5.0e-7*gbar(abs(lap(lap(:rhoA:)))) * :dx4: * :uA: * :drA: * :dx:
:u:         = -sin(meshx/2.0)**2 * sin(meshy)*cos(simtime) * 2.0*pi
:uA:        = :u:
:u1:        = :u:
:v:         = sin(meshy/2.0)**2 * sin(meshx)*cos(simtime) * 2.0*pi
:vA:        = :v:
:v1:        = :v:
# Immersed boundary method used for MM levelset
[:gx:,:gy:,:gz:] = grad( :phi: )
:mgp:  =  sqrt( :gx:**2 + :gy:**2 )
#[:u:,:v:,:w:] = ibmV( [:u:,:v:,0.0], :phi:, [:gx:,:gy:,:gz:], [:u1:,0.0,0.0] )
:rho:         = ibmS( :rho: , :phi:, [:gx:,:gy:,:gz:] )
#[:uA:,:v:,:w:] = ibmV( [:uA:,:v:,0.0], -:phi:, [-:gx:,-:gy:,-:gz:], [:u1:,0.0,0.0] )
:rhoA:         = ibmS( :rhoA: , -:phi:, [-:gx:,-:gy:,-:gz:] )
:uT:        =  where( :phi: > 0.0, :u:, :uA:)
#:rhoB:        =  where( :phi: > 0.0, :rho:, :rhoA:)
#:rhoT:        =  where( abs(:phi:/:dx:*2.0) > 5.0, :rhoB:, .5*(:phi:/:dx:*2.0 + 1.0)*(:rho:-:rhoA:) + (:rhoA:)   )
:rhoT:       =  0.5 * ( tanh( :phi:/(:dx:*1.0) ) + 1.0 ) * (:rho: - :rhoA:) + :rhoA:
"""
ss.EOM(eom)

#:umag:      =  numpy.sqrt(:u:*:u: + :v:*:v:)

# Initialize variables
x = ss.mesh.coords[0].data
y = ss.mesh.coords[1].data
z = ss.mesh.coords[2].data

dx = x[1, 0, 0] - x[0, 0, 0]

rad = numpy.sqrt((x - numpy.pi) ** 2)  # + (y-numpy.pi)**2 ) #+ (z-numpy.pi)**2  )

diffAmp = 1.0
if diffusive:
    diffAmp = 1.0e6

# def swirl(x,y,t):
#    u = - numpy.sin


if not twoD:
    ss.variables["rho"].data += ss.gfilter(
        numpy.where(x < numpy.pi * 0.25, diffAmp, 1.0)
    )
    ss.variables["rhoA"].data += 1.0e6 * (
        1.0 + 0.1 * numpy.exp(-((x - numpy.pi * 1.5) ** 2) / 0.01)
    )
    ss.variables["phi"].data = 3.14159 - x
    ss.variables["uA"].data += 1.0
    ss.variables["u"].data += 1.0
    ss.variables["u1"].data += 1.0

else:
    r = numpy.sqrt((x - numpy.pi / 2.0 * 3.0) ** 2 + (y - numpy.pi / 2.0 * 3.0) ** 2)
    ss.variables["rho"].data += 1.0  # + 1.0e6*numpy.exp( -numpy.abs(r) / (20.0*dx)**2 )
    ss.variables["rhoA"].data += 1.0e6

    r = numpy.sqrt((x - numpy.pi / 2.0) ** 2 + (y - numpy.pi / 2.0) ** 2)
    ss.variables["phi"].data = ss.gfilter(numpy.minimum(r - numpy.pi / 4.0, 10.0 * dx))
    if False:
        ss.variables["vA"].data += 1.0
        ss.variables["v"].data += 1.0
        ss.variables["v1"].data += 1.0
        ss.variables["uA"].data += 1.0
        ss.variables["u"].data += 1.0
        ss.variables["u1"].data += 1.0


# Try circles


# Init momenta
# ss.variables['u1'].data = 1.0

ss.variables["dx"].data += x[1, 0, 0] - x[0, 0, 0]
ss.variables["dx4"].data += (x[1, 0, 0] - x[0, 0, 0]) ** 4


ss.updateVars()

time = 0.0
viz = True

v = 2.0 * numpy.pi

dt_max = L / ss.mesh.nn[0] / v * 0.75

tt = v * 1.0 / 2.0
if not twoD:
    tt *= 0.2


xx = ss.PyMPI.zbar(x)
yy = ss.PyMPI.zbar(y)

dt = dt_max
cnt = 1

v2 = ss.PyMPI.zbar(ss.variables["rhoT"].data)
# plt.figure(2)
# plt.plot( xx[:,0],v2[:,0] )
# plt.contourf( xx,yy,v2 ,64 )
# plt.pause(.001)
# import pdb
# pdb.set_trace()

viz_freq = 5
if twoD:
    viz_freq = 10


viz_vars = ["rhoT", "phi", "u", "v"]
if viz and (not test):
    ss.write(viz_vars)

while tt > time:
    time = ss.rk4(time, dt)
    dt = min(dt_max, (tt - time))

    # Poor man's boundary condition
    if not twoD:
        ss.variables["rho"].data[0, :, :] = diffAmp

    ss.iprint("%s -- %s" % (cnt, time))
    cnt += 1
    if viz and (not test):
        v1 = ss.PyMPI.zbar(ss.variables["rhoA"].data)
        v2 = ss.PyMPI.zbar(ss.variables["rhoT"].data)
        v = ss.PyMPI.zbar(ss.variables["rho"].data)
        vA = ss.PyMPI.zbar(ss.variables["rhoA"].data)
        if cnt % viz_freq == 0:
            ss.write(viz_vars)

        if ss.PyMPI.master and (cnt % viz_freq == 0) and True:
            plt.figure(1)
            plt.clf()
            plt.plot(xx[:, 0], v2[:, 0], "k.-")
            # plt.plot(xx[:,0],v[:,0] ,'r.-')
            # plt.plot(xx[:,0],vA[:,0] ,'b.-')
            plt.pause(0.001)
            if twoD:
                plt.figure(2)
                plt.clf()
                plt.contourf(xx, yy, v2, 8)
                plt.contour(xx, yy, v2, 1)
                plt.pause(0.001)

            # import pdb
            # pdb.set_trace()


if viz and (not test):
    ss.write(viz_vars)


v2 = ss.PyMPI.zbar(ss.variables["rhoT"].data)
error = numpy.min(v2[Npts / 2 :, 0])
ss.iprint(error)
