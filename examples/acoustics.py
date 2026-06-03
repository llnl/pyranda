from __future__ import print_function
import sys
import time
import numpy
import matplotlib.pyplot as plt
from matplotlib import cm

from scipy import fftpack

from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep, pyrandaIBM, pyrandaProbes


# Try to get args for testing
try:
    Npts = int(sys.argv[1])
except Exception:
    Npts = 64

try:
    test = bool(int(sys.argv[2]))
except Exception:
    test = False

try:
    testName = sys.argv[3]
except Exception:
    testName = None


## Define a mesh
# Npts = 32
L = numpy.pi * 2.0
dim = 2
gamma = 1.4

problem = "scattering_test"

Lp = L * (Npts - 1.0) / Npts

from meshTest import zoomMesh_solve

dxf = 4 * Lp / float(Npts) * 0.5
xS = zoomMesh_solve(Npts, -2.0 * Lp, 2.0 * Lp, -5.0, 5.0, 3.0, dxf)


def zoomMesh(i, j, k):
    x = xS[i]
    y = xS[j]
    z = 0.0
    return x, y, z


def spec_1d(y, han=True):
    if han:
        ave = numpy.mean(y)
        hw = numpy.hanning(y.shape[0])
        y = (y - ave) * hw + ave
    f1 = fftpack.fft(y)
    f2 = f1 / (numpy.prod(f1.shape) / 2.0)
    # psd = npy.abs(f2) #*npy.conjugate(f2))
    psd = (f2) * numpy.conjugate(f2)
    psd = numpy.abs(psd)
    return psd[0 : int(y.shape[0] / 2)]


mesh_options = {}
mesh_options["coordsys"] = 3
mesh_options["function"] = zoomMesh
mesh_options["periodic"] = numpy.array([False, False, True])
mesh_options["dim"] = 3
mesh_options["x1"] = [-2 * Lp, -2 * Lp, 0.0]
mesh_options["xn"] = [2 * Lp, 2 * Lp, Lp]
mesh_options["nn"] = [Npts, 1, 1]
if dim == 2:
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
k0 = 1.0  # Wave number of sinewave


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
# Apply constant BCs
[:u:,:v:,:w:] = ibmV( [:u:,:v:,0.0], :phi:, [:gx:,:gy:,:gz:], [:u1:,:u2:,0.0] )
:rho: = ibmS( :rho: , :phi:, [:gx:,:gy:,:gz:] )
:p:   = ibmS( :p:   , :phi:, [:gx:,:gy:,:gz:] )
:inlet: = where( s0*simtime > (meshx+4.*pi) , p0*.01*sin( (simtime-(meshx+4.*pi)/s0) / 3.0 * 2.0*pi*k0 ), 0.0 )
bc.extrap(['rho','p','u'],['xn'])
#bc.const(['u'],['y1','yn'],u0)
bc.const(['v'],['x1','xn','y1','yn'],0.0)
#bc.const(['rho'],['x1','y1','yn'],rho0)
#bc.const(['p'],['x1','y1','yn'],p0)
bc.field(['p'],  ['x1','xn','y1','yn'],p0 + :inlet:)
bc.field(['rho'],['x1','xn','y1','yn'],rho0 + :inlet:/(s0*s0))
bc.field(['u'],  ['x1','xn','y1','yn'],:inlet:/(s0*rho0))
:Et:  = :p: / ( :gamma: - 1.0 )  + .5*:rho:*(:u:*:u: + :v:*:v:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
:dtB: = 0.2* dt.diff(:beta:,:rho:)
:dt: = numpy.minimum(:dt:,:dtB:)
:umag: = sqrt( :u:*:u: + :v:*:v: )
"""
eom = (
    eom.replace("u0", str(u0))
    .replace("p0", str(p0))
    .replace("rho0", str(rho0))
    .replace("s0", str(s0))
    .replace("k0", str(k0))
)


# Add the EOM to the solver
ss.EOM(eom)


# Initialize variables
ic = """
:gamma: = 1.4
:R: = 1.0
:cp: = :R: / (1.0 - 1.0/:gamma: )
:cv: = :cp: - :R:
#rad = sqrt( (meshx-pi)**2  +  (meshy-pi)**2 ) 3
rad = sqrt( meshx**2  +  meshy**2 ) 
:phi: = rad - pi/4.0
:rho: = 1.0 + 3d()
:p:  =  1.0 + 3d() #exp( -(meshx-1.5)**2/.25**2)*.1
#:u: = where( :phi:>0.5, mach * sqrt( :p: / :rho: * :gamma:) , 0.0 )
#:u: = mach * sqrt( :p: / :rho: * :gamma:)
:u: = 0.0 + 3d() #gbar( gbar( :u: ) )
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
tt = 30.0  #

# Start time loop
cnt = 1
viz_freq = 100
pvar = "p"


CFL = 1.0
dt = ss.variables["dt"].data * CFL * 0.01


# Make a probe set for diagnostics
nProbes = Npts
theta = numpy.linspace(0, 2.0 * numpy.pi, nProbes + 1)[:-1]
R0 = 5.0
x = R0 * numpy.cos(theta)
y = R0 * numpy.sin(theta)
probes = pyrandaProbes(ss, x=x, y=y, z=None)

pressure = []

ss.plot.figure(2)
ss.plot.clf()
ss.plot.contourf(pvar, 64, cmap=cm.jet)
ss.plot.contour("phi", [0.0])
ss.plot.showGrid()


dtSample = 1.0e5
while tt > time:
    # Update the EOM and get next dt
    time = ss.rk4(time, dt)
    dt = min(ss.variables["dt"].data * CFL, dt * 1.1)
    dt = min(dt, (tt - time))

    if (time > 18.0) and (time < 20.0):
        dtSample = min(dt, dtSample)

    if time > 20.0:
        dt = dtSample
        prb = probes.get(pvar)
        pressure.append(prb)

    # Print some output
    ss.iprint("%s -- %s --- %f" % (cnt, time, dt))
    cnt += 1
    if viz and (not test):
        v = ss.PyMPI.zbar(ss.variables[pvar].data)
        phi = ss.PyMPI.zbar(ss.variables["phi"].data)
        if cnt % viz_freq == 1:
            ss.write()

        if cnt % viz_freq == 1:  # or True:
            ss.plot.figure(2)
            ss.plot.clf()
            ss.plot.contourf(pvar, 64, cmap=cm.jet)
            ss.plot.contour("phi", [0.0])


ss.writeRestart()


nTimes = len(pressure)
pdata = numpy.array(pressure)

amp = numpy.zeros(nProbes)
for d in range(pdata.shape[1]):
    data = pdata[:, d]
    p1d = spec_1d(data - numpy.mean(data), han=True)
    amp[d] = numpy.sqrt(numpy.max(p1d))


if ss.PyMPI.master:
    plt.figure(3)
    ax = plt.subplot(111, projection="polar")
    ax.plot(theta, amp)

ss.plot.figure(2)
probes.plot()
