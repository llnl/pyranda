################################################################################
# Copyright (c) 2018, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
#
# Written by: Britton J. Olson, olson45@llnl.gov
#
#
# Description: Problem script for GAtech divergent shock tube experiment
#      of Musci and Ranjan.
#
# Notes:
#
#
#
################################################################################
import numpy
import matplotlib.pyplot as plt
from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep
from pyranda import pyrandaProbes, pyrandaRestart
from optparse import OptionParser
################################################################################


##### Command line parser  #####################################################
parser = OptionParser()
parser.add_option(
    "-r",
    "--resolution",
    dest="res",
    default=1,
    type="int",
    help="Resolution multiplier of the mesh",
)
parser.add_option(
    "-d",
    "--dump",
    action="store_true",
    dest="restart",
    default=False,
    help="Run in min",
)
parser.add_option(
    "-k",
    "--kernel",
    dest="problem",
    default="GAtechProblem",
    help="Name of the problem name (kernel) to run",
)
parser.add_option(
    "-s",
    "--stop_time",
    dest="stop_time",
    default=0.006,
    type="float",
    help="Give the stop time (seconds) of the simulation",
)
parser.add_option(
    "-b",
    "--blast",
    dest="blastMag",
    default=60.0,
    type="float",
    help="Magnitude of the blast wave over-pressure at t0",
)


(options, args) = parser.parse_args()


restart = options.restart  # Restart this run
problem = options.problem  # Problem name: used for data output
stop_time = options.stop_time  # Stopping time of simulation (seconds)
################################################################################


# Some useful constants
inToM = 0.0254  # Inches to meters converstion
psiToPa = 6894.76  # PSI to Pascal converstion
Runiv = 8.31446261815324  # Universal gas constant ( J K-1 mol-1 )

## Define a mesh
nx = 50 * options.res  # Points in x
ny = 100 * options.res  # Points in y
dim = 2  # Dimension of the problem
open_angle = 45.0  # Total opening angle (degrees)
L_min = 2.62 * inToM  # Width of tube at inner radius (meters)
h = 40.0 * inToM  # Height is y direction (meters)
intH = 21.0 * inToM  # Height of interface (meters)

## Gas properties
p0 = 1.013e5  # Pressure (Pa)
T0 = 300.0  # Kelvin
mwL = 14.0067 / 1e3  # Molecular weight of light
mwH = 44.01 / 1e3  # Molecular weight of heavy
myGamma = 1.4  # Gamma of gases (single gamma)

## Interface specification; simple sine wave
intAmp = 0.004 * h  # amplitude of the modes in the gas (meters)
intFreq = 10.0  # N-modes at gas interface

## Detonator specification
Pamp = options.blastMag  # Pressure spike multiplier of detonator
detRad = L_min / 4.0  # Radius of the detonator (meters)

## Foam BCs specifications (likely needs tuning to experiment)
off = 0  # Center point of tanh relative to walls (units are dx)
width = 3  # Width of the sponge at the boundaries (units are dx)

## Set up time stepping loop parameters
time = 0.0  # Initial time (seconds)
viz = True  # Will there be real-time vizualization
CFL = 1.0  # CFL number
Nviz = 10  # Number of visualization files


## Compute max length (x-dir) based on angle/height
open_angle *= numpy.pi / 180.0  # Convert to radians
L_max = 2.0 * (L_min / 2.0 + h * numpy.tan(open_angle / 2.0))  # Max width of tube

## Mesh function:  This function specifies a mesh of the divergent
##    shock tube  y(j) ~ a(j+jo)^2 + b, x(i) ~ {uniform}
jo = (ny - 1) * (L_min / L_max)
alpha = h / ((ny - 1 + jo) ** 2 - jo**2)
beta = -alpha * jo**2


def ST_mesh(i, j, k):
    y = alpha * (j + jo) ** 2 + beta
    wgt = y / h
    Lh = L_min * (1.0 - wgt) + L_max * wgt
    x = -Lh / 2.0 + Lh * float(i) / (nx - 1)
    z = 0
    return x, y, z


# Define a python dictionary for the mesh/problem
mesh_options = {}
mesh_options["coordsys"] = 3  # This assumes curvilinear grids
mesh_options["function"] = ST_mesh  # Function takes (i,j,k) args and returns [x,y,z]
mesh_options["periodicGrid"] = False  # Is the grid itself periodic?
mesh_options["dim"] = dim  # Dimension of the problem
mesh_options["nn"] = [nx, ny, 1]  # Grid spacing
mesh_options["periodic"] = [False, False, False]


# Define the equations of motion
eom = """
# Primary Equations of motion here
ddt(:rhoYh:) =  -div(:rhoYh:*:u: - :Jx: , :rhoYh:*:v: - :Jy:  )
ddt(:rhoYl:) =  -div(:rhoYl:*:u: + :Jx: , :rhoYl:*:v: + :Jy:  )
ddt(:rhou:)  =  -div(:rhou:*:u: - :tauxx:, :rhou:*:v: - :tauxy:)
ddt(:rhov:)  =  -div(:rhov:*:u: - :tauxy:, :rhov:*:v: - :tauyy:)
ddt(:Et:)    =  -div( (:Et: - :tauxx:)*:u: - :tauxy:*:v: - :tx:*:kappa:, (:Et: - :tauyy:)*:v: - :tauxy:*:u: - :ty:*:kappa: )
# Conservative filter of the EoM
:rhoYh:     =  fbar( :rhoYh:  )
:rhoYl:     =  fbar( :rhoYl:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:rho:       = :rhoYh: + :rhoYl:
:Yh:        =  :rhoYh: / :rho:
:Yl:        =  :rhoYl: / :rho:
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v:) ) * ( :gamma: - 1.0 )
:mw: = 1.0 / ( :Yh: / mwH + :Yl: / mwL )
:R: = Runiv / :mw:
:cp: = :R: / (1.0 - 1.0/:gamma: )
:cv: = :cp: - :R:
:T:         = :p: / (:rho: * :R: )
# Artificial bulk viscosity (old school way)
:div:       =  div(:u:,:v:)
[:ux:,:uy:,:tz:] = grad( :u: )
[:vx:,:vy:,:tz:] = grad( :v: )
:S:         = sqrt( :ux:*:ux: + :vy:*:vy: + .5*((:uy:+:vx:)**2))
:beta:      =  gbar( ring(:div:) * :rho: ) * 7.0e-2
:mu:        =  gbar( abs(ring(:S:  )) ) * :rho: * 1.0e-3
:taudia:    =  (:beta:-2./3.*:mu:) *:div: - :p:
:tauxx:     =  2.0*:mu:*:ux:   + :taudia:
:tauyy:     =  2.0*:mu:*:vy:   + :taudia:
:tauxy:     =  :mu:*(:uy:+:vx:)
[:tx:,:ty:,:tz:] = grad(:T:)
:kappa:     = gbar( ring(:T:)* :rho:*:cv:/(:T: * :dt: ) ) * 1.0e-3
# Artificial species diffusivities
:Dsgs:      =  ring(:Yh:) * 1.0e-4
:Ysgs:      =  1.0e2*(abs(:Yh:) - 1.0 + abs(1.0-:Yh: ) )*gridLen**2
:adiff:     =  gbar( :rho:*numpy.maximum(:Dsgs:,:Ysgs:) / :dt: )
[:Yx:,:Yy:,:Yz:] = grad( :Yh: )
:Jx:        =  :adiff:*:Yx:
:Jy:        =  :adiff:*:Yy:
# Apply boundary conditions
bc.extrap(['rho','p'],['x1','xn','y1'] ,order=1)
bc.extrap(['Yh'],['x1','xn','y1','yn'], order=1 )
bc.slip( [ ['u','v'] ] , ['x1','xn','y1'] )
# Sponge boundaries/walls
:wgt: =         ( 1.0 + tanh( -(meshi-off)/width ) ) * 0.5
:wgt: = :wgt: + ( 1.0 + tanh( (meshi-((NX-1)-off))/width ) ) * 0.5
:u: = :u:*(1-:wgt:) + gbarx(:u:)*:wgt:
:v: = :v:*(1-:wgt:) + gbarx(:v:)*:wgt:
:rho: = :rho:*(1-:wgt:) + gbarx(:rho:)*:wgt:
:p:   = :p:  *(1-:wgt:) + gbarx(:p:)  *:wgt:
# Sponge outflow
:wgt: = ( 1.0 + tanh( (meshy-myH*(1.0-0.025))/ (.025*myH) ) ) * 0.5
:u: = :u:*(1-:wgt:) + gbary(:u:)*:wgt:
:v: = :v:*(1-:wgt:) + gbary(:v:)*:wgt:
:rho: = :rho:*(1-:wgt:) + gbary(:rho:)*:wgt:
:p:   = :p:  *(1-:wgt:) + gbary(:p:)  *:wgt:
bc.extrap(['u','v','p','rho'],['yn'])
# Update the conserved variables
:Et:  = :p: / ( :gamma: - 1.0 )  + .5*:rho:*(:u:*:u: + :v:*:v:)
:Yl:  = 1.0 - :Yh:
:rhoYh: = :rho:*:Yh:
:rhoYl: = :rho:*:Yl:
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
# Compute some max time steps
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dtC: = dt.courant(:u:,:v:,:w:,:cs:)
:dtB: = 0.1* dt.diff(:beta:,:rho:)
:dt:  = numpy.minimum(:dtC:,:dtB:)
:dtM: = 0.2 * dt.diff(:mu:,:rho:)
:dt: = numpy.minimum(:dt:,:dtM:)
:dtY: = 0.2 * dt.diff(:adiff:,:rho:)
:dt: = numpy.minimum(:dt:,:dtY:)
:umag: = sqrt( :u:*:u: + :v:*:v: )
:op: = ( :p: - p0 ) / toPSI
"""

## Create a dictionary for syntax parsing/replacement for the EOM
eomDict = {}
eomDict["mwH"] = mwH
eomDict["mwL"] = mwL
eomDict["myH"] = h
eomDict["p0"] = p0
eomDict["Runiv"] = Runiv
eomDict["NX"] = nx
eomDict["off"] = off
eomDict["width"] = width
eomDict["p0"] = p0
eomDict["toPSI"] = psiToPa

# Initial conditions of the flow
ic = """
# Interaface perturbations and smootings
sinW = height + AMP*sin( (meshx+intWidth/2.0) / intWidth * 2.0 * pi * FREQ )
:Yh: = 0.5 * (1.0 - tanh( (meshy - sinW ) / thick ) )
:Yl: = 1.0 - :Yh:
# Mixed EOS (partial pressures method)
:gamma: = myGamma
:mw: = 1.0 / ( :Yh: / mwH + :Yl: / mwL )
:R: = Runiv / :mw:
:cp: = :R: / (1.0 - 1.0/:gamma: )
:cv: = :cp: - :R:
# Energy pill for detonator model
rad = sqrt( (meshx-0.0)**2  +  (meshy-0.0)**2 ) 
:p:  =  p0 + (Pamp-1.)*p0* exp( - rad**2 / (expRad)**2 )
:rho: = p0 / (:R: * T0 )
:T: = :p: / (:R: * :rho: )
:u: = 3d()
:v: = 3d()
# Form conservatives from primatives
:rhoYh: = :rho:*:Yh:
:rhoYl: = :rho:*:Yl:
:Et: = :p:/( :gamma: - 1.0 ) + .5*:rho:*(:u:*:u: + :v:*:v:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
# Sound speed and initial time step size
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
"""

## Create a dictionary for syntax parsing/replacement for the ICs
icDict = {}
icDict["myGamma"] = myGamma
icDict["mwH"] = mwH
icDict["mwL"] = mwL
intR = intH / h
icDict["height"] = h * intR
icDict["intWidth"] = L_min * intR + L_max * (1.0 - intR)
icDict["thick"] = h / ny * 2
icDict["AMP"] = intAmp
icDict["FREQ"] = intFreq
icDict["Runiv"] = Runiv
icDict["expRad"] = detRad
icDict["p0"] = p0
icDict["T0"] = T0
icDict["Pamp"] = Pamp


## SETUP A PYRANDA SIMULATION OBJECT: From initial conditions OR restart
if not restart:
    # Initialize a simulation object on a mesh
    ss = pyrandaSim(problem, mesh_options)

    # Add packages to simulation object, EOM, and ICs
    ss.addPackage(pyrandaBC(ss))  # BC package allows for "bc.*" functions
    ss.addPackage(pyrandaTimestep(ss))  # Timestep package allows for "dt.*" functions
    ss.EOM(eom, eomDict)  # Add equations of motion
    ss.setIC(ic, icDict)  # Set the initial conditions
    dt = ss.variables["dt"].data * CFL * 0.01  # Setup the initial 'dt'

    over_pressure = []
    ptime = []


else:
    # Initialize a simulation form a restart file
    [ss, local_vars] = pyrandaRestart(problem)  # Restore simulation object
    locals().update(local_vars)  # Restore local vars saved
    time = ss.time  # Resume time
    dt = ss.deltat  # Last deltat
    import pdb

    pdb.set_trace()


## Make a probe set for diagnostics
xprb = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # X-positions of the probes
yprb = [0.1, 0.2, 0.4, 0.55, 0.7, 0.8]  # Y-positions of the probes
probes = pyrandaProbes(ss, x=xprb, y=yprb, z=None)  # Create a probe set


## Variables to write to visualization
wvars = ["Yh", "rho", "u", "v", "p", "beta", "kappa", "adiff", "mu", "T", "op"]
viz_freq = stop_time / float(Nviz)
viz_dump = time + viz_freq

## Frequency of restarts
dump_freq = 3000

## Frequency of probe queries
probe_freq = 10

## Plot the solution of 'rho'
if viz:
    ss.plot.figure(1)
    ss.plot.contourf("rho", 64)

## Write viz-files
ss.write(wvars)


## Time loop to advance the solution
while time < stop_time:
    # Update the EOM and get next dt
    time = ss.rk4(time, dt)
    dt = min(ss.variables["dt"].data * CFL, dt * 1.01)
    dt = min(dt, (stop_time - time))

    # Time step stability flags; what's limiting the dt?
    stab_type = ""
    if dt == ss.variables["dtB"].data:
        stab_type = "bulk"  # Shock limited
    elif dt == ss.variables["dtY"].data:
        stab_type = "diffusion"  # Species interface limited
    elif dt == ss.variables["dtM"].data:
        stab_type = "shear"  # Vorticity limited
    elif dt == ss.variables["dtC"].data:
        stab_type = "CFL"  # Courant-Friedrichs-Lewy (CFL)
    else:
        stab_type = "ramp"  # Slow rampup

    # Simulation heart-beat
    ss.iprint(
        "Cycle: %5d --- Time: %10.4e --- deltat: %10.4e ( %s )"
        % (ss.cycle, time, dt, stab_type)
    )

    # Sample probes every 'probe_freq' cycles
    if ss.cycle % probe_freq == 0:
        prb = probes.get("op")  # Compute probes' values for variable 'op'
        over_pressure.append(prb)  # Add these probes to a running list
        ptime.append(time)  # Keep track of time

    # Constant time
    if time > viz_dump:
        ss.iprint("Writing viz file........")
        ss.write(wvars)
        viz_dump += viz_freq

        # Native vis (matplotlib)
        if viz:
            # Plot 'p' contours on figure 1
            ss.plot.figure(1)
            ss.plot.clf()
            ss.plot.contourf("p", 64)

            # Plot pressure histories of the probes on figure 2
            if ss.PyMPI.master:
                plt.figure(2)
                op = numpy.array(
                    over_pressure
                )  # Convert to numpy array for plotting/slicing
                for ii in range(op.shape[1]):
                    plt.plot(ptime, op[:, ii])
                plt.pause(0.01)

    # What local data should persist?
    my_local_data = {"over_pressure": over_pressure, "ptime": ptime}

    # Write a full restart file every 'dump_freq' time steps
    if ss.cycle % dump_freq == 0:
        ss.writeRestart(ivars=my_local_data)


# Final IO
ss.write(wvars)  # Viz files for VISIT
ss.writeRestart(ivars=my_local_data)  # Restart data


# Write probes data to file
if ss.PyMPI.master:
    op = numpy.array(over_pressure)
    tt = numpy.array(ptime).reshape(len(ptime), 1)
    data = numpy.append(tt, op, axis=1)
    numpy.savetxt(problem + "/probes.txt", data)

# Simple scalar for max-val of first probe
op = numpy.array(over_pressure)
imax = numpy.argmax(op[:, 0])
val = op[imax, 0]
vtime = ptime[imax]
ss.iprint("Max val:  %s" % val)
ss.iprint("Max time: %s" % vtime)
