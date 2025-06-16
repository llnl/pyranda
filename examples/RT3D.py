import re
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep


# Try to get args for testing
try:
    Npts = int(sys.argv[1])
except:
    Npts = 64

try:
    is2D = bool(int(sys.argv[2]))
except:
    is2D = False

try:
    test = bool(int(sys.argv[3]))
except:
    test = False


## Define a mesh
if is2D:
    problem = 'RT_2D'
    imesh = """
    xdom = (0.0, 2.8*pi , int(Npts*1.4), periodic=False) 
    ydom = (0.0, 2*pi*FF,  Npts, periodic=True)
    zdom = (0.0, 2*pi*FF,  1, periodic=True)
    """.replace('Npts',str(Npts)).replace('pi',str(numpy.pi)).replace('FF',str( float(Npts-1)/Npts ) )
    waveLength = 4
else:
    problem = 'RT_3D'
    imesh = """
    xdom = (0.0, 3.0*pi , int(Npts*1.5), periodic=False) 
    ydom = (0.0, 2*pi*FF,  Npts, periodic=True)
    zdom = (0.0, 2*pi*FF,  Npts, periodic=True)
    """.replace('Npts',str(Npts)).replace('pi',str(numpy.pi)).replace('FF',str( float(Npts-1)/Npts ) )
    waveLength = 4
    
rho_l = 1.0
rho_h = 3.0
mwH  = 3.0
mwL  = 1.0
gx    = -0.01
Runiv = 1.0
CPh   = 1.4
CVh   = 1.0
CPl   = 1.4
CVl   = 1.0


parm_dict = {'gx':gx,
             'CPh':CPh,'CPl':CPl,
             'CVh':CVh,'CVl':CVl,
             'mwH':mwH,'mwL':mwL,
             'Runiv':Runiv,'waveLength':waveLength,
             'rho_l':rho_l,'rho_h':rho_h,
             'delta':2.0* numpy.pi / Npts * 4.0 }

# Initialize a simulation object on a mesh
ss = pyrandaSim(problem,imesh)
ss.addPackage( pyrandaTimestep(ss) )
ss.addPackage( pyrandaBC(ss) )


# User defined function to take mean of data and return it as 3d field, dataBar
def meanXto3d(pysim,data):
    meanX = pysim.PyMPI.yzsum( data ) / (pysim.PyMPI.ny * pysim.PyMPI.nz)
    tmp = pysim.emptyScalar()
    for i in range( tmp.shape[0]):
        ii = int( pysim.mesh.indices[0].data[i,0,0] )
        tmp[i,:,:]  =  meanX[ii]
    return tmp

ss.addUserDefinedFunction("xbar",meanXto3d)


# Define the equations of motion
eom ="""
# Primary Equations of motion here
ddt(:rhoYh:)  =  -ddx(:rhoYh:*:u: - :Jx:)    - ddy(:rhoYh:*:v: - :Jy:)   - ddz(:rhoYh:*:w: - :Jz:)
ddt(:rhoYl:)  =  -ddx(:rhoYl:*:u: + :Jx:)    - ddy(:rhoYl:*:v: + :Jy:)   - ddz(:rhoYl:*:w: + :Jz:)
ddt(:rhou:)   =  -ddx(:rhou:*:u: - :tauxx:)  - ddy(:rhou:*:v: - :tauxy:) - ddz(:rhou:*:w: - :tauxz:) + :rho:*gx
ddt(:rhov:)   =  -ddx(:rhov:*:u: - :tauxy:)  - ddy(:rhov:*:v: - :tauyy:) - ddz(:rhov:*:w: - :tauyz:)
ddt(:rhow:)   =  -ddx(:rhow:*:u: - :tauxz:)  - ddy(:rhow:*:v: - :tauyz:) - ddz(:rhow:*:w: - :tauzz:)
ddt(:Et:)     =  -ddx( (:Et: - :tauxx:)*:u: - :tauxy:*:v: - :tauxz:*:w: - :tx:*:kappa:) - ddy( (:Et: - :tauyy:)*:v: -:tauxy:*:u: - :tauyz:*:w: - :ty:*:kappa:) - ddz( (:Et: - :tauzz:)*:w: - :tauxz:*:u: - :tauyz:*:v: - :tz:*:kappa:) + :rho:*gx*:u:
# Conservative filter of the EoM
:rhoYh:     =  fbar( :rhoYh:  )
:rhoYl:     =  fbar( :rhoYl:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:rhow:      =  fbar( :rhow: )
:Et:        =  fbar( :Et:   )
# Xbar operator
:mybar: = xbar(:rho:*:u:*:u:)
# Update the primatives and enforce the EOS
:rho:       = :rhoYh: + :rhoYl:
:Yh:        =  :rhoYh: / :rho:
:Yl:        =  :rhoYl: / :rho:
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:w:         =  :rhow: / :rho:
:cv:        = :Yh:*CVh + :Yl:*CVl
:cp:        = :Yh:*CPh + :Yl:*CPl
:gamma:     = :cp:/:cv:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v:) ) * ( :gamma: - 1.0 )
:mw:        = 1.0 / ( :Yh: / mwH + :Yl: / mwL )
:R:         = Runiv / :mw:
:T:         = :p: / (:rho: * :R: )
# Artificial bulk viscosity / shear viscosity
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
:Yx:        =  ddx(:Yh:)
:Yy:        =  ddy(:Yh:)
:Yz:        =  ddz(:Yh:)
:enst:      = sqrt( (:uy:-:vx:)**2 + (:uz: - :wx:)**2 + (:vz:-:wy:)**2 )
:S:         = sqrt( :ux:*:ux: + :vy:*:vy: + :wz:*:wz: + .5*((:uy:+:vx:)**2 + (:uz: + :wx:)**2 + (:vz:+:wy:)**2) )
:mu:        = 1.0e-4 * gbar( abs(ring(:S:  )) ) * :rho:
:beta:      = 7.0e-2 * gbar( abs(ring(:div:)) * :rho: )
# Artificial species diffusivities
:Dsgs:      =  1.0e-4 * ring(:Yh:)
:Ysgs:      =  1.0e2  * (abs(:Yh:) - 1.0 + abs(1.0-:Yh: ) )*gridLen**2
:adiff:     =  gbar( :rho:*numpy.maximum(:Dsgs:,:Ysgs:) / :dt: )
:Jx:        =  :adiff:*:Yx:
:Jy:        =  :adiff:*:Yy:
:Jz:        =  :adiff:*:Yz:
:taudia:    =  (:beta:-2./3.*:mu:) *:div: - :p:
:tauxx:     =  2.0*:mu:*:ux:   + :taudia:
:tauyy:     =  2.0*:mu:*:vy:   + :taudia:
:tauzz:     =  2.0*:mu:*:wz:   + :taudia:
:tauxy:     = :mu:*(:uy:+:vx:) 
:tauxz:     = :mu:*(:uz:+:wx:) 
:tauyz:     = :mu:*(:vz:+:wz:) 
[:tx:,:ty:,:tz:] = grad(:T:)
:kappa:     = 1.0e-3 * gbar( ring(:T:)* :rho:*:cv:/(:T: * :dt: ) )
:cs:        = sqrt( :p: / :rho: * :gamma: )
# Time step control routines
:dt:        = dt.courant(:u:,:v:,:w:,:cs:)*1.0
:dt:        = numpy.minimum(:dt:,0.2 * dt.diff(:beta:,:rho:))
:dt:        = numpy.minimum(:dt:,0.2 * dt.diff(:mu:,:rho:))
:dt:        = numpy.minimum(:dt:,0.2 * dt.diff(:adiff:,:rho:))
# Add some BCs in the x-direction
bc.const(['Yh'],['xn'],1.0)
bc.const(['Yh'],['x1'],0.0)
bc.const(['Yl'],['x1'],1.0)
bc.const(['Yl'],['xn'],0.0)
# bc.extrap(['rho','Et'],['x1','xn'])
bc.const(['u','v','w'],['x1','xn'],0.0)
# Sponge BCs
#:leftBC:  = 1.0 - 0.5 * (1.0 + tanh( (meshx-1.0) / .1 ) )
#:rightBC: = 0.5 * (1.0 + tanh( (meshx-8.0) / .1 ) )
#:BC: = numpy.maximum(:leftBC:,:rightBC:)
#:u: = gbar(:u:)*:BC: + :u:*(1.0 - :BC:)
#:rho: = gbar(:rho:)*:BC: + :rho:*(1.0 - :BC:)
#:Et: = gbar(:Et:)*:BC: + :Et:*(1.0 - :BC:)
:rhoYh: = :rho:*:Yh:
:rhoYl: = :rho:*:Yl:
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:rhow: = :rho:*:w:
:mix:  = 4.0*:Yh:*:Yl:
"""

# Add the EOM to the solver
ss.EOM(eom,parm_dict)


# Initialize variables
ic = """
:gamma:= 5./3.
p0     = 1.0
# At     = ( (rho_h - rho_l)/(rho_h + rho_l) )
u0     = 0.0 # sqrt( abs(gx*At/waveLength) ) * 0.1 # velocity oscillations deprecated MDW in favor of perturbed surface
T0     = 1.0
offset = 1.4*pi + 0.5 * sin(2 * meshy) # surface perturbed with sine wave MDW

# Add smooth/diffuse interface
:Yh: = 0.5 * (1.0 + tanh( (meshx-offset) / delta ))
:Yl: = 1.0 - :Yh:

# Define R mixture
:mw:        = 1.0 / ( :Yh: / mwH + :Yl: / mwL )
:R:         = Runiv / :mw:
R_l = Runiv / mwL
R_h = Runiv / mwH

# Define analytical solution for pressure for hydrostatic equilibrium (see RTI_dPdZ_init file for derivation)
A = gx * (R_h + R_l) / (2.0 * Runiv * T0 * R_l * R_h)
B = gx * delta * (R_l - R_h) / (2.0 * Runiv * T0 * R_l * R_h)
:p: = p0 * exp(A * (meshx-offset)) * (numpy.cosh( (meshx-offset) / delta ))**(B)

# If no diffuse interface, Yl,Yh are just constants
# :p: = p0 * exp(( gx * (R_h * :Yl: + R_l * :Yh:) / (2.0 * Runiv * T0 * R_l * R_h)) * (meshx-offset))

# # Pressure initial pressure distribution
# def P_init(meshx,gx,rho_h,rho_l,offset,p0):
#     if ( meshx <= offset):
#         P = p0 + rho_l * gx * (meshx-offset)
#     if ( meshx > offset):
#         P = p0 - rho_h * gx * (meshx-offset)
#     return P
# Define ideal profile
# :p: = P_init(meshx,gx,rho_h,rho_l,offset,p0)

# Add oscillations near interface
wgt = 4*:Yh:*:Yl:
:v: = 3d(0.0)
:w: = 3d(0.0)
:u: = 3d(0.0) # wgt * gbar( (random3D()-0.5)*u0 ) # velocity oscillations deprecated MDW in favor of perturbed surface

# Primitive variables
:T:         = 3d(T0)
# :T:         = :p: / (:rho: * :R: )
:rho:       = :p: / ( :R: * :T: ) # rho_h * :Yh: + rho_l * :Yl: # deprecated MDW
:error:     = ddx(:p:) - :rho:*gx
:cv:        = :Yh:*CVh + :Yl:*CVl
:cp:        = :Yh:*CPh + :Yl:*CPl
:gamma:     = :cp:/:cv:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)

# Form conserved quatities
:rhoYh: = :rho:*:Yh:
:rhoYl: = :rho:*:Yl:
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:rhow: = :rho:*:w:
:Et:  = :p: / (:gamma:-1.0) + 0.5*:rho:*(:u:*:u: + :v:*:v: + :w:*:w:)
"""

# Set the initial conditions
numpy.random.seed(1234)
ss.setIC(ic,parm_dict)
    

# Write a time loop
time = 0.0

# Start time loop
CFL = 1.0
dt = ss.variables['dt'].data * CFL

# Viz/IO
viz_freq = 20
dmp_freq = 200

tstop = 100.0
dtmax = dt * .1

outVars = ['p','u','v','w','rho','Yh']
ss.write(outVars)


mixW  = []
timeW = []

while time < tstop:

    # Update the EOM and get next dt
    time = ss.rk4(time,dt)
    dt = ss.variables['dt'].data * CFL
    dt = min( dt, dtmax*1.1)
    dtmax = dt*1.0

    umax = ss.var('u').mean()
    ss.iprint("%s -- %s --- Umax: %s " % (ss.cycle,time,umax)  ) 

    if ( ss.cycle%viz_freq == 0):
        mixW.append( numpy.trapz( ss.var('mix').mean( axis=[1,2] ) , ss.var('meshx')[:,0,0] )   )
        timeW.append( time )
    
    if not test:
        if (ss.cycle%viz_freq == 0 ) : #True) : # replace with True to plot each timestep

            # 2D contour plots
            ss.plot.figure(1)
            ss.plot.clf()
            ss.plot.contourf( 'rho', 32, cmap='jet')

            ss.plot.figure(2)
            ss.plot.clf() # comment these lines to have each timestep plot overlayed
            ss.plot.plot('error')

            ss.plot.figure(3)
            ss.plot.clf()
            ss.plot.plot('p')

            ss.plot.figure(4)
            ss.plot.clf()
            ss.plot.plot('rho')
            # input('stop here') # uncomment to stop after each cyle

            
        if ( ss.cycle%dmp_freq == 0) :
            ss.write(outVars)


if not test:

    if ss.PyMPI.master:
        plt.figure(3)
        plt.plot(timeW,mixW)
        plt.pause(.1)
        
else:
    fname = problem + '.dat'
    numpy.savetxt( fname  , (timeW,mixW) )
    print(fname)
