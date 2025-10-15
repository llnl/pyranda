import re
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC
from pyranda.pyranda import pyrandaRestart


## Define a mesh
L = numpy.pi * 2.0
Npts = 200
Lp = L * (Npts-1.0) / Npts

imesh = """
xdom = (0.0, Lp, Npts)
""".replace('Lp',str(Lp)).replace('Npts',str(Npts))

# Initialize a simulation object on a mesh
ss = pyrandaSim('sod',imesh)
ss.addPackage( pyrandaBC(ss) )


def minmod(a, b):
    return numpy.where(numpy.abs(a) < numpy.abs(b), a, b) * (numpy.sign(a) == numpy.sign(b))

def van_leer(r):
    return (r + numpy.abs(r)) / (1 + numpy.abs(r))

def superbee(r):
    return numpy.maximum(0, numpy.maximum(numpy.minimum(2*r, 1), numpy.minimum(r, 2)))

def mc(r):
    return numpy.maximum(0, numpy.minimum((1 + r)/2, numpy.minimum(2, 2*r)))

def barth_jespersen(r):
    return numpy.where(r <= 0, 0, numpy.minimum(1, r))

def compute_limiter(r, limiter_type='minmod'):
    if limiter_type == 'minmod':
        return numpy.maximum(0, numpy.minimum(1, r))
    elif limiter_type == 'vanleer':
        return van_leer(r)
    elif limiter_type == 'superbee':
        return superbee(r)
    elif limiter_type == 'mc':
        return mc(r)
    elif limiter_type == 'barth-jespersen':
        return barth_jespersen(r)
    elif limiter_type == 'none':
        return 1.0    
    elif limiter_type == 'low':
        return 0.0    
    else:
        raise ValueError("Unknown limiter type: {}".format(limiter_type))


def upwind_flux_old(flux, vel):
    """
    Compute upwind flux at cell faces from zone-centered velocities.
    flux: numpy array of cell-centered values, shape (nx, ...)
    vel: numpy array of cell-centered velocities, shape (nx, ...)
    Returns: numpy array of upwinded flux at faces, shape (nx+1, ...)
    """
    nx = flux.shape[0]

    # Interpolate velocity to faces
    vel_face = numpy.zeros((nx+1,) + vel.shape[1:], dtype=vel.dtype)
    vel_face[1:-1,:,:] = 0.5 * (vel[:-1,:,:] + vel[1:,:,:])
    # Optionally, set boundary faces to zero or extrapolate:
    vel_face[0,:,:] = vel[0,:,:]
    vel_face[-1,:,:] = vel[-1,:,:]

    # Upwind flux at faces
    up_flux = numpy.empty((nx+1,) + flux.shape[1:], dtype=flux.dtype)
    # For each face, choose upwind value
    left = flux[:-1,:,:]
    right = flux[1:,:,:]
    # For interior faces
    up_flux[1:-1,:,:] = numpy.where(vel_face[1:-1] > 0, left, right)

    
    # For boundary faces, use boundary cell
    up_flux[0,:,:] = flux[0,:,:]
    up_flux[-1,:,:] = flux[-1,:,:]

    return up_flux

def upwind_flux(flux, vel, U):
    """
    Compute Lax-Friedrichs (Rusanov) numerical flux at cell faces from
    cell-centered state U, its flux f(U), and advection speed vel.

    Parameters
    ----------
    flux : numpy.ndarray
        Physical flux at cell centers, f(U), shape (nx, ...).
    vel : numpy.ndarray
        Cell-centered advection speed a(x), shape (nx, ...).
    U : numpy.ndarray
        Cell-centered state, shape (nx, ...).

    Returns
    -------
    F_face : numpy.ndarray
        Face-centered numerical flux, shape (nx+1, ...), computed as
          F_{i+1/2} = 0.5*(f_L + f_R) - 0.5*alpha_{i+1/2}*(U_R - U_L),
        with alpha_{i+1/2} = max(|a_i|, |a_{i+1}|).
    """
    import numpy as numpy

    nx = flux.shape[0]
    F_face = numpy.empty((nx + 1,) + flux.shape[1:], dtype=flux.dtype)

    # Interior left/right states and fluxes
    UL = U[:-1, ...]
    UR = U[1:, ...]
    fL = flux[:-1, ...]
    fR = flux[1:, ...]

    # Per-face maximum wave speed alpha = max(|a_L|, |a_R|)
    alpha_int = numpy.maximum(numpy.abs(vel[:-1, ...]), numpy.abs(vel[1:, ...]))

    # Lax-Friedrichs flux on interior faces
    F_face[1:-1, ...] = 0.5 * (fL + fR) - 0.5 * alpha_int * (UR - UL)

    # Boundary faces: one-sided flux (use adjacent cell)
    F_face[0, ...] = flux[0, ...]
    F_face[-1, ...] = flux[-1, ...]

    return F_face



def myDiv(pysim, flux, vel, U, limiter_type='barth-jespersen'):

    # Low-order flux (upwind) at faces
    loflux = upwind_flux_old(flux, vel )   # shape (nx+1, ...)
    loflux = upwind_flux(flux, vel, U )   # shape (nx+1, ...)
    

    
    # High-order flux at faces
    hoflux = numpy.zeros_like(loflux)
    hoflux[:-1,:,:] = pysim.interp_z2fx(flux)  # shape (nx+1, ...)

    #plt.figure(11)
    #plt.clf()
    #plt.plot(loflux[:,0,0],'k-')
    #plt.plot(hoflux[:,0,0],'b-')
    #plt.pause(.1)

    
    # Compute differences for limiter at interior faces
    dF_forward = hoflux[1:-1,:,:] - hoflux[0:-2,:,:]
    dF_backward = hoflux[2:,:,:] - hoflux[1:-1,:,:]
    epsilon = 1e-12
    r = dF_forward / (dF_backward + epsilon)

    phi = compute_limiter(r, limiter_type=limiter_type)

    #import pdb
    #pdb.set_trace()
    
    sensor = pysim.var('sensor').data
    phi = 1.0 - (sensor[1:,:,:] + sensor[0:-1,:,:] ) / 2.0
    

    # Blend high- and low-order fluxes at interior faces
    flux_limited = numpy.empty_like(hoflux)
    flux_limited[0,:,:] = loflux[0,:,:]
    flux_limited[-1,:,:] = loflux[-1,:,:]
    flux_limited[1:-1,:,:] = loflux[1:-1,:,:] + phi * (hoflux[1:-1,:,:] - loflux[1:-1,:,:])

    #import pdb
    #pdb.set_trace()

    
    # Compute divergence at cell centers
    div = numpy.empty_like(flux)
    div[:,:,:] = (flux_limited[1:,:,:] - flux_limited[:-1,:,:]) / pysim.dx
    #div[1:-1,:,:] = (flux_limited[2:,:,:] - flux_limited[1:-1,:,:]) / pysim.dx
    #div[0,:,:] = 0.0
    #div[-1,:,:] = 0.0

    return div

def myDivOld( pysim, flux ):

    hoflux = pysim.interp_z2fx( flux )

    # Make a low-order flux here?

    
    div = numpy.empty_like( flux )
    div[1:-1,:,:] = ( hoflux[2:,:,:] - hoflux[1:-1,:,:] ) / pysim.dx

    # BCs
    div[0 ,:,:] = 0.0
    div[-1,:,:] = 0.0

    return div




ss.addUserDefinedFunction('fvd',myDiv)


# Define the equations of motion
eom ="""
# Primary Equations of motion here
:c: = sqrt( :p:*:gamma: / :rho: )
:vel: = abs(:u:) + :c:
ddt(:rho:)  =  -fvd(:rho:*:u:, :vel:, :rho:)
ddt(:rhou:) =  -fvd(:rhou:*:u: + :p: - :tau: , :vel:, :rhou: )
ddt(:Et:)   =  -fvd( (:Et: + :p: - :tau:)*:u: , :vel:, :Et: )
# Conservative filter of the EoM
#:rho:       =  fbar( :rho:  )
#:rhou:      =  fbar( :rhou: )
#:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u:) ) * ( :gamma: - 1.0 )
# Artificial bulk viscosity (old school way)
:div:       =  ddx(:u:) 
#:div:       =  fvd(:u:,:u:) 
:beta:      =  gbar( ring(:div:) * :rho:) * 7.0e-2
:tau:       =  :beta:*:div:
:sensor:    = gbar( where( abs( :tau: / :p: ) > 0.005, 1.0, 0.0 ) )
:sensor1:    = gbar( where( abs( :tau: / :p: ) > 0.001, 1.0, 0.0 ) )
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
dt_max = v / ss.mesh.nn[0] * 0.25
tt = L/v * .25 #dt_max

# Start time loop
dt = dt_max
cnt = 1
viz_freq = 50
pvar = 'u'
viz = True

while tt > time:

    # Update the EOM and get next dt
    time = ss.rk4(time,dt)
    dt = min(dt_max, (tt - time) )
    
    # Print some output
    ss.iprint("%s -- %s" % (cnt,time)  )
    cnt += 1
    if viz:

        if (cnt%viz_freq == 0):
            ss.plot.figure(1)
            plt.clf()
            ss.plot.plot('p','b.-')
            ss.plot.plot('u','r.-')
            ss.plot.plot('rho','g.-')
            ss.plot.plot('sensor','k.-')
            input("paused")

        
ss.writeGrid()
ss.write()



