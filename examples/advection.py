import sys
import time
import numpy 
import matplotlib.pyplot as plt

from pyranda import pyrandaSim

# Try to get args
try:
    Npts = int(sys.argv[1])
except:
    Npts = 100

try:
    test = bool(int(sys.argv[2]))
except:
    test = False

## Define a mesh
L = numpy.pi * 2.0  
Lp = L * (Npts-1.0) / Npts


# Define the domain/mesh
imesh = """
Lp = %s
Npts = %d
xdom = (0.0, Lp,  Npts, periodic=True)
""" % ( Lp, Npts)

# Initialize a simulation object on a mesh
ss = pyrandaSim('advection',imesh)


def interp( pysim, val):

    flux = numpy.zeros_like(val)

    flux[1:-1,:,:] = (val[0:-2,:,:] + val[1:-1,:,:] ) / 2.0
    flux[0 ,:,:] = (val[-1,:,:] + val[0 ,:,:])/2.0
    flux[-1,:,:] = (val[-2,:,:] + val[-1,:,:])/2.0

    return flux

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


def upwind_flux(flux, vel):
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
    up_flux[-1,:,:] = up_flux[0,:,:] #flux[-1,:,:]

    return up_flux

def myDiv(pysim, flux, vel, limiter_type='low'):

    # Low-order flux (upwind) at faces
    loflux = upwind_flux(flux, vel)   # shape (nx+1, ...)

    # High-order flux at faces
    hoflux = numpy.zeros_like(loflux)
    hoflux[:-1,:,:] = pysim.interp_z2fx(flux)  # shape (nx+1, ...)
    hoflux[-1,:,:] = hoflux[0,:,:]
    
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

    # Blend high- and low-order fluxes at interior faces
    flux_limited = numpy.empty_like(hoflux)
    flux_limited[0,:,:] = loflux[0,:,:]
    flux_limited[-1,:,:] = loflux[-1,:,:]
    flux_limited[1:-1,:,:] = loflux[1:-1,:,:] + phi * (hoflux[1:-1,:,:] - loflux[1:-1,:,:])

    #import pdb
    #pdb.set_trace()

    
    # Compute divergence at cell centers
    div = numpy.empty_like(flux)
    #div[:,:,:] = (flux_limited[1:,:,:] - flux_limited[:-1,:,:]) / pysim.dx
    div[1:-1,:,:] = (flux_limited[2:-1,:,:] - flux_limited[1:-2,:,:]) / pysim.dx
    #div[1:-1,:,:] = (flux_limited[1:-2,:,:] - flux_limited[0:-3,:,:]) / pysim.dx
    div[0 ,:,:] = ( flux_limited[ 1,:,:] - flux_limited[0 ,:,:]  ) / pysim.dx
    div[-1,:,:] = ( flux_limited[-1,:,:] - flux_limited[-2,:,:] )  / pysim.dx

    return div


def myDivOld( pysim, flux ):

    iflux = pysim.interp_z2fx( flux )
    #iflux = interp( pysim, flux )

    div = numpy.zeros_like( flux )

    #div[1:-1,:,:] = ( iflux[1:-1,:,:] - iflux[0:-2,:,:] ) / pysim.dx
    div[1:-1,:,:] = ( iflux[2:,:,:] - iflux[1:-1,:,:] ) / pysim.dx

    
    
    # BCs - Periodic
    #div[0 ,:,:] = ( iflux[ 0,:,:] - iflux[-1,:,:]  ) / pysim.dx
    #div[-1,:,:] = ( iflux[-1,:,:] - iflux[-2,:,:] ) / pysim.dx
    div[0 ,:,:] = ( iflux[ 1,:,:] - iflux[0,:,:]  ) / pysim.dx
    div[-1,:,:] = ( iflux[ 0,:,:] - iflux[-1,:,:] ) / pysim.dx

    return div


def ddxlow( pysim, var ):
    dx = pysim.dx
    return (numpy.roll(var, -1, axis=0) - numpy.roll(var, 1,axis=0)) / (2 * dx)
    


ss.addUserDefinedFunction('fvd',myDiv)
ss.addUserDefinedFunction('ddxlow',ddxlow)
ss.addUserDefinedFunction('interp',interp)



# Define the equations of motion
ss.EOM("""
:flux: = ddx( :phi: )
:flux2: = fvd( :phi:, :c: )
:flux3: = ddxlow( :phi: )
:iphi: = iz2fx( :phi: )
#:iphi: = interp( :phi: )
ddt(:phi:)  =  -:c: * :flux2: 
""")


# Initialize variables
ic = """
r   = sqrt( (meshx-pi)**2  )
:phi: = 1.0 + 0.1 * exp( -(r/(pi/4.0))**2 )
:phi2: = 1.0*:phi:
:c:   = 3d(1.0)
"""
ss.setIC(ic)

#ss.variables["u"].data += 1.0

x  = ss.mesh.coords[0].data
xx =  ss.PyMPI.zbar( x )

# Time step size
v = 1.0
dt_max = v / ss.mesh.nn[0] * L * .90
tt = L/v * 1.0 



# Test new interp routines
rho = ss.variables['phi'].data


irho = ss.interp_z2fx( rho )
#irho = interp( ss,  rho )

drho = ss.ddx( rho )

x  = ss.mesh.coords[0].data
xf = x - ss.dx/2.0

plt.figure(2)
plt.plot(x[:,0,0], rho[:,0,0],'kx-')
#plt.figure(3)
#plt.plot(drho[:,0,0])
#plt.figure(4)
plt.plot(xf[:,0,0], irho[:,0,0] , 'bo--' )
plt.pause(.1)

input('wait')




# Main time loop for physics
dt = dt_max
cnt = 1
time = 0.0
viz = True
while tt > time:

    #raw_input('Pause...')
    
    time = ss.rk4(time,dt)
    dt = min(dt_max, (tt - time) )

    if not test:
        ss.iprint("%s -- %s ----- %s ---- %s" % (cnt,time,ss.variables['phi'].data.max() , ss.variables['flux'].sum()  ) )

    # Plot animation of advection
    cnt += 1
    if viz:
        phi = ss.PyMPI.zbar( ss.variables['phi'].data )
        iphi = ss.PyMPI.zbar( ss.variables['iphi'].data )        
        f1 = ss.PyMPI.zbar( ss.variables['flux'].data )
        f2 = ss.PyMPI.zbar( ss.variables['flux2'].data )
        f3 = ss.PyMPI.zbar( ss.variables['flux3'].data )
        if 1: #(ss.PyMPI.master and (cnt%5 == 0)) and (not test):
            plt.figure(1)
            plt.clf()
            plt.plot(xx[:,0],phi[:,0] ,'kx-')
            plt.plot(xf[:,0,0], iphi[:,0], 'bo--')
            plt.figure(2)
            plt.clf()
            plt.plot(xx[:,0],f1[:,0] ,'k.-',label='Native 10th FD')
            plt.plot(xx[:,0],f2[:,0] ,'b.-',label='HO interp')
            plt.plot(xx[:,0],f3[:,0] ,'g.-',label='LO FD')
            plt.legend()
            plt.pause(.001)
            input("pause")


phi = ss.variables['phi'].data
phi2 = ss.variables['phi2'].data
error = numpy.sum( (phi-phi2)**2  )
ss.iprint( error ) 
            


