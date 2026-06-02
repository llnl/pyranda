from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep
import numpy as np

name = 'shockbubble'

# Fluid and ambient properties
rhoA = 0.138
rhoB = 1.0
gammaA = 1.67
gammaB = 1.4
mwA = 0.138
mwB = 1.0
p0 = 1.0   # Ambient pressure
Ma = 1.2   # Shock Mach number
st = 10    # smearing thickness in cells
eps = 1e-8 # Epsilon to avoid dividing by zero

# Domain parameters
xl = -2.0
xh = 4.0
yl = -1.0
yh = 1.0
Nx = 1200
Ny = 400

# Bubble parameters
r0 = 25.0 / 89.0
bx0 = -0.5
by0 =  0.0

# Shock location
x0_shock = -1.0

# Calculate normal shock relations
a0 = float(np.sqrt(gammaB * p0 / rhoB))
p_shock = p0 * (2.0 * gammaB * Ma**2.0 - (gammaB - 1.0)) / (gammaB + 1.0)
rho_shock = rhoB * (gammaB + 1.0) * Ma**2.0 / ((gammaB - 1.0)*Ma**2.0 + 2.0)
u_shock = a0 * (2.0 / (gammaB + 1.0)) * (Ma - 1.0 / Ma)
igammaB = 1.0 / (gammaB - 1.0)
Et_shock = p_shock*igammaB + 0.5*rho_shock*u_shock**2

print("p_shock = ", p_shock)
print("rho_shock = ", rho_shock)
print("u_shock = ", u_shock)

# Define the mesh
mesh_options = {}
mesh_options['dim'] = 2
mesh_options['x1'] = [xl, yl, 0.0] # Left/Bottom/Back
mesh_options['xn'] = [xh, yh, 1.0] # Right/Top/Front
mesh_options['nn'] = [Nx, Ny, 1]  # Nx/Ny/Nz
mesh_options['periodic'] = [False, False, False]

dx = (xh - xl)/Nx

# Initialize the pyranda sim with the mesh options
ss = pyrandaSim(name, mesh_options)

# Add packages to simulation object
ss.addPackage(pyrandaBC(ss))        # Allows for bc.* functions
ss.addPackage(pyrandaTimestep(ss))  # Allows for "dt.*" functions

eom = f"""
# Primary Equations of Motion (4-eqn model)
ddt(:rhoYA:) =  -div(:rhoYA:*:u: - :Jx: , :rhoYA:*:v: - :Jy:  )
ddt(:rhoYB:) =  -div(:rhoYB:*:u: + :Jx: , :rhoYB:*:v: + :Jy:  )
ddt(:rhou:)  =  -div(:rhou:*:u: - :tauxx:, :rhou:*:v: - :tauxy:)
ddt(:rhov:)  =  -div(:rhov:*:u: - :tauxy:, :rhov:*:v: - :tauyy:)
ddt(:Et:)    =  -div( (:Et: - :tauxx:)*:u: - :tauxy:*:v: - :tx:*:kappa:, (:Et: - :tauyy:)*:v: - :tauxy:*:u: - :ty:*:kappa: )
# Conservative filter of the EOM
:rhoYA: = fbar(:rhoYA:)
:rhoYB: = fbar(:rhoYB:)
:rhou: = fbar(:rhou:)
:rhov: = fbar(:rhov:)
:Et: = fbar(:Et:)
# Primitive variables, mixture rules, and EOS
:rho: = :rhoYA: + :rhoYB:
:YA: = :rhoYA: / :rho:
:YB: = :rhoYB: / :rho:
:u: = :rhou: / :rho:
:v: = :rhov: / :rho:
:gamma: = :YA: * gammaA + :YB: * gammaB
:p: = (:Et: - 0.5*:rho:*(:u:*:u: + :v:*:v:)) * (:gamma: - 1.0)
:mw: = 1.0 / (:YA: / mwA + :YB: / mwB)
:R: = 1.0 / :mw:
:cp: = :R: / (1.0 - 1.0 / :gamma:)
:cv: = :cp: - :R:
:T: = :p: / (:rho: * :R:)
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
:Dsgs:      =  ring(:YA:) * 2.0e-4
:Ysgs:      =  1.0e1*(abs(:YA:) - 1.0 + abs(1.0-:YA: ) )*gridLen**2
:adiff:     =  gbar( :rho:*numpy.maximum(:Dsgs:,:Ysgs:) / :dt: )
[:Yx:,:Yy:,:Yz:] = grad( :YA: )
:Jx:        =  :adiff:*:Yx:
:Jy:        =  :adiff:*:Yy:
# Apply BCs
# left (x1): constant post-shock inflow
bc.const(['rhoYA'], ['x1'], 0.0)
bc.const(['rhoYB'], ['x1'], rho_shock)
bc.const(['rhou'], ['x1'], rho_shock*u_shock)
bc.const(['rhov'], ['x1'], 0.0)
bc.const(['Et'], ['x1'], Et_shock)
# outflow / far-field elsewhere
bc.extrap(['rhoYA','rhoYB','rhou','rhov','Et'], ['xn','y1','yn'])
# Time step computations
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
:dtY: = 0.2 * dt.diff(:adiff:,:rho:)
:dt: = numpy.minimum(:dt:,:dtY:)
:umag: = sqrt( :u:*:u: + :v:*:v: )
"""

EOMDict = {}
EOMDict['mwA'] = mwA
EOMDict['mwB'] = mwB
EOMDict['gammaA'] = gammaA
EOMDict['gammaB'] = gammaB
EOMDict['rho_shock'] = rho_shock
EOMDict['u_shock'] = u_shock
EOMDict['Et_shock'] = Et_shock

# Add the EOM to the solver
ss.EOM(eom, EOMDict)

# Define the initial condition
ic = f"""
s1 = (1 - 2*eps) * (0.5 * (1.0 - tanh( (meshx - (Vx0_shock)) / (Vst*Vdx) ))) + eps
r  = sqrt( (meshx - Vbx0)**2 + (meshy - Vby0)**2 )
s2 = (1 - 2*eps) * (0.5 * (1.0 - tanh( (r - Vr0) / (Vst*Vdx) ))) + eps
:rhoYA: = VrhoA * s2
:rhoYB: = (VrhoB * (1 - s1) + Vrho_shock * s1) * (1.0 - s2)
:rho: = :rhoYA: + :rhoYB:
:YA: = :rhoYA: / :rho:
:YB: = :rhoYB: / :rho:
:rhoA: = :rhoYA: / :YA:
:rhoB: = :rhoYB: / :YB:
# Mixed EOS (partial pressures method)
:gamma: = :YA:* VgammaA + :YB: * VgammaB
:mw: = 1.0 / ( :YA: / VmwA + :YB: / VmwB )
:R: = 1.0 / :mw:
:cp: = :R: / (1.0 - 1.0/:gamma: )
:cv: = :cp: - :R:
# Velocity (post-shock behind the front, quiescent elsewhere)
:u: = Vu_shock * s1
:v: = 0.0
:rhou: = :rho: * :u:
:rhov: = :rho: * :v:
# Pressure (post-shock behind front, ambient ahead; bubble is pressure-matched)
:p: = Vp0*(1.0 - s1) + Vp_shock*s1
# Total energy and temperature
:Et: = :p: * (1.0 / (:gamma: - 1.0))+ 0.5*:rho:*(:u:*:u: + :v:*:v:)
:T: = :p: / (:rho: * :R:)
# Sound speed and initial time step size
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
"""

icDict = {}
icDict['Vx0_shock'] = x0_shock
icDict['Vu_shock'] = u_shock
icDict['Vp_shock'] = p_shock
icDict['Vrho_shock'] = rho_shock
icDict['Vp0'] = p0
icDict['Vst'] = st
icDict['Vdx'] = dx
icDict['Vr0'] = r0
icDict['VrhoA'] = rhoA
icDict['VrhoB'] = rhoB
icDict['Vbx0'] = bx0
icDict['Vby0'] = by0
icDict['VmwA'] = mwA
icDict['VmwB'] = mwB
icDict['VgammaA'] = gammaA
icDict['VgammaB'] = gammaB
icDict['eps'] = eps

# Set the initial condition
ss.setIC(ic, icDict, verbose=True)

time = 0.0 # Start time
tt = 3.0   # Stop time
out_vars = ['YA','rhoYA','rhoYB','u','v','p']

# Start time loop
dump_freq = 200
CFL = 1.0
dt = ss.variables['dt'].data * CFL*.8

viz_freq = tt / 200.0
viz_dump = viz_freq

def myIO():
    ss.write(out_vars)

myIO()

while time < tt :

    # Update the EOM and get next dt
    time = ss.rk4(time,dt)
    dt = min( ss.variables['dt'].data * CFL, dt*1.1)
    dt = min(dt, (tt - time) )

    ss.iprint("(%.2f%%) Cycle: %5d --- Time: %10.4e --- deltat: %10.4e" % (100*time/tt,ss.cycle,time,dt)  )

    # Constant time
    if time > viz_dump:
        viz_dump += viz_freq
        myIO()

    if (ss.cycle%dump_freq == 0) :
        ss.writeRestart()

myIO()
ss.writeRestart()

