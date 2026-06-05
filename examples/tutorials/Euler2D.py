################################################################################
# Copyright (c) 2018, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
#
# Written by: Britton J. Olson, olson45@llnl.gov
################################################################################
from pyranda import pyrandaSim, pyrandaBC
from numpy import pi

# Domain is specified with a "mesh dictionary"
mesh_options = {}
mesh_options["x1"] = [0.0, 0.0, 0.0]
mesh_options["xn"] = [2.0 * pi, 2.0 * pi, 1.0]
mesh_options["nn"] = [64, 64, 1]

ss = pyrandaSim("sod", mesh_options)


ss.addPackage(pyrandaBC(ss))


# Define the equations of motion
eom = """
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)                  - ddy(:rho:*:v:)
ddt(:rhou:) =  -ddx(:rhou:*:u: + :p: - :tau:)   - ddy(:rhou:*:v:)
ddt(:rhov:) =  -ddx(:rhov:*:u:)                 - ddy(:rhov:*:v: + :p: - :tau:)
ddt(:Et:)   =  -ddx( (:Et: + :p: - :tau:)*:u: ) - ddy( (:Et: + :p: - :tau:)*:v: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v:) ) * ( :gamma: - 1.0 )
# Artificial bulk viscosity (old school way)
:div:       =  ddx(:u:) + ddy(:v:)
:beta:      =  gbar(abs(ring(:div:))) * :rho: * 7.0e-2
:tau:       =  :beta:*:div:
# Apply wall BC's
bc.extrap(['rho','Et'],['x1','xn','y1','yn'])
bc.const(['u','v'],['x1','xn','y1','yn'],0.0)
:umag: = sqrt( :u:*:u: + :v:*:v: )
"""

# Add the EOM to the solver
ss.EOM(eom)

# Initial conditions SOD-like problem in 2d
ic = """
rad = sqrt( (meshx-pi)**2  +  (meshy-pi)**2 )
:gamma: = 1.4
:Et:  = gbar( where( rad < pi/2.0, 1.0/(:gamma:-1.0) , .1 /(:gamma:-1.0) ) )
:rho: = gbar( where( rad < pi/2.0, 1.0    , .125 ) )
"""

# Set the initial conditions
ss.setIC(ic)

# Approx a max dt and stopping time
dt = 0.01
tt = 0.75
time = 0.0
while tt > time:
    time = ss.rk4(time, dt)

ss.plot.figure(1)
ss.plot.clf()
ss.plot.contourf("umag", 16)
