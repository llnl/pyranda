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

ss = pyrandaSim("sod", "xdom=(0.0,6.0,200)")

ss.addPackage(pyrandaBC(ss))

# Define the equations of motion
eom = """
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)
ddt(:rhou:) =  -ddx(:rhou:*:u: + :p: - :tau:)
ddt(:Et:)   =  -ddx( (:Et: + :p: - :tau:)*:u: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u:) ) * ( :gamma: - 1.0 )
# Artificial bulk viscosity (old school way)
:div:       =  ddx(:u:) 
:beta:      =  gbar( ring(:div:) * :rho:) * 7.0e-2
:tau:       =  :beta:*:div:
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

# Approx a max dt and stopping time
dt = 1.0 / 200.0 * 0.75
tt = 1.5

# Start time loop
time = 0.0
while tt > time:
    time = ss.rk4(time, dt)

ss.plot.figure(1)
ss.plot.plot("rho", "b.-")
