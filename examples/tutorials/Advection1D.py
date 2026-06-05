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

from pyranda import pyrandaSim

# Define the domain/mesh
domain = "xdom = (0.0 , 1.0 , 100 )"

# Initialize a simulation object on a mesh
pysim = pyrandaSim("advection", domain)

# Define the equations of motion
pysim.EOM("ddt(:phi:)  =  - ddx(:phi:)")

# Initialize variables
pysim.setIC(":phi: = 1.0 + 0.1 * exp( -(abs(meshx-.5)/.1 )**2 )")

# Integrate in time
dt = 0.001
time = 0.0
while time < 0.1:
    time = pysim.rk4(time, dt)

# Plot solution
pysim.plot.plot("phi")
