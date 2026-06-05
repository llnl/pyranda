import matplotlib.pyplot as plt
from pyranda import pyrandaSim

# Define the domain/mesh
domain = "xdom = (0.0 , 1.0 , 100 )"

# Initialize a simulation object on a mesh
pysim = pyrandaSim("advection", domain)

# Define the equations of motion
pysim.EOM(" ddt(:phi:)  =  -:c: * ddx(:phi:) ")

# Initialize variables
ic = """
:phi: = 1.0 + 0.1 * exp( -(abs(meshx-.5)/.1 )**2 )
:phi0: = :phi:
:c:   = 1.0
"""
pysim.setIC(ic)

# Integrate in time
dt = 0.001
time = 0.0
while time < 0.1:
    time = pysim.rk4(time, dt)


# Plot the initial/final solution
x = pysim.mesh.coords[0].data
phi = pysim.variables["phi"].data
phi0 = pysim.variables["phi0"].data
plt.plot(x[:, 0, 0], phi[:, 0, 0], "k.-")
plt.plot(x[:, 0, 0], phi0[:, 0, 0], "b-")
plt.show()
