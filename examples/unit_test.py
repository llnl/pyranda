import numpy
from pyranda import pyrandaSim

ERROR = 0.0
##############################
########  DERIVATIVES ########
##############################

# Define the domain/mesh
Nx = 32
L = numpy.pi * 2.0 * (Nx - 1) / Nx
domain = """
xdom = (0.0 , L , Nx ,periodic=True)
ydom = (0.0 , L , Nx ,periodic=True)
zdom = (0.0 , L , Nx ,periodic=True)
""".replace("L", str(L)).replace("Nx", str(Nx))

# Periodic function in X-Y-Z
pysim = pyrandaSim("advection", domain)

expr = "sum( abs( ddx( sin(meshx) ) - cos(meshx) ))"
errorX = pysim.eval(expr)
errorY = pysim.eval(expr.replace("x", "y"))
errorZ = pysim.eval(expr.replace("x", "z"))
ERROR += errorX + errorY + errorZ

print("Error X = %s" % errorX)
print("Error Y = %s" % errorY)
print("Error Z = %s" % errorZ)

# Non-periodic case
domain = """
xdom = (0.0 , 1.0 , Nx )
ydom = (0.0 , 1.0 , Nx )
zdom = (0.0 , 1.0 , Nx )
""".replace("Nx", str(Nx))
pysim2 = pyrandaSim("nonperiodic", domain)

# Zero
print("Zero ")
expr = "mean( abs(ddx(meshx*0.0)) )"
errorX = pysim2.eval(expr)
errorY = pysim2.eval(expr.replace("x", "y"))
errorZ = pysim2.eval(expr.replace("x", "z"))
ERROR += errorX + errorY + errorZ

print("Error X = %s" % errorX)
print("Error Y = %s" % errorY)
print("Error Z = %s" % errorZ)

# Constant
print("Constant ")
expr = "mean( abs(ddx(meshx*0.0 + 1.0)) )"
errorX = pysim2.eval(expr)
errorY = pysim2.eval(expr.replace("x", "y"))
errorZ = pysim2.eval(expr.replace("x", "z"))
ERROR += errorX + errorY + errorZ

print("Error X = %s" % errorX)
print("Error Y = %s" % errorY)
print("Error Z = %s" % errorZ)

# Powers
powers = [1, 2, 3, 4]
expr = "mean( abs(ddx(meshx**power) - power*meshx**(power-1) )  )"
for power in powers:
    print("Power = %s" % power)
    exp = expr.replace("power", str(float(power)))
    errorX = pysim2.eval(exp)
    errorY = pysim2.eval(exp.replace("x", "y"))
    errorZ = pysim2.eval(exp.replace("x", "z"))
    ERROR += errorX + errorY + errorZ
    print("Error X = %s" % errorX)
    print("Error Y = %s" % errorY)
    print("Error Z = %s" % errorZ)


##############################
##########  FILTERS ##########
##############################
print("Filters-fbar")
# Zero
print("Zero ")
expr = "mean( abs(fbar(meshx*0.0)) )"
errorX = pysim2.eval(expr)
print("Error X = %s" % errorX)
ERROR += errorX

# Constant
print("Constant ")
expr = "mean( abs(fbar(meshx*0.0 + 1.0) - (meshx*0.0 + 1.0)  ) )"
errorX = pysim2.eval(expr)
print("Error X = %s" % errorX)
ERROR += errorX

# Linear
print("Linear ")
expr = "mean( abs(fbar(meshx) - (meshx)  ) )"
errorX = pysim2.eval(expr)
print("Error X = %s" % errorX)
ERROR += errorX

# Periodic
print("sine wave")
expr = "mean( abs(fbar(sin(meshx)) - sin(meshx)  ) )"
pysim.PyMPI.setPatch()
errorX = pysim.eval(expr)
print("Error X = %s" % errorX)
ERROR += errorX

print("periodic constant")
expr = "mean( abs(fbar( meshx*0.0 + 1.0) - 1.0  ) )"
errorX = pysim.eval(expr)
print("Error X = %s" % errorX)
ERROR += errorX

print("periodic zero")
expr = "mean( abs(fbar( meshx*0.0 )  ) )"
errorX = pysim.eval(expr)
print("Error X = %s" % errorX)
ERROR += errorX


print("Filters-gbar")
# Zero
print("Zero ")
expr = "mean( abs(gbar(meshx*0.0)) )"
pysim2.PyMPI.setPatch()
errorX = pysim2.eval(expr)
print("Error X = %s" % errorX)
ERROR += errorX

# Constant
print("Constant ")
expr = "mean( abs(gbar(meshx*0.0 + 1.0) - (meshx*0.0 + 1.0)  ) )"
errorX = pysim2.eval(expr)
print("Error X = %s" % errorX)
ERROR += errorX

# Linear
print("Linear ")
expr = "mean( abs(gbar(meshx) - (meshx)  ) )"
errorX = pysim2.eval(expr)
print("Error X = %s" % errorX)
ERROR += errorX

# Periodic
print("sine wave")
expr = "mean( abs(gbar(sin(meshx)) - sin(meshx)  ) )"
pysim.PyMPI.setPatch()
errorX = pysim.eval(expr)
print("Error X = %s" % errorX)
ERROR += errorX

print("periodic constant")
expr = "mean( abs(gbar( meshx*0.0 + 1.0) - 1.0  ) )"
errorX = pysim.eval(expr)
print("Error X = %s" % errorX)
ERROR += errorX

print("periodic zero")
expr = "mean( abs(gbar( meshx*0.0 )  ) )"
errorX = pysim.eval(expr)
print("Error X = %s" % errorX)
ERROR += errorX


# fid = open('unit_test.out','w')
# fid.write(out)
# fid.close()
print(ERROR)
