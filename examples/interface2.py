import numpy
import matplotlib.pyplot as plt
from pyranda import pyrandaSim, pyrandaIBM


## Define a mesh
Npts = 128
L = numpy.pi * 2.0
Lp = L * (Npts - 1.0) / Npts

mesh_options = {}
mesh_options["type"] = "cartesian"
mesh_options["periodic"] = numpy.array([False, True, True])
mesh_options["dim"] = 3
mesh_options["x1"] = [0.0, 0.0, 0.0]
mesh_options["xn"] = [Lp, Lp, Lp]
mesh_options["nn"] = [Npts, 1, 1]
# mesh_options['nn'] = [ Npts, Npts ,  1  ]


# Initialize a simulation object on a mesh
ss = pyrandaSim("advection", mesh_options)
ss.addPackage(pyrandaIBM(ss))


# Define the equations of motion
eom = """
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)                  - ddy(:rho:*:v:)
ddt(:rhou:) =  -ddx(:rhou:*:u: + :p: - :tau:)   - ddy(:rhou:*:v:)
ddt(:rhoA:)  =  -ddx(:rhoA:*:uA:)  
ddt(:rhoAuA:) =  -ddx(:rhoAuA:*:uA: + :pA: - :tauA:)
ddt(:EtA:)   =  -ddx( (:EtA: + :pA: - :tauA:)*:uA: ) 
ddt(:rhov:) =  -ddx(:rhov:*:u:)                 - ddy(:rhov:*:v: + :p: - :tau:)
ddt(:Et:)   =  -ddx( (:Et: + :p: - :tau:)*:u: ) - ddy( (:Et: + :p: - :tau:)*:v: )
ddt(:phi:)  =  - :gx: * :u1:
#:u1:        =  .001 * sin( simtime * 0.0 )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhoA:      =  fbar( :rhoA:  )
:rhoAuA:    =  fbar( :rhoAuA: )
:EtA:       =  fbar( :EtA:   )
:rhov:      =  fbar( :rhov: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho: * 0.0 + 0.7
:v:         =  :rhov: / :rho: * 0.0 + 0.7
:uA:        =  :rhoAuA: / :rhoA:
# Artificial bulk viscosity (old school way)
:div:       =  ddx(:u:) + ddy(:v:)
:beta:      =  gbar(abs(lap(lap(:div:))))*:dx6: * :rho: * 0.2
:tau:       =  :beta:*:div:
:divA:       =  ddx(:uA:) 
:betaA:      =  gbar(abs(lap(lap(:divA:))))*:dx6: * :rhoA: * 0.2
:tauA:       =  :betaA:*:divA:
# Immersed boundary method used for MM levelset
[:gx:,:gy:,:gz:] = grad( :phi: )
[:u:,:v:,:w:] = ibmV( [:u:,:v:,0.0], :phi:, [:gx:,:gy:,:gz:], [:u1:,0.0,0.0] )
:Et:          = ibmS( :Et:  , :phi:, [:gx:,:gy:,:gz:] )
:rho:         = ibmS( :rho: , :phi:, [:gx:,:gy:,:gz:] )
[:uA:,:v:,:w:] = ibmV( [:uA:,:v:,0.0], -:phi:, [-:gx:,-:gy:,-:gz:], [:u1:,0.0,0.0] )
:EtA:          = ibmS( :EtA:  , -:phi:, [-:gx:,-:gy:,-:gz:] )
:rhoA:         = ibmS( :rhoA: , -:phi:, [-:gx:,-:gy:,-:gz:] )
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v:) ) * ( 1.4 - 1.0 )
:pA:        =  ( :EtA: - .5*:rhoA:*(:uA:*:uA: ) ) * ( 1.4 - 1.0 )
:uT:        =  where( :phi: > 0.0, :u:, :uA:)
:pT:        =  where( :phi: > 0.0, :p:, :pA:)
:rhoB:        =  where( :phi: > 0.0, :rho:, :rhoA:)
:rhoT:        =  where( abs(:phi:/:dx:*2.0) > 1.0, :rhoB:, .5*(:phi:/:dx:*2.0 + 1.0)*(:rho:-:rhoA:) + (:rhoA:)   )
"""
ss.EOM(eom)

#:umag:      =  numpy.sqrt(:u:*:u: + :v:*:v:)

# Initialize variables
x = ss.mesh.coords[0].data
y = ss.mesh.coords[1].data
z = ss.mesh.coords[2].data


rad = numpy.sqrt((x - numpy.pi) ** 2)  # + (y-numpy.pi)**2 ) #+ (z-numpy.pi)**2  )
ss.variables["Et"].data += 2.5
ss.variables["rho"].data += 1.0

ss.variables["EtA"].data += 2.5
ss.variables["rhoA"].data += 10.0e1 * (
    1.0 + 0.1 * numpy.exp(-((x - numpy.pi * 1.5) ** 2) / 0.01)
)

# Init momenta
ss.variables["rhou"].data = ss.variables["rho"].data * 0.70
ss.variables["rhoAuA"].data = ss.variables["rhoA"].data * 0.70
ss.variables["u1"].data = 0.70

ss.variables["dx"].data += x[1, 0, 0] - x[0, 0, 0]
ss.variables["dx6"].data += (x[1, 0, 0] - x[0, 0, 0]) ** 6

# Define the interface location
ss.variables["phi"].data = 3.14159 - x


time = 0.0
viz = True

v = 1.0

dt_max = v / ss.mesh.nn[0] * 1.0

tt = L / v * 0.25  # dt_max

xx = ss.PyMPI.zbar(x)
yy = ss.PyMPI.zbar(y)

dt = dt_max
cnt = 1


while tt > time:
    time = ss.rk4(time, dt)
    dt = min(dt_max, (tt - time))

    ss.iprint("%s -- %s" % (cnt, time))
    cnt += 1
    if viz:
        # v1 = ss.PyMPI.zbar( ss.variables['rhoA'].data )
        v2 = ss.PyMPI.zbar(ss.variables["rhoT"].data)
        if ss.PyMPI.master and (cnt % 5 == 0) and True:
            plt.figure(2)
            plt.clf()
            # plt.contourf( xx,yy,phi ,64 )
            # plt.plot(xx[:,0],v1[:,0] ,'b.-')
            plt.plot(xx[:, 0], v2[:, 0], "k.-")
            plt.pause(0.001)
            # import pdb
            # pdb.set_trace()
