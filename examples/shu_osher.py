import sys
import numpy
import matplotlib.pyplot as plt

from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep


# Try to get args
try:
    test = bool(int(sys.argv[1]))
except Exception:
    test = False


## Define a mesh


def getShu(Npts, viz=False):

    L = 10.0
    Lp = L * (Npts - 1.0) / Npts
    imesh = "xdom = (-Lp/2.0, Lp/2.0, Npts)"
    imesh = imesh.replace("Lp", str(Lp))
    imesh = imesh.replace("Npts", str(Npts))

    # Initialize a simulation object on a mesh
    ss = pyrandaSim("sod", imesh)
    ss.addPackage(pyrandaBC(ss))
    ss.addPackage(pyrandaTimestep(ss))

    # Define the equations of motion
    from equation_library import euler_1d

    eom = euler_1d
    eom += """
# Apply constant BCs
bc.const(['u'],['xn'],0.0)
bc.const(['u'],['x1'],2.629369)
bc.const(['rho'],['x1'],3.857143)
bc.const(['p'],['x1'],10.33333)
bc.const(['Et'],['x1'],39.166585)
:rhou: = :rho:*:u:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,0.0,0.0,:cs:)
:dt: = numpy.minimum(:dt:,0.2 * dt.diff(:beta:,:rho:))
"""

    # Add the EOM to the solver
    ss.EOM(eom)

    # Initial conditions Shu-Osher test problem
    ic = """
:gamma: = 1.4
eps = 2.0e-1
:tmp: = (meshx+4.0)/%s
:dum: = tanh(:tmp:)
:dum: = (:dum:+1.0)/2.0
:rho: = 3.857143 + :dum:*(1.0+eps*sin(5.0*meshx) - 3.857143)
:u: = 2.629369*(1.0-:dum:)
:p: = 10.33333 + :dum:*(1.0-10.33333)
:rhou: = :rho: * :u:
:Et:  = :p:/(:gamma: -1.0) + .5*:rho:*:u:*:u:
"""

    # Set the initial conditions
    ss.setIC(ic % ss.dx)

    # Write a time loop
    time = 0.0

    # Approx a max dt and stopping time
    CFL = 0.5

    dt = ss.variables["dt"].data * CFL * 0.01

    # Mesh for viz on master
    x = ss.mesh.coords[0].data
    xx = ss.PyMPI.zbar(x)

    # Start time loop
    cnt = 1
    viz_freq = 50
    pvar = "rho"
    # viz = True
    tt = 1.8
    while tt > time:
        # Update the EOM and get next dt
        time = ss.rk4(time, dt)
        dt = min(dt * 1.1, ss.variables["dt"].data * CFL)
        dt = min(dt, (tt - time))

        # Print some output
        ss.iprint("%s -- %s" % (cnt, time))
        cnt += 1
        v = ss.PyMPI.zbar(ss.variables[pvar].data)
        if viz:
            if (ss.PyMPI.master and (cnt % viz_freq == 0)) and True:
                # raw_input('Poop')
                plt.figure(1)
                plt.clf()
                plt.plot(xx[:, 0], v[:, 0], "k.-")
                plt.title(pvar)
                plt.pause(0.001)

    return [ss, xx[:, 0], v[:, 0]]


Npts = [100, 200, 400, 800]
solution = []
for Npt in Npts:
    [ss, x, rho] = getShu(Npt)
    solution.append(numpy.interp(1.85, x, rho))
    if not test:
        plt.plot(x, rho, label="N=%s" % Npt)

solution = numpy.array(solution)
diff = numpy.abs(solution[1:] - solution[:-1])
cr1 = numpy.log(diff[-1] / diff[-2]) / numpy.log(0.5)
cr2 = numpy.log(diff[-2] / diff[-3]) / numpy.log(0.5)

print((cr1 + cr2) / 2.0)

if not test:
    # Make a tex file
    ss.setupLatex()

    ss.latex.tMap[":rho:"] = r"\rho "
    ss.latex.tMap[":rhou:"] = r"\rho u"
    ss.latex.tMap[":Et:"] = r"E_t "

    intro = ss.latex.addSection("Introduction")
    intro.body = """
    This brief document shows the ability to integrate the simulation with formal 
    documentation.  The advection equations is solved and plotted below.
    """

    equations = ss.latex.addSection("Equations of Motion")
    equations.body = """
    Pyranda solves the following equations using a 10th order compact finite difference
    stencil and a 5 stage RK4 temporal integration scheme.
    """
    equations.body += ss.latex.renderEqu()

    equations.body += """
    where terms are closed following:
    """

    equations.body += ss.latex.renderEqu("ALG")

    results = ss.latex.addSection("Results")
    plt.figure(1)
    results.addFigure("Profiles", size=0.45)
    results.body = "Here in the figure, we see the grid convergence"

    ss.latex.makeDoc()
    ss.latex.showPDF()

    plt.legend()
    plt.show()
