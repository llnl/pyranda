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
import numpy


class pyrandaPlot:
    def __init__(self, pyranda):

        self.pyranda = pyranda  # point to main obj

        # Onlt import on master
        if self.pyranda.PyMPI.master:
            import matplotlib as mpl
            import matplotlib.pyplot as plt
            import matplotlib.style

            mpl.style.use("classic")

    # Some master overloaded functions
    def figure(self, fignum):
        if self.pyranda.PyMPI.master:
            import matplotlib.pyplot as plt

            plt.figure(fignum)

    def show(self):
        if self.pyranda.PyMPI.master:
            import matplotlib.pyplot as plt

            plt.show()

    def clf(self):
        if self.pyranda.PyMPI.master:
            import matplotlib.pyplot as plt

            plt.clf()

    def title(self, name):
        if self.pyranda.PyMPI.master:
            import matplotlib.pyplot as plt

            plt.title(name)
            plt.pause(0.01)

    def pause(self, val=0.01):
        if self.pyranda.PyMPI.master:
            import matplotlib.pyplot as plt

            plt.pause(val)

    def plot(self, var, style="k-", slice2d="j=0;k=0", **kwargs):
        """
        Plots 1D line plots
        """
        nx = self.pyranda.nx
        ny = self.pyranda.ny
        nz = self.pyranda.nz

        # 1D - figure out directions
        if ny == nz == 1:
            slice2d = "j=0;k=0"
        if ny == nx == 1:
            slice2d = "j=0;i=0"
        if nx == nz == 1:
            slice2d = "i=0;k=0"

        vv = self.pyranda.variables[var].data
        vdata = self.getLine(vv, slice2d)
        xdata = self.getGrid1d(slice2d)

        if self.pyranda.PyMPI.master:
            import matplotlib.pyplot as plt

            plt.plot(xdata, vdata, style, **kwargs)
            plt.title(
                "Plot of %s (t=%4.4e and cycle=%d)"
                % (var, self.pyranda.time, self.pyranda.cycle)
            )
            plt.pause(0.01)

    def contour(self, var, levels, slice3d="k=0", **kwargs):
        import matplotlib.pyplot as plt

        self.contourf(var, levels, slice3d=slice3d, filled=False, **kwargs)

    def contourf(
        self, var, levels, slice3d="k=0", filled=True, noPause=False, **kwargs
    ):
        """
        Plots 2d contour plot.  Only works for 2d and 3d problems.
        2d - figure our directions and plot
        3d - assumes slice at k=0 unless specified
        """
        nx = self.pyranda.nx
        ny = self.pyranda.ny
        nz = self.pyranda.nz
        # 1D - Reject
        if (ny == nz == 1) or (nx == nz == 1) or (ny == nx == 1):
            self.pyranda.iprint("pyranda contourf only works for 2D & 3D problems")

        # 2D - figure out directions
        if nz == 1:
            slice3d = "k=0"
        if ny == 1:
            slice3d = "j=0"
        if nx == 1:
            slice3d = "i=0"

        # Get the sliced data
        vv = self.pyranda.variables[var].data
        vdata = self.getSlice(vv, slice3d)
        [xdata, ydata] = self.getGrid2d(slice3d)

        if self.pyranda.PyMPI.master:
            import matplotlib.pyplot as plt

            if filled:
                plt.contourf(xdata, ydata, vdata, levels, **kwargs)
                plt.colorbar()
            else:
                if not (type(levels) == type([])):
                    levels = numpy.linspace(vdata.min(), vdata.max(), levels)
                plt.contour(xdata, ydata, vdata, levels, **kwargs)
            plt.axis("equal")
            plt.title(
                "Contours of %s (t=%4.4e and cycle=%d)"
                % (var, self.pyranda.time, self.pyranda.cycle)
            )
            if not noPause:
                plt.pause(0.01)

    def showGrid(self, slice3d="k=0"):

        [xx, yy] = self.getGrid2d(slice3d)

        if self.pyranda.PyMPI.master:
            import matplotlib.pyplot as plt

            plt.plot(xx, yy, "k-", lw=0.5, alpha=0.5)
            plt.plot(xx.T, yy.T, "k-", lw=0.5, alpha=0.5)
            plt.axis("equal")
            plt.pause(0.01)

    def getGrid2d(self, slice3d, om=None):

        nx = self.pyranda.nx
        ny = self.pyranda.ny
        nz = self.pyranda.nz

        x = self.pyranda.mesh.coords[0].data[:, :, :]
        y = self.pyranda.mesh.coords[1].data[:, :, :]
        z = self.pyranda.mesh.coords[2].data[:, :, :]

        # 2D - figure out directions
        if nz == 1:
            slice3d = "k=0"
        if ny == 1:
            slice3d = "j=0"
        if nx == 1:
            slice3d = "i=0"

        # get directions
        if "k" in slice3d:
            xx = self.getSlice(x, slice3d)
            yy = self.getSlice(y, slice3d)
        if "j" in slice3d:
            xx = self.getSlice(x, slice3d)
            yy = self.getSlice(z, slice3d)
        if "i" in slice3d:
            xx = self.getSlice(y, slice3d)
            yy = self.getSlice(z, slice3d)

        return [xx, yy]

    def getSlice(self, data, slice3d):

        nx = self.pyranda.nx
        ny = self.pyranda.ny
        nz = self.pyranda.nz

        exec(slice3d, globals())

        # get directions - k slice
        if "k" in slice3d:
            gdata = self.pyranda.PyMPI.getIJK(data, [0, nx], [0, ny], [k, k + 1])[
                :, :, 0
            ]
        if "j" in slice3d:
            gdata = self.pyranda.PyMPI.getIJK(data, [0, nx], [j, j + 1], [0, nz])[
                :, 0, :
            ]
        if "i" in slice3d:
            gdata = self.pyranda.PyMPI.getIJK(data, [i, i + 1], [0, ny], [0, nz])[
                0, :, :
            ]

        return gdata

    def getGrid1d(self, slice2d):

        nx = self.pyranda.nx
        ny = self.pyranda.ny
        nz = self.pyranda.nz

        x = self.pyranda.mesh.coords[0].data
        y = self.pyranda.mesh.coords[1].data
        z = self.pyranda.mesh.coords[2].data

        # 1D - figure out directions
        if ny == nz == 1:
            slice2d = "j=0;k=0"
        if ny == nx == 1:
            slice2d = "j=0;i=0"
        if nx == nz == 1:
            slice2d = "i=0;k=0"

        # get directions
        if ("k" in slice2d) and ("j" in slice2d):
            xx = self.getLine(x, slice2d)
        if ("j" in slice2d) and ("i" in slice2d):
            xx = self.getLine(z, slice2d)
        if ("i" in slice2d) and ("k" in slice2d):
            xx = self.getLine(y, slice2d)

        return xx

    # ruff :noqa (i, j, k defined implicitly)
    def getLine(self, data, slice2d):

        nx = self.pyranda.nx
        ny = self.pyranda.ny
        nz = self.pyranda.nz

        exec(slice2d, globals())

        # get directions - k slice
        if ("k" in slice2d) and ("j" in slice2d):
            gdata = self.pyranda.PyMPI.getIJK(data, [0, nx], [j, j + 1], [k, k + 1])[
                :, 0, 0
            ]
        if ("j" in slice2d) and ("i" in slice2d):
            gdata = self.pyranda.PyMPI.getIJK(data, [i, i + 1], [j, j + 1], [0, nz])[
                0, 0, :
            ]
        if ("i" in slice2d) and ("k" in slice2d):
            gdata = self.pyranda.PyMPI.getIJK(data, [i, i + 1], [0, ny], [k, k + 1])[
                0, :, 0
            ]

        return gdata
