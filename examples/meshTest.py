import numpy

import matplotlib.pyplot as plt  # noqa: F401


def zoomMesh(npts, x1, xn, xa, xb, tt, dxc, dxf):

    x = numpy.zeros((npts))

    x[0] = x1
    dx = dxc
    for i in range(1, npts):
        x[i] = x[i - 1] + dx

        xmin = (xa + xb) / 2.0
        xh = (xb - xa) / 2.0
        xpr = numpy.abs(x[i] - xmin)
        w = 0.5 * (numpy.tanh((xpr - xh) / tt) + 1.0)
        dx = w * dxc + (1.0 - w) * dxf

    return x


def sfunc(npts, x1, xn, xa, xb, tt, dxc, dxf):
    xx = zoomMesh(npts, x1, xn, xa, xb, tt, dxc, dxf)
    error = xx[-1] - xn
    return error


def zoomMesh_solve(npts, x1, xn, xa, xb, tt, dxf):

    # take initial guess at dxc and solve
    nmin = (xb - xa) / dxf
    if nmin > npts:
        print("Not enough points given")
        exit()
    dxc = ((xn - x1) - (xb - xa)) / (npts - nmin)

    # NR solvers
    delta = 1.0001
    f1 = xn
    cnt = 1
    while (abs(f1) > 0.01 * xn) and cnt < 100:
        f1 = sfunc(npts, x1, xn, xa, xb, tt, dxc, dxf)
        f2 = sfunc(npts, x1, xn, xa, xb, tt, dxc * delta, dxf)
        dxc = dxc * (1.0 - f1 / (f2 - f1) * (delta - 1.0))
        cnt += 1

        # print(f1)
        # print(zoomMesh(npts,x1,xn,xa,xb,tt,dxc,dxf)[-1])
    return zoomMesh(npts, x1, xn, xa, xb, tt, dxc, dxf)


npts = 100
x1 = 0
xn = 1
xa = 0.4
xb = 0.6
tt = 0.02
dxf = (xn - x1) / float(npts) * 0.75


# xx = zoomMesh_solve(npts,x1,xn,xa,xb,tt,dxf)
# plt.plot( xx[:-1], numpy.diff(xx) )
# plt.figure()
# plt.plot( xx )
# plt.show()
