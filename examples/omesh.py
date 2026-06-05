import matplotlib.pyplot as plt
import numpy as npy
from scipy import interpolate
import naca


def extend_normal(
    xpts, ypts, dist, sign=1, xcent=0.0, ycent=0.0, feath=False, Rmax=0.0, blend=1.0
):

    npts = xpts.shape[0]
    xnew = npy.zeros(npts)
    ynew = npy.zeros(npts)

    mag = []
    dx = []
    dy = []
    for i in range(npts):
        idx = xpts[(i + 1) % npts] - xpts[i - 1]
        idy = ypts[(i + 1) % npts] - ypts[i - 1]
        mag.append(npy.sqrt(idx * idx + idy * idy))

        dx.append(idx)
        dy.append(idy)

    for i in range(npts):
        idx = dx[i]  # xpts[(i+1)%npts]-xpts[i-1]
        idy = dy[i]  # ypts[(i+1)%npts]-ypts[i-1]
        imag = mag[i]  # npy.sqrt( dx*dx + dy*dy )

        nx = -idy / imag
        ny = idx / imag

        if feath:
            xouter = xcent + Rmax * npy.cos(float(npts - i) / npts * 2.0 * npy.pi)
            youter = ycent + Rmax * npy.sin(float(npts - i) / npts * 2.0 * npy.pi)
            rr = npy.sqrt((xpts[i] - xouter) ** 2 + (ypts[i] - youter) ** 2)
            wgt = min(rr / Rmax, 1.0)
            wgt = blend
            nx = nx * (1.0 - wgt) + wgt * (xouter - xpts[i]) / rr
            ny = ny * (1.0 - wgt) + wgt * (youter - ypts[i]) / rr

        xnew[i] = xpts[i] + nx * dist * sign
        ynew[i] = ypts[i] + ny * dist * sign

    return [xnew, ynew]


def smoothMesh(x, y):

    nx = x.shape[0]
    ny = x.shape[1]

    xnew = x * 1.0
    ynew = y * 1.0
    eps = 0.1
    for i in range(nx):
        # J=0 is Diriclet
        # J=1 to force Neuman

        for j in range(0, ny - 1):
            ip = (i + 1) % nx

            Lx = (x[i - 1, j] + x[ip, j] + x[i, j - 1] + x[i, j + 1]) - 4.0 * x[i, j]
            Ly = (y[i - 1, j] + y[ip, j] + y[i, j - 1] + y[i, j + 1]) - 4.0 * y[i, j]

            wgt = min(float(j) / 5.0, 1.0)

            xnew[i, j] = (
                x[i, j] + wgt * eps * Lx
            )  # ( x[i-1,j]+x[ip,j]+x[i,j-1]+x[i,j+1] ) / 4.0
            ynew[i, j] = (
                y[i, j] + wgt * eps * Ly
            )  # ( y[i-1,j]+y[ip,j]+y[i,j-1]+y[i,j+1] ) / 4.0

    return [xnew, ynew]


def roundTrailingEdge(xnaca, ynaca, teRad):
    X = xnaca
    Y = ynaca
    npts = int((len(X) - 1) / 2)
    gap = 0.0
    te_in = -1
    while gap < teRad:
        te_in += 1
        xl = X[te_in]
        xu = X[2 * npts - te_in]
        yl = Y[te_in]
        yu = Y[2 * npts - te_in]
        gap = npy.sqrt((xu - xl) ** 2 + (yu - yl) ** 2)

    # This chops off the trailing edge
    X = X[te_in:-te_in]
    Y = Y[te_in:-te_in]

    # Now add a parabola for smooth TE
    # USE POLYFIT
    xpoly = []
    xpoly.append(Y[1])
    xpoly.append(Y[0])
    xpoly.append(Y[-1])
    xpoly.append(Y[-2])
    xpoly = npy.array(xpoly)

    ypoly = []
    ypoly.append(X[1])
    ypoly.append(X[0])
    ypoly.append(X[-1])
    ypoly.append(X[-2])
    ypoly = npy.array(ypoly)

    # Rotate to make ortho coords
    theta = -npy.arctan2(-(X[0] - X[-1]), -(Y[0] - Y[-1]))
    R = npy.array([[npy.cos(theta), -npy.sin(theta)], [npy.sin(theta), npy.cos(theta)]])

    xpr = npy.matmul(R, npy.array([xpoly, ypoly]))
    z = npy.polyfit(xpr[0], xpr[1], 2)
    p = npy.poly1d(z)
    xx = npy.linspace(xpr[0][1], xpr[0][-2], npts)
    # xx = [
    yy = p(xx)
    xedge = npy.matmul(npy.linalg.inv(R), [xx, yy])
    xtip = xedge[1]
    ytip = xedge[0]

    # Walk around parabola to find dtheta
    itip = npy.argmax(xtip)
    tipB = [xtip[:itip][::-1], ytip[:itip][::-1]]
    tipT = [xtip[itip:][::-1], ytip[itip:][::-1]]

    # Add bottom to front
    X[:0] = list(tipB[0])[:-1]
    Y[:0] = list(tipB[1])[:-1]

    # Add top to end
    X += list(tipT[0])[1:]
    Y += list(tipT[1])[1:]

    return [X, Y]


def naca_omesh(
    NACA, nx, ny, te=0.05, dratio=10, teSig=40, dr0=0.001, fact=1.05, iS=20, PLOT=True
):

    # NACA = '2412' # NACA airfoil code
    # te = .05      # Radius of Trailing Edge
    # dratio = 10   # TE grid spacing ratio (Coord:TE)
    # teSig = 40    # Denom in guassian (thickness) [Pts/teSig]

    Fpts = nx  # 400    # nx
    Mpts = ny  # 100    # ny
    # dr0 = .001    # dr at airfoil
    # fact = 1.05   # Zoom ratio

    # iS = 20       # Laplacian smooth iterations

    npts = int(nx * dratio)
    X, Y = naca.naca(NACA, npts)

    # Reverse order
    X = X[::-1]
    Y = Y[::-1]

    [X, Y] = roundTrailingEdge(X, Y, te)

    # Wrap around (closed loop)
    X += [X[0]]
    Y += [Y[0]]

    x = npy.array(X)
    y = npy.array(Y)
    tck, u = interpolate.splprep([x, y], s=0)

    # Uniform spacing around foil
    unew = npy.arange(0, 1.0, 1.0 / Fpts)

    # Add stretched meshing.. 0 and 0.5
    wgt = []
    wsum = 0
    for i in range(Fpts):
        wgt.append(wsum)
        dd = min(i, ((Fpts - 1) - i))
        exp0 = npy.exp(-dd * dd / (Fpts / float(teSig)) ** 2)  # 1 @ refine, 0 else
        exp0 *= 1.0 - 1.0 / dratio
        dr = 1.0 - exp0
        wsum += dr

    unew = npy.array(wgt) / wsum
    out = interpolate.splev(unew, tck)

    X = out[0]
    Y = out[1]

    if PLOT:
        plt.figure(1)
        plt.plot(X, Y, "b-")
        plt.axis("equal")

    xgrid = npy.zeros((Fpts, Mpts))
    ygrid = npy.zeros((Fpts, Mpts))
    x1 = npy.array(X)
    y1 = npy.array(Y)
    xgrid[:, 0] = x1
    ygrid[:, 0] = y1
    dr = dr0 * 1.0
    for i in range(1, Mpts):
        # [x1,y1] = extend_normal( x1, y1, dr0 )
        [x1, y1] = extend_normal(
            x1,
            y1,
            dr,
            xcent=0.5,
            ycent=0.0,
            feath=True,
            Rmax=15.0,
            blend=1.0 - dr0 / dr,
        )
        dr *= fact
        xgrid[:, i] = x1
        ygrid[:, i] = y1

    # Smooth the candidate mesh
    for i in range(iS):
        [xgrid, ygrid] = smoothMesh(xgrid, ygrid)

    if PLOT:
        plt.plot(xgrid, ygrid, "k-", lw=0.5, alpha=0.5)
        plt.plot(xgrid.T, ygrid.T, "k-", lw=0.5, alpha=0.5)
        plt.axis("equal")
        plt.show()

    return [xgrid, ygrid]


if __name__ == "__main__":
    NACA = "2412"
    nx = 400
    ny = 100
    naca_omesh(NACA, nx, ny)
