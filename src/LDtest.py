#!/usr/bin/env python3
'''
Program: Practice Lagrangian Descriptors with model systems:
           1. Analytical potential from CPL 2020, 754, 137610.
           2. Fitting potential of methoxyl radical isomerization
           3. Fitting potential of NCH1 (asymmeric PTSB) system

History
    2021/01/12, Grace, generate trajectories
'''

# import libraries
import sys
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import random

# set up Global Variables
# Boundary condition
ymin = -1.5
ymax = 1.5
# reduced mass
m_x = 1
m_y = 1
totTraj = 1  # total amount of trajectories
totStep = 10  # total time steps
Dt = 0.5


def main():
    # 1. Define potential energy function
    # PES(): 1 = analytical PES, 2 = fitting PES
    X, Y, V = PES(1)

    # 2. Generate trajectories w/ fixed total energy
    # Constraint
    H = 0.1  # microcanonical ensemble [NVE]
    x = 0.05  # c.f. eqn 5 @ CPL 2020, 754, 137610.

    # traj is a 3-D list;TODO: sampling I.C.
    # traj = genTraj(H, x, totStep, totTraj)

    # 3. LD? what the hell?
    LD()

    # 4. Visualize results
    # 4-1. Plot PES only
    plt = plotV(X, Y, V)

    # 4-2. Add trajectories on the above figure
    plt = plotTraj(plt, traj)

    plt.show()


def BC():
    x = np.linspace(0, 2, 50)
    y = np.linspace(ymin, ymax, 50)
    X, Y = np.meshgrid(x, y)
    return X, Y


def totline(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    f.close()
    return i + 1


def polyfit2d(x, y, z, kx=3, ky=3, order=None):
    # c.f. https://stackoverflow.com/questions/33964913/equivalent-of-polyfit-for-a-2d-polynomial-in-python
    '''
    Two dimensional polynomial fitting by least squares.
    Fits the functional form f(x,y) = z.

    Notes
    -----
    Resultant fit can be plotted with:
    np.polynomial.polynomial.polygrid2d(x, y, soln.reshape((kx+1, ky+1)))

    Parameters
    ----------
    x, y: array-like, 1d
        x and y coordinates.
    z: np.ndarray, 2d
        Surface to fit.
    kx, ky: int, default is 3
        Polynomial order in x and y, respectively.
    order: int or None, default is None
        If None, all coefficients up to maxiumum kx, ky, ie. up to and including x^kx*y^ky, are considered.
        If int, coefficients up to a maximum of kx+ky <= order are considered.

    Returns
    -------
    Return paramters from np.linalg.lstsq.


    soln: np.ndarray
        Array of polynomial coefficients.
    residuals: np.ndarray
    rank: int
    s: np.ndarray

    '''

    # grid coords
    x, y = np.meshgrid(x, y)
    # coefficient array, up to x^kx, y^ky
    coeffs = np.ones((kx+1, ky+1))

    # solve array
    a = np.zeros((coeffs.size, x.size))

    # for each coefficient produce array x^i, y^j
    for index, (j, i) in enumerate(np.ndindex(coeffs.shape)):
        # do not include powers greater than order
        if order is not None and i + j > order:
            arr = np.zeros_like(x)
        else:
            arr = coeffs[i, j] * x**i * y**j
        a[index] = arr.ravel()

    # do leastsq fitting and return leastsq result
    return np.linalg.lstsq(a.T, np.ravel(z), rcond=None)


def V1(X, Y):
    # c.f. eqn. 2 @ CPL 2020, 754, 137610.
    return (8/3)*X**3 - 4*X**2 + (1/2)*Y**2 + X*Y**2*(Y**2-2)


def a(X, Y):
    # c.f. eqn. 3 @ CPL 2020, 754, 137610.
    # FIXME: negative sign?
    ax = (8*X*(1-X)+Y**2*(2-Y**2))/m_x
    ay = (Y*(4*X*(1-Y**2)-1))/m_y
    return ax, ay


def totE(traj):
    # x=traj[0], y=traj[1], px=traj[2], py=traj[3]
    v = V1(traj[0], traj[1])  # CPL 2020, 754, 137610.
    return traj[2]**2/(2*m_x)+traj[3]**2/(2*m_y)+v


def PES(case):

    if case == 1:  # Analytical PES
        # Model 1. Analytical potential from CPL 2020, 754, 137610. (c.f. eqn 2)
        X, Y = BC()  # define boundary of X and Y
        V = V1(X, Y)
    elif case == 2:  # Fitting PES
        # 1. Import rawdata
        # Model 2 and 3. Fitting potential of methoxyl radical isomerization/NCH1 (asymmeric PTSB) system
        tline = totline(str(sys.argv[1]))
        dim = int(math.sqrt(tline))
        x, y, v = np.loadtxt(str(sys.argv[1]), delimiter=' ', unpack=True)
        X = np.reshape(x, (dim, dim))
        Y = np.reshape(y, (dim, dim))
        V = np.reshape(v, (dim, dim))

        # Fitting

    return X, Y, V


def genTraj(H, xi, totStep, totTraj):
    '''
    Reference of sampling: https://pubs.acs.org/doi/abs/10.1021/acs.jpca.8b07205

    Hamilton's eqn of motion; H(q,p,t)
    '''
    traj = list()
    for nTraj in range(totTraj):
        # 1. Initial condition
        vi = 2*H
        # TODO: DS sampling
        while H-vi < 0:
            random.seed(random.randint(0, 100))
            yi = ymin+(random.random()*(ymax-ymin))  # random = [0,1]
            vi = V1(xi, yi)

        # TODO: momentum sampling
        # 1. Assume initial momentum are equal
        pxi = math.sqrt((H - vi) * m_x)
        pyi = pxi
        # initial velocity
        vel_xi = pxi/m_x
        vel_yi = pyi/m_y
        ax, ay = a(xi, yi)

        # SingleTraj.append([xi, yi, pxi, pyi])
        if nTraj == 0:
            traj = [[[xi, yi, pxi, pyi]]]
        else:
            traj.append([[xi, yi, pxi, pyi]])

        for nStep in range(totStep):
            x = xi + vel_xi * Dt
            y = yi + vel_yi * Dt
            vi = V1(x, y)
            # TODO: other way to continue?
            if H-vi < 0:
                break
            vel_x = vel_xi + ax*Dt
            vel_y = vel_yi + ay*Dt
            traj[nTraj].append([x, y, vel_x*m_x, vel_y*m_y])

            # replace the initial value
            ax, ay = a(x, y)
            xi = x
            yi = y
            vel_xi = vel_x
            vel_yi = vel_y
            # print(x, y, vel_x*m_x, vel_y*m_y)
    # m = np.array(traj)
    # print(len(traj))
    # print(m.shape)
    return traj


def LD():
    '''
    Reference: 
    theory: https://www.chemicalreactions.io/prologue/LDs/LDs_victor-jekyll.html
    code:
       1. https://github.com/broncio123/cirque_project_champs/blob/master/notebooks/compute_LD_total_variable_kn_tau.py
       2. https://github.com/broncio123/cirque_project_champs/blob/1d1e58ab29e92a171ec0be19c7017379261179e8/notebooks/pylds/base.py#L272
    '''
    print('test')


def plotV(X, Y, V):
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    ct = plt.contour(X, Y, V,  50, colors='k')
    # plt.axis('square') # square figure
    plt.clabel(ct, inline=1, fmt='%1.0f', fontsize=8)  # value of contour
    return plt


def plotTraj(plt, traj):
    for nTraj in range(totTraj):
        for i in range(totStep):
            plt.plot(list(zip(*traj[nTraj]))[0],
                     list(zip(*traj[nTraj]))[1], 'r-')
    return plt


main()
