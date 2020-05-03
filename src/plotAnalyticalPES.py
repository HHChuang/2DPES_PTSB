#!/usr/bin/env python3
#############################################################################
# Program:                                                                  #
#   use analytical 2D-PES to plot PTSB model potential.                     #
#                                                                           #
# History:                                                                  #
#   2019/08/16, Grace                                                       #
#   2020/04/24, Grace, add sympy (calculate analytical gradient/Hessian)    #
#   2020/04/25, Grace, formulate Hessian matrix and then diagonalize it.    #
#   2020/04/27, Grace, formulate gradient descent path; (IRC/IRP).          #
#   2020/04/30, Grace, add contour line and modify stationaryPts().         #
#                                                                           #
#############################################################################
import sys
import numpy as np
import sympy as sp
from sympy.solvers.solveset import solveset, solveset_real, nonlinsolve
import math
# from sympy.plotting.plot import plot_contour
import matplotlib.pyplot as plt


def sp_symPES(X, Y):
    # Theor Chem Acc 1998, 100, 285-299.
    # https://f1000.com/work/item/3809905/resources/2861282/pdf
    # sp_symPES = 2*Y + Y**2 + X**2*(0.4*X**2+Y)

    # Theor Chem Acc 2004, 112, 40-51.
    # https://f1000.com/work/item/5637171/resources/4634218/pdf
    sp_symPES = (1/2.0) * (X*Y**2 - Y*X**2 - 2.0*X + 2.0*Y) \
        + (1/30.0) * (X**4 + Y**4)

    #
    # sp_symPES = 44730.4129 \
    #     - 66786.5363 * sp.cos(Y) \
    #     + 26352.6908 * sp.cos(2*Y) \
    #     - 3117.3613 * sp.cos(4*Y) \
    #     - 659.3217 * sp.cos(6*Y) \
    #     + 621.9640 * sp.sin(3*X) * sp.sin(Y) \
    #     - 138.3050 * sp.sin(3*X) * sp.sin(2*Y) \
    #     - 111.5488 * sp.cos(8*Y) \
    #     + 41.8227 * sp.sin(3*X) * sp.sin(4*Y) \
    #     - 7.7979 * sp.sin(3*X) * sp.sin(6*Y) \
    #     + 9.9258 * sp.cos(6*X)\
    #     - 19.0681 * sp.cos(6*X) * sp.cos(Y) \
    #     + 600 * sp.cos(6*X) * sp.cos(2*Y) \
    #     - 500 * sp.sin(3*Y)

    return sp_symPES


def np_symPES(X, Y):

    # np_symPES = 2*Y + Y**2 + X**2*(0.4*X**2+Y)

    np_symPES = (1/2) * (X*Y**2 - Y*X**2 - 2*X + 2*Y) \
        + (1/30) * (X**4 + Y**4)

    #
    # np_symPES = 44730.4129 \
    # - 66786.5363 * np.cos(Y) \
    # + 26352.6908 * np.cos(2*Y) \
    # - 3117.3613 * np.cos(4*Y) \
    # - 659.3217 * np.cos(6*Y) \
    # + 621.9640 * np.sin(3*X) * np.sin(Y) \
    # - 138.3050 * np.sin(3*X) * np.sin(2*Y) \
    # - 111.5488 * np.cos(8*Y) \
    # + 41.8227 * np.sin(3*X) * np.sin(4*Y) \
    # - 7.7979 * np.sin(3*X) * np.sin(6*Y) \
    # + 9.9258 * np.cos(6*X)\
    # - 19.0681 * np.cos(6*X) * np.cos(Y) \
    # + 600 * np.cos(6*X) * np.cos(2*Y) \
    # - 500 * np.sin(3*Y)

    return np_symPES


def sp_asymPES(X, Y):
    # np_asymPES = np_symPES(X, Y) + 500*X*Y
    np_asymPES = (1/2) * (X*Y**2 - Y*X**2 - 1.5*X + 2*Y) \
        + (1/30) * (X**4 + Y**4)

    return np_asymPES


def np_asymPES(X, Y):
    # np_asymPES = np_symPES(X, Y) + 500*X*Y
    np_asymPES = (1/2) * (X*Y**2 - Y*X**2 - 1.5*X + 2*Y) \
        + (1/30) * (X**4 + Y**4)

    return np_asymPES


def stationaryPts(E, x, y):
    # Calculate gradient
    # G[0] = Gx, G[1] = Gy
    G = [sp.diff(E, x), sp.diff(E, y)]

    # Find stationary points; i.e. G = 0
    # find a valid solver (solve, solveset, nonlinsolve)
    # ref.
    # 1. sympy library: https://docs.sympy.org/latest/modules/solvers/solveset.html

    # roots = sp.solveset([G[0], G[1]], [x, y], domain=sp.S.Reals)
    # roots = sp.solve([G[0], G[1]], [x, y], manual=True, set=True)
    # force=True , manual=True, set=True)
    # roots = nonlinsolve([sp.Eq(G[0], 0), sp.Eq(G[1], 0)], [x, y])

    # 2. online solver: https://www.wolframalpha.com/input/?i=solve%7B%281%2F2%29%28y%5E2-2xy-2%29%2B%284%2F30%29%28x%5E3%29%3D%3D0%2C%281%2F2%29%282xy-x%5E2%2B2%29%2B%284%2F30%29%28y%5E3%29%3D%3D0%7D

    # symPES
    roots = [
        [-0.849174, 0.849174],
        [0.789279, -0.789279],
        [-1.14427, -3.07012],
        [3.07012, 1.14427]
    ]
    # asymPES
    # mu = 1.75
    # roots = [
    #     [-0.879504, 0.765168],
    #     [0.884098, -0.648009],
    #     [-1.19982, -3.11036],
    #     [2.98781, 1.09984]
    # ]
    # mu = 1.5
    # roots = [
    #     [-0.916175, 0.678958],
    #     [-1.25225, -3.14745],
    #     [0.987403, -0.501976],
    #     [2.89612, 1.04955]
    # ]
    # mu = 1.0
    # roots = [
    #     [-1.0080658794313, 0.50500000616268894],
    #     [-1.34948, -3.2140],
    #     [1.24222, -0.183244],
    #     [2.66755, 0.919982]
    # ]

    # print(roots)
    # exit()

    nroots = len(roots)
    # for i in range(nroots):
    #     sp.pprint(roots[i])

    # Calculate Hessian
    Hessian = sp.Matrix(2, 2,   [sp.diff(G[0], x), sp.diff(G[0], y),
                                 sp.diff(G[0], y), sp.diff(G[1], y)
                                 ]
                        )
    # sp.pprint(Hessian)

    # Get the sign of Hessian eigenvalue to define minimum/maximum
    rootcol = [None]*nroots
    Hvec = [None]*nroots
    for i in range(nroots):
        Htest = Hessian.subs([(x, roots[i][0]), (y, roots[i][1])])
        # sp.pprint(Htest)

        # check the smallest eigenvalue for a 2D Hessian
        if Htest.eigenvects()[0][0] > Htest.eigenvects()[1][0]:
            index = 1
        else:
            index = 0

        # print(Htest.eigenvects()[index][0])
        if Htest.eigenvects()[index][0] > 0:
            rootcol[i] = 'bo'  # minimum
        else:
            rootcol[i] = 'ro'  # saddle point

        # store the selected eigenvector for each point
        Hvec[i] = [Htest.eigenvects()[index][2][0][0],
                   Htest.eigenvects()[index][2][0][1]]

    return G, Hvec, roots, rootcol


def GD(firstPts, E, G, Hvec):
    # ref: https://zhuanlan.zhihu.com/p/84702727

    ds = 0.1

    path = [[]]  # declear a 2D list
    path[0] = [firstPts[0], firstPts[1],
               E.subs([(x, firstPts[0]), (y, firstPts[1])])
               ]

    # move a little bit along the eigenvector direction
    secondPts = [
        path[0][0] - Hvec[0] * ds,
        path[0][1] - Hvec[1] * ds
    ]
    path.append(
        [secondPts[0], secondPts[1],
         E.subs([(x, secondPts[0]), (y, secondPts[1])])
         ]
    )
    for i in range(2, 1000):
        # calculate gradient
        Gx = G[0].subs([(x, path[i-1][0]), (y, path[i-1][1])])
        Gy = G[1].subs([(x, path[i-1][0]), (y, path[i-1][1])])
        # normalize gradient
        Gnorm = sp.sqrt(Gx**2 + Gy**2)
        Gx = Gx/Gnorm
        Gy = Gy/Gnorm
        nPts = [
            path[i-1][0] - Gx * ds,
            path[i-1][1] - Gy * ds
        ]
        path.append(
            [nPts[0], nPts[1],
             E.subs([(x, nPts[0]), (y, nPts[1])])
             ]
        )
        # print(Gnorm)
        if Gnorm <= 10e-2:
            break

    return path


def plotPES(np_func):

    # 1.
    # x = np.linspace(-2.5, 2.5, 100)
    # y = np.linspace(-3.0, 0.5, 100)

    # 2.
    x = np.linspace(-2.5, 3.2, 100)
    y = np.linspace(-3.3, 2.0, 100)

    # 3.
    # x = np.linspace(0.0, 1.0, 100)
    # y = np.linspace(-1.2, -0.2, 100)

    X, Y = np.meshgrid(x, y)
    E = np_func(X, Y)

    # aspect ratio of x and y is equal
    plt.gca().set_aspect('equal')
    symlevel = [-2.0, -0.8, 0, 0.8, 1.2]
    asymlevel = [-2.5, -1.5, -0.5, 0, 0.5, 1.0]
    fig = plt.contour(X, Y, E, 100, colors='black',
                      linestyles='-', levels=symlevel)
    # label the contour lines
    plt.clabel(fig, inline=True, fmt='%1.1f', fontsize=8)


def plotPts(roots, rootcol):
    # plot stationary points
    for i in range(len(roots)):
        plt.plot(roots[i][0], roots[i][1], rootcol[i])


def plotPath(path):
    # plot gradient descent path
    path_x = []
    path_y = []
    path_e = []
    for i in range(len(path)):
        path_x.append(path[i][0])
        path_y.append(path[i][1])
        path_e.append(path[i][2])
    plt.plot(path_x, path_y, 'r-', lw=2)


if __name__ == "__main__":

    # Define an analytical PES
    x, y = sp.symbols('x y')
    E = sp_symPES(x, y)

    # Find stationary points
    G, Hvec, roots, rootcol = stationaryPts(E, x, y)
    # print(roots)
    # print(rootcol)

    # Calculate gradient descent path
    # assign the index of the first point
    path = GD(roots[0], E, G, Hvec[0])
    # path = GD(roots[1], E, G, Hvec[1])

    # Plot figure
    plotPES(np_symPES)
    plotPts(roots, rootcol)
    plotPath(path)
    plt.show()
