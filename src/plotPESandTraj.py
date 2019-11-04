#!/usr/bin/env python3
#############################################################
#                                                           #
# Input :                                                   #
#       $1 = PES; e.g. NCH2_2DPES.15.dat                    #
#       $2 = coord. of important points                     #
#           format:                                         #
#           $(name of point) $(x) $(y) $(z)                 #
#                                                           #
# History:                                                  #
# 2018/10/23, Grace, revised; add input arguments           #
# 2018/12/22, Grace, add mapping traj. on 2D and also       #
#   classify traj. into 2 groups, RP1 and RP2.              #
# 2019/07/11, Grace, change the unit of reaction coord from #
#   the amount of grid points to mass-weighted coordinate.  #
#   And also modify it from square figure to recentangular. #
# 2019/09/24, Grace, add VRIregion(); zoom-in the region    #
#   between TSS1 and TSS2, and also colorcode this area.    #
#############################################################

import sys
import os
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D
# from itertools import cycle

workDir = os.getcwd()
PATH_RP1 = workDir + '/RP1.coord/'
PATH_RP2 = workDir + '/RP2.coord/'

PATH_fig = workDir
# title = '61*61=3721 pts'
title = ''
xaxis = "IRC of TSS1 ($\sqrt{amu}\cdot$Bohr)"
yaxis = "IRC of TSS2 ($\sqrt{amu}\cdot$Bohr)"
zaxis = 'kcal/mol'
plt.rcParams.update({'font.size': 14})


def main():
    # 1. Input
    X, Y, E, pts_name, pts_coord = getInput()

    # 2. plot 2D PES with important points
    twoDwPts(X, Y, E, pts_name, pts_coord)
    # VRIregion(X, Y, E, pts_name, pts_coord)

    # 3. plot 2D PES with mapping trajectories
    # twoDwTraj(X, Y, E, PATH_RP1, 'r-')  # red trajectories
    # twoDwTraj(X, Y, E, PATH_RP2, 'b-')  # blue trajectories

    # 4. plot 3D PES with important points
    # threeDwPts(X, Y, E, pts_name, pts_coord)

    # 5. plot 3D PES with mapping trajectories
    # threeDwTraj(X, Y, E, PATH_RP1, 'r')
    # threeDwTraj(X, Y, E, PATH_RP2, 'b')


def totline(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    f.close()
    return i + 1


def totTraj(PATH):
    path, dirs, files = next(os.walk(PATH))
    name = np.array(files)
    return len(files), name


def getInput():
    # import PES w/ x, y and z coord
    # TODO: change square to recentagular
    tline = totline(str(sys.argv[1]))
    dim = int(math.sqrt(tline))
    x, y, e = np.loadtxt(str(sys.argv[1]), delimiter=' ', unpack=True)
    X = np.reshape(x, (dim, dim))
    Y = np.reshape(y, (dim, dim))
    E = np.reshape(e, (dim, dim))

    # import name and coord of pts
    # FIXME: pts_coord may have datatype problem
    pts_name = list()
    pts_coord = list()
    with open(str(sys.argv[2])) as f:
        for i, line in enumerate(f):
            pts_name.append(line.split()[0])
            pts_coord.append(line.split()[1:4])
    f.close()
    pts_name = np.array(pts_name)
    pts_coord = np.array(pts_coord)
    pts_coord = pts_coord.astype(np.float)
    return X, Y, E, pts_name, pts_coord


def getTraj(path, fname):
    x = []
    y = []
    e = []

    with open(path + str(fname)) as f:
        for i, line in enumerate(f):
            x.append(line.split()[0])
            y.append(line.split()[1])
            e.append(line.split()[2])
    f.close()

    x = np.array(x)
    y = np.array(y)
    e = np.array(e)
    x = x.astype(np.float)
    # y = y.astype(np.int)
    y = y.astype(np.float)
    e = e.astype(np.float)
    return x, y, e


def twoDwPts(X, Y, E, pts_name, pts_coord):
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    fig = plt.figure()
    plt.gca().set_aspect('equal')  # aspect ratio of x and y is equal

    # Calculate relative energy in kcal/mol, and use reactant
    # energy as the energy reference
    E_R = pts_coord[0][2]
    dim = int(math.sqrt(np.size(E)))
    relE = E
    for i in range(dim):
        for j in range(dim):
            relE[i][j] = (E[i][j] - E_R) * 627.5095

    ct = plt.contour(X, Y, relE, 20, colors='k')
    plt.clabel(ct, inline=1, fmt='%1.0f', fontsize=8)  # value of contour
    plt.title(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    # plt.xticks([0.47, 0.97, 1.47])
    # plt.yticks([98, 108, 118])

    # plot important points
    npts = totline(str(sys.argv[2]))
    # NCH1
    shift_x = [0, -5, 1, 0, 0]
    shift_y = [3, 1, 1, 3, -1]
    # H3CO
    # shift_x = [0.0, -0.3, -0.3, -0.1, -0.1]
    # shift_y = [0.3, 0.3, 0.3, 0.3, -0.2]
    for n in range(npts):
        plt.plot(pts_coord[n][0], pts_coord[n][1], 'ro')
        plt.text(pts_coord[n][0] + shift_x[n], pts_coord[n][1]+shift_y[n],
                 pts_name[n], weight='bold', backgroundcolor='white',
                 verticalalignment='top', multialignment='right', fontsize=10)
    # order of savefig() and show() is important
    # fig.savefig(PATH_fig + '2DwPts.png', dpi=100)

    # label the VRI region within a red box
    VRIheight = [3,6,3,3,6]
    m = 5 # 1 = NCH1, 2 = NCH2, 3 = NCH3, 4 = NCH4 and 5 = NCH5
    plt.plot( [ pts_coord[1][0], pts_coord[1][0] ], [ -VRIheight[m-1],  VRIheight[m-1] ], '-r')
    plt.plot( [ pts_coord[1][0], pts_coord[2][0] ], [  VRIheight[m-1],  VRIheight[m-1] ], '-r')
    plt.plot( [ pts_coord[2][0], pts_coord[2][0] ], [ -VRIheight[m-1],  VRIheight[m-1] ], '-r')
    plt.plot( [ pts_coord[1][0], pts_coord[2][0] ], [ -VRIheight[m-1], -VRIheight[m-1] ], '-r')

    fig = plt.show()
    plt.close(fig)


def VRIregion(X, Y, E, pts_name, pts_coord):
    # Calculate relative energy in kcal/mol, and use reactant
    # energy as the energy reference
    E_R = pts_coord[0][2]
    dim = int(math.sqrt(np.size(E)))
    relE = E
    for i in range(dim):
        for j in range(dim):
            relE[i][j] = (E[i][j] - E_R) * 627.5095

    fig = plt.figure()
    ax = plt.gca()  # .set_aspect('equal')
    ax.set_aspect('equal')  # aspect ratio of x and y is equal
    # set x and y limit
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xlim((pts_coord[1][0], pts_coord[2][0]))
    ax.set_ylim((-3, 3))

    vmin = 50
    vmax = 120
    vint = 8
    plt.contour(X, Y, relE, levels=np.linspace(vmin, vmax, vint), colors='k')
    ct = plt.contourf(X, Y, relE, levels=np.linspace(
        vmin, vmax, vint), cmap="RdBu_r")

    plt.title(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)

    # plot important points
    npts = totline(str(sys.argv[2]))
    shift_x = [0.0, 0.1, -0.9, 0.0, 0.0]
    shift_y = [0.1, 0.5, 0.5, 0.1, -0.1]
    for n in range(npts):
        plt.plot(pts_coord[n][0], pts_coord[n][1], 'ro')
        plt.text(pts_coord[n][0]+shift_x[n], pts_coord[n][1]+shift_y[n],
                 pts_name[n], weight='bold', backgroundcolor='white', verticalalignment='top', multialignment='right', fontsize=12)
    # order of savefig() and show() is important
    # fig.savefig(PATH_fig + '2DwPts.png', dpi=100)
    plt.colorbar()
    plt.show()
    plt.close(fig)


def twoDwTraj(X, Y, E, PATH, color):
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    fig = plt.figure()
    plt.gca().set_aspect('equal')  # aspect ratio of x and y is equal

    ct = plt.contour(X, Y, E, 20, colors='k')
    plt.title(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)

    # Plot trajectories
    numRP, nameTraj = totTraj(PATH)

    for i in range(numRP):
        traj_X, traj_Y, traj_E = getTraj(PATH, nameTraj[i])
        plt.plot(traj_X, traj_Y, color)

    # fig.savefig(PATH_fig + '2DwTraj.png', dpi=100)
    fig = plt.show()
    plt.close(fig)

# https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
# https://github.com/matplotlib/matplotlib/issues/1077/


def set_axes_radius(ax, origin, radius):
    ax.set_xlim3d([origin[0] - radius, origin[0] + radius])
    ax.set_ylim3d([origin[1] - radius, origin[1] + radius])
    ax.set_zlim3d([origin[2] - radius, origin[2] + radius])


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])

    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    set_axes_radius(ax, origin, radius)


def set_aspect_equal_3d(ax):
    """Fix equal aspect bug for 3D plots."""

    xlim = ax.get_xlim3d()
    ylim = ax.get_ylim3d()
    zlim = ax.get_zlim3d()

    from numpy import mean
    xmean = mean(xlim)
    ymean = mean(ylim)
    zmean = mean(zlim)

    plot_radius = max([abs(lim - mean_)
                       for lims, mean_ in ((xlim, xmean),
                                           (ylim, ymean),
                                           (zlim, zmean))
                       for lim in lims])

    ax.set_xlim3d([xmean - plot_radius, xmean + plot_radius])
    ax.set_ylim3d([ymean - plot_radius, ymean + plot_radius])
    ax.set_zlim3d([zmean - plot_radius, zmean + plot_radius])
#########################################################


def threeDwPts(X, Y, E, pts_name, pts_coord):
    # change hartree to relative energy
    npts = totline(str(sys.argv[2]))
    #   energy reference : reactant
    Renergy = pts_coord[0][2]
    dim2 = np.size(E)
    dim = int(math.sqrt(dim2))
    for i in range(dim):
        for j in range(dim):
            E[i][j] = (E[i][j] - Renergy) * 627.5095

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.pbaspect = [2.0, 0.6, 0.25]
    # plt.title(title)
    # plt.xticks([0.47, 0.97, 1.47])
    # plt.yticks([98, 108, 118])
    # ax = fig.add_subplot((111), aspect='equal', projection='3d')

    ax.set_aspect('equal')  # ratio of x and y axes, 'auto', 'equal'

    # ratio of 3D box
    ax.pbaspect = [1.0, 2.0, 1.0]

    # set_axes_equal(ax)
    # set_aspect_equal_3d(ax)
    # ax.auto_scale_xyz([-16, 5], [-16, 16], [0, 100])

    surf = ax.plot_surface(X, Y, E,  cmap=mpl.cm.binary, alpha=0.6,
                           linewidth=1, edgecolors='black')

    # ct = plt.contour(X, Y, E, 25, colors='k')

    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)
    ax.set_zlabel(zaxis)
    # range of z-axis; range of energy profile
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    ax.set_zlim(-20, 100)

    # plot important points
    npts = totline(str(sys.argv[2]))
    # H3CO
    shift_x = [-0.5, -0.6, -0.6, -0.8, -1.0]
    shift_y = [0.3, 0.3, 0.3, 0.3, 0.4]
    for i in range(npts):
        pts_coord[i][2] = (pts_coord[i][2] - Renergy) * 627.5095

    for n in range(npts):
        ax.scatter(pts_coord[n][0], pts_coord[n][1],
                   pts_coord[n][2], marker='o', c='r', edgecolors='k', zorder=10)
        ax.text(pts_coord[n][0] + shift_x[n], pts_coord[n][1] + shift_y[n], pts_coord[n]
                [2] + 2, pts_name[n], weight='bold', backgroundcolor='white', zorder=10, fontsize=12)
    # fig.savefig(PATH_fig + '3DwPts.png', dpi=100)
    plt.show()
    plt.close(fig)


def threeDwTraj(X, Y, E, PATH, color):
    # change hartree to relative energy
    #   energy reference : reactant
    dim2 = np.size(E)
    dim = int(math.sqrt(dim2))
    E_R = E[0][int((dim - 1) / 2)]
    # E_R = pts_coord[0][2]
    for i in range(dim):
        for j in range(dim):
            E[i][j] = (E[i][j] - E_R) * 627.5095

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, E,  cmap=mpl.cm.binary, alpha=0.6,
                           linewidth=1, edgecolors='black', antialiased=True)
    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)
    ax.set_zlabel(zaxis)
    # range of z-axis; range of energy profile
    plt.set_zlim([0, 100])

    # Plot trajectories with kcal/mol
    numRP, nameTraj = totTraj(PATH)

    for i in range(numRP):
        traj_X, traj_Y, traj_E = getTraj(PATH, nameTraj[i])
        traj_E = (traj_E - Renergy) * 657.5095
        ax.plot(traj_X, traj_Y, traj_E, '-', c=color, zorder=10)

    plt.show()


main()
