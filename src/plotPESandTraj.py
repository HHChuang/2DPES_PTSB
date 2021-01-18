#!/usr/bin/env python3
#############################################################
#                                                           #
# Input :                                                   #
#       $1 = PES; e.g. NCH2_2DPES.15.dat                    #
#       $2 = coord. of important points                     #
#           format:                                         #
#           $(name of point) $(x) $(y) $(z)                 #
#           R1  x1  y1  z1                                  #
#           TS1 x2  y2  z2                                  #
#           (other important pts...)                        #
#                                                           #
#           *energy reference = energy of reactant = z1*    #
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
# 2020/05/21, Grace, debug of the size of figures.          #
# 2021/01/04, Grace, rewrite twoDwPts(); divide it into     #
#               different small functions.                  #
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
RP1dir = 'RP1.coord'
RP2dir = 'RP2.coord'
PATH_RP1 = workDir + '/' + RP1dir + '/'
PATH_RP2 = workDir + '/' + RP2dir + '/'

PATH_fig = workDir
# title = '61*61=3721 pts'
title = ''
# xaxis = "Modified IRC of TSS1 ($\sqrt{amu}\cdot$Bohr)"
xaxis = ''
# yaxis = "IRC of TSS2 ($\sqrt{amu}\cdot$Bohr)"
yaxis = ''
zaxis = 'kcal/mol'
# plt.rcParams.update({'font.size': 14})
axis_fontsize = 15

# sysIndex = 5  # 1 = NCH1, 2 = NCH2, 3 = NCH3, 4 = NCH4 and 5 = NCH5


def main():
    # 1. Input
    X, Y, E, pts_name, pts_coord = getInput()

    # 2. plot 2D PES with important points
    # plt = twoDwPts(X, Y, E, pts_name, pts_coord)
    # VRIregion(X, Y, E, pts_name, pts_coord)

    # 3. plot 2D PES with mapping trajectories
    # twoDwTraj(X, Y, E, PATH_RP1, '-r', RP1dir)  # red trajectories
    # twoDwTraj(X, Y, E, PATH_RP2, '-b', RP2dir)  # blue trajectories

    # 4. plot 3D PES with important points
    threeDwPts(X, Y, E, pts_name, pts_coord)

    # 5. plot 3D PES with mapping trajectories
    # threeDwTraj(X, Y, E, PATH_RP1, 'r')
    # threeDwTraj(X, Y, E, PATH_RP2, 'b')

    plt.show()
    plt.close()


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
    # FIXME: change square to recentagular
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


def twoD_Hartree(X, Y, E):
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    ct = plt.contour(X, Y, E,  20, colors='k')
    plt.clabel(ct, inline=1, fmt='%1.0f', fontsize=8)  # value of contour
    return plt


def twoD_kcal(X, Y, E, E_R):
    relE = E
    dim = int(math.sqrt(np.size(E)))
    for i in range(dim):
        for j in range(dim):
            relE[i][j] = (E[i][j] - E_R) * 627.5095
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    ct = plt.contour(X, Y, relE, 20, colors='k')
    plt.clabel(ct, inline=1, fmt='%1.0f', fontsize=8)  # value of contour
    return plt


def twoDwPts(X, Y, E, pts_name, pts_coord):

    ### Step 1. Select different energy unit TODO: ###
    #   Case 1.  energy in Hartree
    plt = twoD_Hartree(X, Y, E)
    #   Case 2. energy in kcal/mol
    # plt = twoD_kcal(X, Y, E, pts_coord[0][2])

    # General setting of figure 'plt'
    plt.gca().set_aspect('equal')  # aspect ratio of x and y is equal
    plt.title(title)
    plt.xlabel(xaxis, fontsize=axis_fontsize)
    plt.ylabel(yaxis, fontsize=axis_fontsize)
    plt.xticks(fontsize=axis_fontsize)
    plt.yticks(fontsize=axis_fontsize)

    ### Step 2. plot important point on this 2D-PES ###
    shift = 1.5  # shift x and y value by 1.5 times
    npts = totline(str(sys.argv[2]))
    for n in range(npts):
        plt.plot(pts_coord[n][0], pts_coord[n][1], 'ro')
        plt.text(shift*pts_coord[n][0], shift*pts_coord[n][1], pts_name[n], weight='bold', backgroundcolor='white',
                 verticalalignment='top', multialignment='right', fontsize=axis_fontsize)

    # FIXME:
    # label the VRI region within a red box
    # VRIheight = [3, 6, 3, 3, 6]
    # # m = sysIndex
    # plt.plot([pts_coord[1][0], pts_coord[1][0]],
    #          [-VRIheight[sysIndex-1],  VRIheight[sysIndex-1]], '-r')
    # plt.plot([pts_coord[1][0], pts_coord[2][0]], [
    #          VRIheight[sysIndex-1],  VRIheight[sysIndex-1]], '-r')
    # plt.plot([pts_coord[2][0], pts_coord[2][0]],
    #          [-VRIheight[sysIndex-1],  VRIheight[sysIndex-1]], '-r')
    # plt.plot([pts_coord[1][0], pts_coord[2][0]],
    #          [-VRIheight[sysIndex-1], -VRIheight[sysIndex-1]], '-r')

    # # fig.set_size_inches(5, 7.5)
    # fig.tight_layout()
    # fig = plt.show()
    # plt.close(fig)
    return plt


def VRIregion(X, Y, E, pts_name, pts_coord):
    # Calculate relative energy in kcal/mol, and use reactant
    # energy as the energy reference
    E_R = pts_coord[0][2]
    dim = int(math.sqrt(np.size(E)))
    relE = E
    for i in range(dim):
        for j in range(dim):
            relE[i][j] = (E[i][j] - E_R) * 627.5095

    # print(E_R)
    fig = plt.figure()
    ax = plt.gca()
    ax.set_aspect('equal')  # aspect ratio of x and y is equal

    # set x and y limit
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xlim((pts_coord[1][0], pts_coord[2][0]))

    yrange = ([-3, 3], [-6, 6], [-3, 3], [-3, 3], [-6, 6])
    ax.set_ylim((yrange[sysIndex-1][0], yrange[sysIndex-1][1]))

    vmin = 50
    vmax = 120
    vint = 8
    plt.contour(X, Y, relE, levels=np.linspace(vmin, vmax, vint), colors='k')
    ct = plt.contourf(X, Y, relE, levels=np.linspace(
        vmin, vmax, vint), cmap="RdBu_r")

    plt.title(title)
    plt.xlabel(xaxis, fontsize=axis_fontsize)
    plt.ylabel(yaxis, fontsize=axis_fontsize)
    plt.xticks(fontsize=axis_fontsize)
    plt.yticks(fontsize=axis_fontsize)

    # only plot TSS1 and TSS2; TSS1 = 1, TSS2 = 2 in the pts list.
    # values for TSS1/TSS2 in NCHn systems, n=1-5.
    shift_x = [[0.2, -1.1], [0.3, -2.1],
               [0.2, -1.1], [0.2, -1.05], [0.3, -2.1]]
    shift_y = [[0.6, 0.6], [1.1, 1.1], [0.6, 0.6], [0.6, 0.6], [1.1, 1.1]]
    for n in [1, 2]:
        print(n, pts_name[n], shift_x[sysIndex-1]
              [n-1], shift_y[sysIndex-1][n-1])
        plt.plot(pts_coord[n][0], pts_coord[n][1], 'ro')
        plt.text(pts_coord[n][0] + shift_x[sysIndex-1][n-1], pts_coord[n][1] + shift_y[sysIndex-1][n-1],
                 pts_name[n], weight='bold', backgroundcolor='white', verticalalignment='top', multialignment='right', fontsize=axis_fontsize)
    # order of savefig() and show() is important
    # fig.savefig(PATH_fig + '2DwPts.png', dpi=100)
    # fig.set_size_inches(5, 7.5)
    fig.tight_layout()
    # plt.colorbar()  # TODO:
    plt.show()
    plt.close(fig)


def twoDwTraj(X, Y, E, PATH, color, dirname):
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    fig = plt.figure()
    plt.gca().set_aspect('equal')  # aspect ratio of x and y is equal

    ct = plt.contour(X, Y, E, 20, colors='k')
    plt.title(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    plt.xticks(fontsize=axis_fontsize)
    plt.yticks(fontsize=axis_fontsize)
    # remove the indices on axis
    # plt.xticks([])
    # plt.yticks([])

    # Plot trajectories
    numRP, nameTraj = totTraj(PATH)

    for i in range(numRP):
        traj_X, traj_Y, traj_E = getTraj(PATH, nameTraj[i])
        plt.plot(traj_X, traj_Y, color)

    # fig.set_size_inches(5, 7.5)
    # remove the blank part
    fig.tight_layout()
    # plt.subplots_adjust(hspace=0,wspace=0)
    fig.savefig(PATH_fig + dirname + '.png', dpi=100)
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


def threeD_Hartree(X, Y, E, ptE):
    ax = plt.gca(projection='3d')
    surf = ax.plot_surface(X, Y, E,  cmap=mpl.cm.binary,
                           alpha=0.6, linewidth=1, edgecolors='black')
    return ax, ptE


def threeD_kcal(X, Y, E, refE, pts_coord):
    npts = totline(str(sys.argv[2]))
    dim2 = np.size(E)
    dim = int(math.sqrt(dim2))
    relE = E
    # change energy into kcal/mol for PES
    for i in range(dim):
        for j in range(dim):
            relE[i][j] = (E[i][j] - refE) * 627.5095
    # change energy into kcal/mol for selected points
    nPts = np.size(pts_coord, 0)
    ptE = pts_coord[:, 2]
    for i in range(nPts):
        ptE[i] = (pts_coord[i][2]-refE)*627.5095

    ax = plt.gca(projection='3d')
    surf = ax.plot_surface(X, Y, relE,  cmap=mpl.cm.binary,
                           alpha=0.6, linewidth=1, edgecolors='black')
    return ax, ptE


def threeDwPts(X, Y, E, pts_name, pts_coord):
    # Step 1. Set up the unit of energy
    #   Case 1. Hartree
    ax, ptE = threeD_Hartree(X, Y, E, pts_coord[:, 2])
    #   Case 2. kcal/mol
    # ax, ptE = threeD_kcal(X, Y, E, pts_coord[0][2], pts_coord)
    #   General setting
    ax.pbaspect = [2.0, 0.6, 0.25]
    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)
    ax.set_zlabel(zaxis)
    # range of z-axis; range of energy profile
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    # ax.set_zlim(-20, 100)

    # Archive setting for test ###
    # plt.title(title)
    # plt.xticks([0.47, 0.97, 1.47])
    # plt.yticks([98, 108, 118])
    # ax = fig.add_subplot((111), aspect='equal', projection='3d')

    # ax.set_aspect('equal')  # ratio of x and y axes, 'auto', 'equal'

    # # ratio of 3D box
    # ax.pbaspect = [1.0, 2.0, 1.0]

    # set_axes_equal(ax)
    # set_aspect_equal_3d(ax)
    # ax.auto_scale_xyz([-16, 5], [-16, 16], [0, 100])
    # # ct = plt.contour(X, Y, E, 25, colors='k')
    ### End of archive setting ###

    # Step 2. Plot important points
    npts = totline(str(sys.argv[2]))
    shift = 1.0

    for n in range(npts):
        ax.scatter(pts_coord[n][0], pts_coord[n][1],
                   ptE[n], marker='o', c='r', edgecolors='k', zorder=10)
        ax.text(shift*pts_coord[n][0], shift*pts_coord[n][1], shift*ptE[n],
                pts_name[n], weight='bold', backgroundcolor='white', zorder=10, fontsize=12)


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
