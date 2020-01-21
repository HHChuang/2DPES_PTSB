#!/usr/bin/env python3
#
#   Program:
#       For presentation, combine 2D-PES, 3D-PES and molecular evolution of one dynamic trajectory into one animation.
#
#   Input:
#       $1 = PES
#       $2 = Important point on this surface
#       $3 = Coordinate of trajectory
#       $4 = xyz structure file of trajectory
#
#   Output:
#
#
# Reference of tutorial blog or youtube
#   1. display static image of plotly lib
#   https://plot.ly/python/static-image-export/
#   2. basic example of creating animation
#   https://learndataanalysis.org/a-basic-example-how-to-create-animation-matplotlib-tutorial/
#   3. combine two 2D animations in one fig
#   https://pythonmatplotlibtips.blogspot.com/2018/01/combine-two-2d-animations-in-one-figure-matplotlib-artistanimation.html
#   4. combine 3D and two 2D animations in one figure
#   https://pythonmatplotlibtips.blogspot.com/2018/01/combine-3d-two-2d-animations-in-one-figure-artistdanimation.html
#   5. convert xyz file into molelcular structure
#   https://github.com/zotko/xyz2graph
#
#   History:
#       2020/01/19, Grace
#

import sys
import os
import math
import numpy as np
from copy import deepcopy
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D
import plotly
from xyz2graph import MolGraph, to_plotly_figure
import plotly.graph_objects as go
import plotly.io as pio
# from IPython.display import Image
# from plotly.offline import offline, iplot
# To display the figure defined by this dict, use the low-level plotly.io.show function

import matplotlib.animation as animation

workDir = os.getcwd()
PATH_fig = workDir

title = ''
xaxis = ''
yaxis = ''
zaxis = 'kcal/mol'
axis_fontsize = 6

fig = plt.figure()
threeD = fig.add_subplot(1, 2, 2,projection='3d') # top and bottom right 
mol = fig.add_subplot(2,2,1)
twoD = fig.add_subplot(2, 2, 3)  # bottom left

def main():
    # 1. Import rawdata
    X, Y, E, pts_name, pts_coord = getPES()
    traj = getTraj(pts_coord) # format: x,y,E
    
    # 2. Plot the background: 2D-PES and 3D-PES
    twoDwPts(twoD, X, Y, E, pts_name, pts_coord)
    threeDwPts(threeD, X, Y, E, pts_name, pts_coord)
    # FIXME:
    molGeo(mol)
    
    # 3. Combine above subroutines and export them into one video
    # ani(fig,twoD,threeD,traj)
    plt.show()


def totline(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    f.close()
    return i + 1

def getTraj(pts_coord):
    # read the coordinate of trajectory
    tline=totline(str(sys.argv[3]))
    traj=[]
    with open(str(sys.argv[3])) as f:
        for i,line in enumerate(f):
            traj.append(line.split())
    f.close()

    traj = np.array(traj)
    traj = traj.astype(np.float)

    # Calculate relative energy in kcal/mol, and use reactant energy as the energy 
    ncol,nrow = np.shape(traj)
    E_R = pts_coord[0][2]
    relE = deepcopy(traj) #copy by values, not by reference
    for i in range(ncol-1):
        relE[i][2] = (traj[i][2] - E_R) * 627.5095

    return relE

def getPES():
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

def twoDwPts(twoD, X, Y, E, pts_name, pts_coord):
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    plt.gca().set_aspect('equal')  # aspect ratio of x and y is equal

    # Calculate relative energy in kcal/mol, and use reactant energy as the energy reference
    E_R = pts_coord[0][2]
    dim = int(math.sqrt(np.size(E)))
    relE = np.arange(dim*dim).reshape(dim,dim)
    for i in range(dim):
        for j in range(dim):
            relE[i][j] = (E[i][j] - E_R) * 627.5095

    ct = twoD.contour(X, Y, relE, 10, colors='k')
    # plt.clabel(ct, inline=1, fmt='%1.0f', fontsize=8)  # value of contour
    # plt.title(title)
    # twoD.set_xlabel(xaxis, fontsize=axis_fontsize)
    # twoD.set_ylabel(yaxis, fontsize=axis_fontsize)
    # # plt.xticks([0.47, 0.97, 1.47])
    # # plt.yticks([98, 108, 118])

    # plot important points
    npts = totline(str(sys.argv[2]))
    # NCH1
    shift_x = [0, -5, 1, 0, 0]
    shift_y = [3, 1, 1, 3, -1]
    for n in range(npts):
        twoD.plot(pts_coord[n][0], pts_coord[n][1], 'ro')
        twoD.text(pts_coord[n][0] + shift_x[n], pts_coord[n][1]+shift_y[n],
                 pts_name[n], weight='bold', backgroundcolor='white',
                 verticalalignment='top', multialignment='right', fontsize=6)
    
def threeDwPts(threeD, X, Y, E, pts_name, pts_coord):
    # change hartree to relative energy
    npts = totline(str(sys.argv[2]))
    #   energy reference : reactant
    Renergy = pts_coord[0][2]
    dim2 = np.size(E)
    dim = int(math.sqrt(dim2))
    relE = np.arange(dim*dim).reshape(dim,dim)
    for i in range(dim):
        for j in range(dim):
            relE[i][j] = (E[i][j] - Renergy) * 627.5095

    # ratio of 3D box FIXME:
    # threeD.pbaspect = [2.0, 0.6, 0.25]
    # threeD.set_aspect('equal')  # ratio of x and y axes, 'auto', 'equal'

    threeD.set_xlabel(xaxis)
    threeD.set_ylabel(yaxis)
    threeD.set_zlabel(zaxis)
    # range of z-axis; range of energy profile
    threeD.set_xlim(X.min(), X.max())
    threeD.set_ylim(Y.min(), Y.max())
    threeD.set_zlim(-20, 100) 

    


    surf = threeD.plot_surface(X, Y, relE,  rstride=1, cstride=1,cmap=mpl.cm.coolwarm, linewidth=1, antialiased=False)
    # fig.colorbar(surf, shrink=0.5, aspect=10)

    # ct = plt.contour(X, Y, E, 25, colors='k')

    # plot important points
    npts = totline(str(sys.argv[2]))
    # H3CO
    shift_x = [-0.5, -0.6, -0.6, -0.8, -1.0]
    shift_y = [0.3, 0.3, 0.3, 0.3, 0.4]
    for i in range(npts):
        pts_coord[i][2] = (pts_coord[i][2] - Renergy) * 627.5095

    for n in range(npts):
        threeD.scatter(pts_coord[n][0], pts_coord[n][1],
                   pts_coord[n][2], marker='o', c='r', edgecolors='k', zorder=10)
        threeD.text(pts_coord[n][0] + shift_x[n], pts_coord[n][1] + shift_y[n], pts_coord[n]
                [2] + 2, pts_name[n], weight='bold', backgroundcolor='white', zorder=10, fontsize=6)

def molGeo(mol):
    mg = MolGraph()
    mg.read_xyz(str(sys.argv[4]))
    # Create the Plotly figure object
    mol = to_plotly_figure(mg)
    pio.renderers.default = 'png'
    mol.show()
    # test = go.Figure(data=mol) 
    # test.show(renderer='png')

    # Plot the figure
    # Image(pio.to_image(mol,format='png'))
    # mol.to_image(format='png')
    # mol.write_image('test.png') #work

def ani(fig,twoD,threeD,traj):
    # test for printing full trajectory 
    # twoD.plot(traj[:,0],traj[:,1],'r-')
    # threeD.plot(traj[:,0],traj[:,1],traj[:,2],'r-',zorder=10)

    ncol,nrow = np.shape(traj)
    t = np.linspace(0, 80, ncol)
    x = traj[:,0]
    y = traj[:,1]
    z = traj[:,2]
    
    lines = []
    for i in range(ncol):
        head = i - 1
        head_slice = (t > t[i] - 1.0) & (t < t[i])
        #print(x[head],y[head],z[head])
        # 3-D
        line1,  = threeD.plot(x[:i], y[:i], z[:i],color='black')
        line1a, = threeD.plot(x[head_slice], y[head_slice], z[head_slice],color='red', linewidth=2)
        line1e, = threeD.plot([x[head]], [y[head]], [z[head]],color='red', marker='o', markeredgecolor='r')
        # 2-D
        line2,  = twoD.plot(x[:i], y[:i],color='black')
        line2a, = twoD.plot(x[head_slice], y[head_slice],color='red', linewidth=2)
        line2e, = twoD.plot(x[head], y[head],color='red', marker='o', markeredgecolor='r')
        lines.append([line1,line1a,line1e,line2,line2a,line2e])

    ani = animation.ArtistAnimation(fig, lines, interval=100, blit=True)

main()
