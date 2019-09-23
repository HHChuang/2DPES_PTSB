#!/usr/bin/env python3

# analytical formula:

# 2019/08/16, Grace


import matplotlib.pyplot as plt
import numpy as np


def main():
    xlist = np.linspace(0.0, 1.0, 100)
    ylist = np.linspace(-1.2, -0.2, 100)
    X, Y = np.meshgrid(xlist, ylist)

    # Symmetric PES
    Z = symPES(X, Y)
    pt_X = [0.52, 0.52, 0.52, 0.10, 0.95]
    pt_Y = [-0.22, -0.85, -1.0, -1.1, -1.1]
    plotPES(X, Y, Z, pt_X, pt_Y)

    # Asymmetric PES
    Z = symPES(X, Y) + + 500*X*Y
    plotPES(X, Y, Z, pt_X, pt_Y)


def symPES(X, Y):
    symPES = 44730.4129-66786.5363*np.cos(Y)+26352.6908*np.cos(2*Y) \
        - 3117.3613*np.cos(4*Y)-659.3217*np.cos(6*Y) \
        + 621.9640*np.sin(3*X)*np.sin(Y)-138.3050*np.sin(3*X)*np.sin(2*Y) \
        - 111.5488*np.cos(8*Y) + 41.8227*np.sin(3*X)*np.sin(4*Y) \
        - 7.7979*np.sin(3*X)*np.sin(6*Y) + 9.9258*np.cos(6*X)\
        - 19.0681*np.cos(6*X)*np.cos(Y)+600*np.cos(6*X) * np.cos(2*Y) \
        - 500*np.sin(3*Y)

    return symPES


def plotPES(X, Y, Z, pt_X, pt_Y):
    plt.figure()
    plt.gca().set_aspect('equal')  # aspect ratio of x and y is equal
    cp = plt.contour(X, Y, Z, 50, colors='black', linestyles='-')
    # for n in range(5):
    #     plt.plot(pt_X[n], pt_Y[n], 'ro')

    plt.show()


main()
