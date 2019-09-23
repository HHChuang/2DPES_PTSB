import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.xlabel('x')
plt.ylabel('y')
# plt.clabel('E(x,y)')

delta = 0.02
x = np.arange(-2.5, 2.5, delta)
y = np.arange(-3, 0.5, delta)
x, y = np.meshgrid(x, y)
# model PES
# E = 2 * y + y**2 + (y + 0.4 * x**2) * x**2
E = 0.5 * (x * y**2 - y * x**2 - 2 * x + 2 * y) + (1 / 30) * (x**4 + y**4)

# customize the z axis
ax.set_zlim(-5, 5)

surf = ax.plot_surface(x, y, E, cmap=cm.jet, alpha=0.3,  # alpha for the transparency
                       linewidth=0.1, antialiased=False, shade=True)
ax.contour(x, y, E, 50, lw=1, linestyles="dashed", offset=-5)
# Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)

# add path #################################################################
prefix = 'IRC.'  # Modified this part ########################
# prefix = 'ModifiedIRC.'
fileex = '.dat'  # file extension
for i in range(1, 51):
    name = prefix + str(i) + fileex
    rawdata0 = [r.split() for r in open(name).read().strip().split(' ')]
    rawdata0 = [[float(c) for c in r.split()]
                for r in open(name).read().strip().split('\n')]
    x_path0 = []
    y_path0 = []
    E_path0 = []
    for j in range(0, 1000):
        x_path0.append(rawdata0[j][0])
        y_path0.append(rawdata0[j][1])
        E_path0.append(rawdata0[j][2])
    path0 = ax.plot(x_path0, y_path0, E_path0, '-', color='red', linewidth=2)
plt.show()
