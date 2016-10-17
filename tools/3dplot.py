from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(-10, 10, 1)
Y = np.arange(-10, 10, 1)
Z = np.array([200, 162, 128, 98, 72, 50, 32, 18, 8, 2, 2, 8, 18, 32, 50, 72, 98, 128,162,200])
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)


plt.show()
