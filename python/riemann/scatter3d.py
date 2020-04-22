'''
==============
3D scatterplot
==============

Demonstration of a basic scatterplot in 3D.
'''

# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
import numpy as np

# Fixing random state for reproducibility
np.random.seed(19680801)


def randrange(n, vmin, vmax):
    '''
    Helper function to make an array of random numbers having shape (n, )
    with each number distributed Uniform(vmin, vmax).
    '''
    return (vmax - vmin)*np.random.rand(n) + vmin

fig = plt.figure()
#ax = fig.add_subplot(211, projection='3d')

n = 100

# For each set of style and range settings, plot n random points in the box
# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
i = 1
for m, zlow, zhigh in [('o', -50, -25), ('^', -30, -5)]:
    xs = randrange(n, 0, 10)
    ys = randrange(n, 0, 100)
    print(type(xs),xs.shape)
    #zs = randrange(n, zlow, zhigh)
    zs = (0.2*xs + 0.1*ys)*(zhigh-zlow)/12 + zlow
    ax = fig.add_subplot(2,1,i, projection='3d')
    ax.scatter(xs, ys, zs, marker=m)
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    i = i +1


plt.show()
