## file to make exponential matrix grid

import math
import numpy as np
import matplotlib.pyplot as plt

def exponential_correlation(xg,yg,xg0=0,yg0=0,dl=1,var=1):

    XG0 = xg0*np.ones(xg.shape)

    YG0 = yg0*np.ones(yg.shape)

    DX = xg - XG0

    DY = yg - YG0

    L = DX**2 + DY**2

    return var*np.exp(-L/dl**2)



# min and max x and y coordinates of domain
x_min = 0
x_max = 200

y_min = 0
y_max = 200

# spacing between lattice points
delta_x = delta_y = 1

npoints_x = int( np.ceil( (x_max-x_min) / delta_x + 1 ) )
npoints_y = int( np.ceil( (y_max-y_min) / delta_y + 1 ) )

xs = np.linspace(x_min, x_max, npoints_x)
ys = np.linspace(y_min, y_max, npoints_y)

X, Y = np.meshgrid(xs,ys)

Decay_length = 10

var= 2

G = exponential_correlation(X, Y, xs[100], ys[100], Decay_length,var)

np.save('exponential_correlation_matrix_decaylength_{}'.format(Decay_length),G)

# plot the avg correlation
fig, ax = plt.subplots()
ax.imshow(G)
plt.show()


