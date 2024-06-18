
import math
import numpy as np
import gstools as gs # gstools for correlated noise models
import matplotlib.pyplot as plt
import random

# min and max x and y coordinates of domain
x_min = 0
x_max = 500

y_min = 0
y_max = 500

# spacing between lattice points
delta_x = delta_y = 1

# Make sure random seeds are the same for all iterations with the same lattice
seed = 11111
random.seed(seed)
np.random.seed(seed)

npoints_x = int( np.ceil( (x_max-x_min) / delta_x + 1 ) )
npoints_y = int( np.ceil( (y_max-y_min) / delta_y + 1 ) )

xs = np.linspace(x_min, x_max, npoints_x)
ys = np.linspace(y_min, y_max, npoints_y)

# Gaussian covariance model
# covariance is of the form exp(-r**2 / len_scale**2)
class Gau(gs.CovModel):

    # h = r/len_scale
    def cor(self, h):
        return np.exp(-(h)**2)

# Define the length scale
len_scale = 10

# define the noise variance
var = 2

# model with given variance and covariance
model = Gau(dim=2, var=var, len_scale=len_scale)

# create a model seed
s = random.randint(0, 1000000)

# generate random noise surface
noise_surface_model = gs.SRF(model, seed=s)
noise_surface = noise_surface_model.structured([xs, ys])

print(noise_surface.shape)

# plot the noise surface
fig, ax = plt.subplots()
ax.imshow(noise_surface)
plt.show()