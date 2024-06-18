import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import sys
import os

# Add the parent directory to sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SAMPLE_NUM = 500

Decay_length = 10

# Load noise surface from the array 
noise_surface_V0VN_avg = np.load('avg_disorder_correlation_matrix_over_{}_realizations.npy'.format(SAMPLE_NUM))

G = np.load('exponential_correlation_matrix_decaylength_{}.npy'.format(Decay_length))


P  = noise_surface_V0VN_avg - G

#P  = noise_surface_V0VN_avg

#P  = G

print(P)

# plot the avg correlation
#fig, ax = plt.subplots()
plt.imshow(P)
plt.colorbar()
plt.show()
