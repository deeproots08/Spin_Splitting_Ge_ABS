import numpy as np
#from numpy.linalg import eigvals, eigvalsh
#from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh,eigs
import matplotlib.pyplot as plt
import math, cmath
import sys
import os
import joblib
from pathlib import Path
import pickle
import time
from mpl_toolkits.mplot3d import Axes3D

import Minimal_Model_1D as M

###defining constants##

hbar = 1.055*10**(-34);   #in J-s
ee = 1.602*10**(-19);     #in C
hbar_mev = hbar/ee*1000   # in meV-s
me = 9.11*10**(-31);      #in kg
hbm = hbar**2*1000/2/me/ee/(1e-18) #in meV (nm)**2

### Plotting spin splitting while varying width(omega) and chemical potential ###

NSM = 300          # in nm
MY = 4             # no of subbands
AX = 1             # lattice spacing in nm
DELTA = 0.165      # in meV
MEFF_SC = 0.05     # proximitized Ge 
MEFF_SM = 0.05     # Ge
NSC = 2000
N = 2*NSC + NSM
PHI = 0.5*np.pi 
ALPHAT = 20  #meV-nm^3 # 20 is within perturbative regime for widths atleast up to 100 nm

W = np.linspace(20, 100, 3)

OMEGA = 2*4*hbm/W**2/MEFF_SC

OMEGA = OMEGA[::-1]

MU_ex = np.linspace(0.001, OMEGA[-1]*2, 3)

MU =    OMEGA + MU_ex

OMEGA_G, MU_G = np.meshgrid(OMEGA, MU) 


K = 8

Ep_Split_cu = np.empty((OMEGA.shape[0],MU.shape[0]))

root = Path.cwd()
data_folder1 = 'result_dataset' 
folder_loc1 = os.path.join(root, data_folder1 )
os.makedirs(folder_loc1, exist_ok=True)
data_folder2 = f'Linear_spin_Split_map_soc_{ALPHAT}_nsm_{NSM}'
folder_loc2 = os.path.join( folder_loc1, data_folder2 )
os.makedirs(folder_loc2, exist_ok=True)


for i in range(OMEGA.shape[0]):

    for j in range(MU.shape[0]):

        t1 = time.time()

        Fcu = M.oneD_Minimal_Model(N, NSM, NSC, MY, AX, OMEGA[i], ALPHAT, DELTA, MEFF_SC, MEFF_SM, MU[j], MU[j])
        
        Hcu = Fcu.H_sys(PHI, so = 'li')
        
        Ecu = eigsh(Hcu ,k=K,sigma=0, return_eigenvectors=False)
        
        Ecu = np.sort(Ecu)
        
        Ep_Split_cu[i,j] = (Ecu[K//2+1] - Ecu[K//2])/DELTA

        t2 = time.time()

        print(i,j, '\t', round(t2-t1) )

        if i%1 == 0 and j%1 == 0:

            data_folder3 = f'data_till_omega_{OMEGA[i]:.1f}_mu_{MU[j]:.1f}'
            folder_loc3 = os.path.join(folder_loc2, data_folder3 )
            os.makedirs(folder_loc3, exist_ok=True)

            varfile_path = os.path.join(folder_loc3, 'variables.pkl')
            # Save all variables to a file using joblib
            joblib.dump(Ep_Split_cu, varfile_path)

print(f"Variables saved to '{varfile_path}'")

# Set fontsize for the plot
plt.rc('font', size=16)          

# Create a figure and a 3D axis
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')

# Create a 3D surface plot
ax.plot_surface(OMEGA_G, MU_G, Ep_Split_cu, cmap='jet')

# Label the axes
ax.set_xlabel(r'$\Omega$(in meV)', fontsize=24)
ax.set_ylabel(r'$\mu$ (in meV)', fontsize=24)
ax.set_zlabel(r'$\dfrac{\delta E}{\Delta}$', fontsize=24)
ax.set_title(r'Linear Spin splitting map', fontsize=24)

# Add grid
plt.grid(True)

# Save the plot as an image file
plt.savefig(os.path.join(folder_loc2, 'linear_spin_split_map.png'))

# Show the plot
plt.show()