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
ALPHAT = 400  #meV-nm^3 # 20 is within perturbative regime for widths atleast up to 100 nm

'''
MU0 = 0.001
#W0 = 20 
THETA = np.pi/5
W = np.linspace(20, 100, 51)
OMEGA = 2*4*hbm/W**2/MEFF_SC
MU = OMEGA + MU0 + math.tan(THETA)*OMEGA

MU_SC = MU
MU_SM = MU
'''
W = np.linspace(20, 100, 2)

OMEGA = 2*4*hbm/W**2/MEFF_SC

MU_ex = np.linspace(0, OMEGA[0]*2, 2)

MU =    OMEGA + MU_ex

OMEGA_G, MU_G = np.meshgrid(OMEGA, MU) 


K = 8

Ep_Split_cu = np.empty((OMEGA.shape[0],MU.shape[0]))

root = Path.cwd()
data_folder1 = 'result_dataset' 
folder_loc1 = os.path.join(root, data_folder1 )
os.makedirs(folder_loc1, exist_ok=True)

data_folder2 = f'cubic_spin_Split_map_soc_{ALPHAT}_nsm_{NSM}'
folder_loc2 = os.path.join( folder_loc1, data_folder2 )
os.makedirs(folder_loc2, exist_ok=True)

for i in range(OMEGA.shape[0]):

    for j in range(MU.shape[0]):

        Fcu = M.oneD_Minimal_Model(N, NSM, NSC, MY, AX, OMEGA[i], ALPHAT, DELTA, MEFF_SC, MEFF_SM, MU[j], MU[j])
        
        Hcu = Fcu.H_sys(PHI, so = 'cu')
        
        Ecu = eigsh(Hcu ,k=K,sigma=0, return_eigenvectors=False)
        
        Ecu = np.sort(Ecu)
        
        Ep_Split_cu[i,j] = (Ecu[K//2+1] - Ecu[K//2])/DELTA

        print(i,j )

        if i%1 == 0 and j%1 == 0:

            data_folder3 = f'data_till_omega_{OMEGA[i]:.1f}_mu_{MU[j]:.1f}'
            folder_loc3 = os.path.join(folder_loc2, data_folder3 )
            os.makedirs(folder_loc3, exist_ok=True)

            varfile_path = os.path.join(folder_loc3, 'variables.pkl')
            # Save all variables to a file using joblib
            joblib.dump(Ep_Split_cu, varfile_path)

print(f"Variables saved to '{varfile_path}'")

plt.rc('font', size=18)

plt.figure(figsize=(10, 8))

plt.contourf(OMEGA_G, MU_G, Ep_Split_cu)

# Increase the font size of the axes labels and title
plt.xlabel(r'$\Omega$(in meV)', fontsize=24)
plt.ylabel(r'$\mu$ (in meV)', fontsize=24)
plt.title(r'Cubic Spin splitting map', fontsize=24)

#plt.grid(True)
plt.grid(True)
#plt.legend(['E1', 'E2', 'E3', 'E4', 'E5'])

# Set x-axis limits to -1 and +1
#plt.xlim(-1, 1)

plt.tight_layout()

plt.colorbar()

# Increase the number of ticks on the x-axis
#plt.xticks(np.linspace(-1, 1, num=10))

# Increase the number of ticks on the y-axis
#plt.yticks(np.linspace(np.min(Ep2), np.max(Ep2), num=10))


# Adjust subplot parameters to remove white space
#plt.subplots_adjust(left=0.08, bottom=0.08, right=0.95, top=0.92)

# Save the plot as an image file
plt.savefig(os.path.join(folder_loc2, 'cubic_spin_split_map.png'))


plt.show()
