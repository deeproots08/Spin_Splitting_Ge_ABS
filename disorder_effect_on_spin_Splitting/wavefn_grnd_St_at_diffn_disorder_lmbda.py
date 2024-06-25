
import sys
import os

# Add the parent directory to sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config import *
#import Minimal_Model_1D as mm1D
import Minimal_Model_2D as mm2D
#import Operators2D as Op
import time



NSM = 60            # number of sm x lattice points
M = 60              # should be atleast thrice confinement length
AX = 5              # lattice spacing in nm
AY = 5              # lattice spacing in nm
DELTA = 0.165       # in meV
MEFF_SC = 0.023     # proximitized InAs
MEFF_SM = 0.023     # InAs
NSC = 800
NTOT = 2*NSC + NSM
PHI = 0.5*np.pi 
ALPHAT = 400        # meV-nm^3 # 20 is within perturbative regime for widths atleast up to 100 nm
K_NUM = 8           # No of eigenvalues needed to be found 
#W_MIN = 50

LY = 100           # Width in nm
OMEGA = 2*4*hbm/LY**2/MEFF_SC


MU_SC =  0.8 + OMEGA/2   # in between two subbands
MU_SM =  0.8 + OMEGA/2   # in between two subbands

UNCORR_VAR = 1167  ### uncorrelated variance of disorder

LMBDA = 0.21
SEED = 42  ## seed for the disorder 

t1 = time.time()

Fcu = mm2D.twoD_Minimal_Model(NTOT, NSM, NSC, M, AX, AY, OMEGA, ALPHAT, DELTA, MEFF_SC, MEFF_SM, MU_SC, MU_SM)
            
#Hcu = Fcu.H_sys(PHI, so = 'cu', uncorr_var=1, lmbda=1)
Hcu = Fcu.H_sys(PHI, so = 'cu', uncorr_var=UNCORR_VAR, lmbda=LMBDA, seed=SEED)

Ecu, Vcu = eigsh(Hcu ,k=K_NUM ,sigma=0, return_eigenvectors=True)

# Get the indices that would sort Ecu
sorted_indices = np.argsort(Ecu)

# Sort Ecu using the sorted indices
sorted_Ecu = Ecu[sorted_indices]

# Sort the columns of V according to the sorted indices
sorted_V = Vcu[:, sorted_indices]

print(sorted_Ecu, sorted_V.shape)

WV_ref_spinu = sorted_V[:, K_NUM//2]

WV_ref_spind = sorted_V[:, K_NUM//2+1]

WV_ex_spinu = sorted_V[:, K_NUM//2+2]

WV_ex_spind = sorted_V[:, K_NUM//2+3]

# Save the array
np.save('GS_spinU_wv_disorder_lmbda_{}.npy'.format(LMBDA), WV_ref_spinu)

np.save('GS_spinD_wv_disorder_lmbda_{}.npy'.format(LMBDA), WV_ref_spind)

np.save('FE_spinU_wv_disorder_lmbda_{}.npy'.format(LMBDA), WV_ex_spinu)

np.save('FE_spinD_wv_disorder_lmbda_{}.npy'.format(LMBDA), WV_ex_spind)

#Ep_cu[i,:] = Ecu/DELTA

t2 = time.time()

print( "time: ", round(t2-t1), "lmbda: ", LMBDA)

print('GS spin up energy: ', sorted_Ecu[K_NUM//2]/DELTA)

print('GS spin down energy: ', sorted_Ecu[K_NUM//2+1]/DELTA)

print('EX spin up energy: ', sorted_Ecu[K_NUM//2+2]/DELTA)

print('EX spin down energy: ', sorted_Ecu[K_NUM//2+3]/DELTA)

