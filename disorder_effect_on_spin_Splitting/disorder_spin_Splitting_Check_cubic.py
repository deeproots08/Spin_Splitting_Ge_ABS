
import sys
import os

# Add the parent directory to sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config import *
#import Minimal_Model_1D as mm1D
import Minimal_Model_2D as mm2D
#import Operators2D as Op
import time

### checking spin splitting for different magnitudes of disorder LSOC###

NSM = 60          # number
M = 60             # should be atleast thrice confinement length
AX = 5             # lattice spacing in nm
AY = 5             # lattice spacing in nm
DELTA = 0.165      # in meV
MEFF_SC = 0.023     # proximitized InAs
MEFF_SM = 0.023     # InAs
NSC = 800
NTOT = 2*NSC + NSM
PHI = 0.5*np.pi 
ALPHAT = 400       #meV-nm^3 # 20 is within perturbative regime for widths atleast up to 100 nm
K_NUM = 8          # No of eigenvalues needed to be found 
#W_MIN = 50

LY = 100           # Width in nm
OMEGA = 2*4*hbm/LY**2/MEFF_SC

MU_SC =  0.8 + OMEGA/2   #in between two subbands
MU_SM =  0.8 + OMEGA/2   #in between two subbands

UNCORR_VAR = 1167  ### uncorrelated variance of disorder

LMBDA = 1E-3  ## scaling of disorder

#SEED = 42  ## seed for the disorder
 
### Sample 100 times for each value of uncorr_Var

SAMPLE_NUM = 3

Ep_Split_cu = np.empty((SAMPLE_NUM,1))

for i in range(SAMPLE_NUM):

    #print(time.time())

    Fcu = mm2D.twoD_Minimal_Model(NTOT, NSM, NSC, M, AX, AY, OMEGA, ALPHAT, DELTA, MEFF_SC, MEFF_SM, MU_SC, MU_SM)
                
    #Hcu = Fcu.H_sys(PHI, so = 'cu', uncorr_var=1, lmbda=1)
    Hcu = Fcu.H_sys(PHI, so = 'cu', uncorr_var=UNCORR_VAR, lmbda=LMBDA)

    Ecu = eigsh(Hcu ,k=K_NUM ,sigma=0, return_eigenvectors=False)

    Ecu = np.sort(Ecu)

    Ep_Split_cu[i,0] = (Ecu[K_NUM//2+1] - Ecu[K_NUM//2])/DELTA

    #print(time.time())

    print("sample no. {}".format(i))


# Calculating mean
splitting_mean_value = Ep_Split_cu.mean()

# Calculating variance
splitting_variance_value = Ep_Split_cu.var()

print("uncorrvar: {}".format(UNCORR_VAR))

print("lmbda: {}".format(LMBDA))

print("Mean:", splitting_mean_value)
print("Variance:", splitting_variance_value)

print("spin splitting array for cubic: {}".format(Ep_Split_cu))

print("for spin splitting array w/o disorder refer to end line of of file cubic_spin_splitting_wdisorder.txt")

'''
Fcu = mm2D.twoD_Minimal_Model(NTOT, NSM, NSC, M, AX, AY, OMEGA, ALPHAT, DELTA, MEFF_SC, MEFF_SM, MU_SC, MU_SM)
            
#Hcu = Fcu.H_sys(PHI, so = 'cu', uncorr_var=1, lmbda=1)
Hcu = Fcu.H_sys(PHI, so = 'cu', uncorr_var=1, lmbda=0)

Ecu = eigsh(Hcu ,k=K_NUM ,sigma=0, return_eigenvectors=False)

Ecu = np.sort(Ecu)

Ep_Split_cu_nd = (Ecu[K_NUM//2+1] - Ecu[K_NUM//2])/DELTA

print("spin splitting array for cubic w/o disorder: {}".format(Ep_Split_cu_nd))

'''