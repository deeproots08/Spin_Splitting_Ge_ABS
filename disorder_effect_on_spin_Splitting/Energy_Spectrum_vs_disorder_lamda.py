
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

#LMBDA = np.array([0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1,1])
#LMBDA = np.linspace(0.075, 0.2, 51)
LMBDA = np.linspace(0.15, 0.21, 81)
SEED = 42  ## seed for the disorder 

### Sample only once for each value of lmbda

#SAMPLE_NUM = 3

Ep_cu = np.empty((LMBDA.shape[0],K_NUM))

for i in range(LMBDA.shape[0]):

    t1 = time.time()

    Fcu = mm2D.twoD_Minimal_Model(NTOT, NSM, NSC, M, AX, AY, OMEGA, ALPHAT, DELTA, MEFF_SC, MEFF_SM, MU_SC, MU_SM)
                
    #Hcu = Fcu.H_sys(PHI, so = 'cu', uncorr_var=1, lmbda=1)
    Hcu = Fcu.H_sys(PHI, so = 'cu', uncorr_var=UNCORR_VAR, lmbda=LMBDA[i], seed=SEED)

    Ecu = eigsh(Hcu ,k=K_NUM ,sigma=0, return_eigenvectors=False)

    Ecu = np.sort(Ecu)

    Ep_cu[i,:] = Ecu/DELTA

    t2 = time.time()

    print("sample no. {}".format(i), "time: ", round(t2-t1), "lmbda: ", LMBDA[i])

    print(Ep_cu[i,:])

print(Ep_cu, LMBDA, UNCORR_VAR)

# Save the array
np.save('Ep_vs_lmbda_{:.2f}_to_{:.2f}.npy'.format(LMBDA[0],LMBDA[-1]), Ep_cu)