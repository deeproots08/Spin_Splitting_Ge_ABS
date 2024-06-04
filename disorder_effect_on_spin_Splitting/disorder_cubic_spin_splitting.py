
'''
from config import *
#import Minimal_Model_1D as mm1D
import Minimal_Model_2D as mm2D
#import Operators2D as Op
import time
'''

import sys
import os

# Add the parent directory to sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Now you can import config
from config import *
import Minimal_Model_2D as mm2D
import time


### checking spin splitting for different magnitudes of disorder LSOC###

NSM = 60          # number
M = 20             # 
AX = 5             # lattice spacing in nm
AY = 5             # lattice spacing in nm
DELTA = 0.165      # in meV
MEFF_SC = 0.023     # proximitized InAs
MEFF_SM = 0.023     # InAs
NSC = 500
NTOT = 2*NSC + NSM
PHI = 0.5*np.pi 
ALPHAT = 400       #meV-nm^3 # 20 is within perturbative regime for widths atleast up to 100 nm
ALPHAT_REF = 1e-5
K_NUM = 4          # No of eigenvalues needed to be found 

LY = 100           # Width in nm
OMEGA = 2*4*hbm/LY**2/MEFF_SC

MU_SC =  0.8 + OMEGA/2   #in between two subbands
MU_SM =  0.8 + OMEGA/2   #in between two subbands

#### First find wavefunction of SNS with 0 spin orbit coupling

t1 = time.time()

print(t1)

Fcu_ref = mm2D.twoD_Minimal_Model(NTOT, NSM, NSC, M, AX, AY, OMEGA, ALPHAT_REF, DELTA, MEFF_SC, MEFF_SM, MU_SC, MU_SM)
            
#Hcu = Fcu.H_sys(PHI, so = 'cu', uncorr_var=1, lmbda=1)
Hcu = Fcu_ref.H_sys(PHI, so = 'cu')
Ecu, V = eigsh(Hcu ,k=K_NUM ,sigma=0, return_eigenvectors=True)

print(Ecu, V.shape)

WV_ref = V[:, K_NUM//2]

print(np.sum(np.abs(WV_ref)**2))

t2 = time.time()

print(t2-t1)

### Finding equivalent variance fo correlated disorder based on ref Wavefunction of SNS 

### Disorder parameters 

A = 1.3 ### in meV

DL = 15 ### in nm

UNCORR_VAR =  Fcu_ref.variance_of_equivalent_uncorrelated_disorder(WV_ref, A, DL)

t3 = time.time()

print('uncorr_Var = ', UNCORR_VAR)

print(t3-t2)

print('success')