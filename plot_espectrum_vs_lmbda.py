import sys
import os

# Add the parent directory to sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config import *
#import Minimal_Model_1D as mm1D
#import Minimal_Model_2D as mm2D
#import Operators2D as Op
#import time


'''
LMBDA = np.array([0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1,1])
LMBDA = np.linspace(0, 0.2, 10)
LMBDA = np.linspace(0.075, 0.2, 51)
LMBDA = np.linspace(0.15, 0.21, 81)

Ep_cu = np.load('Ep_vs_lmbda_{:.2f}_to_{:.2f}.npy'.format(LMBDA[0],LMBDA[-1]))

print(Ep_cu.shape)

plt.plot(LMBDA, Ep_cu)

plt.show()
'''

LMBDA = np.array([0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1,1])
LMBDA = np.linspace(0, 0.2, 10)
LMBDA = np.linspace(0.075, 0.2, 51)
LMBDA = np.linspace(0.15, 0.21, 81)
LMBDA = np.linspace(0.15, 0.21, 31)

Ep_li = np.load('Ep_vs_lmbda_{:.2f}_to_{:.2f}_li.npy'.format(LMBDA[0],LMBDA[-1]))

print(Ep_li.shape)

plt.plot(LMBDA, Ep_li)

plt.show()