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

## Adress of file to load ###
root = Path.cwd()
data_folder = f'cubic_spin_Split_map_soc_{ALPHAT}_nsm_{NSM}'
folder_loc = os.path.join(root, f'data_till_omega_{OMEGA[i]:.1f}_mu_{MU[j]:.1f}',data_folder )
os.makedirs(folder_loc, exist_ok=True)
varfile_path = os.path.join(folder_loc, 'variables.pkl')  
Ep_Split_cu =  joblib.load(varfile_path)                       ## load variables
