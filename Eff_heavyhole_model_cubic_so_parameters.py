import cmath
import numpy as np
from scipy.sparse import csr_matrix, coo_matrix,  bmat, diags, identity, issparse, vstack, hstack
from scipy.sparse.linalg import eigsh,eigs


###defining universal constants##

hbar = 1.055*10**(-34);   #in J-s
ee = 1.602*10**(-19);     #in C
hbar_mev = hbar/ee*1000   # in meV-s
me = 9.11*10**(-31);      #in kg
hbm = hbar**2*1000/2/me/ee/(1e-18) #in meV (nm)**2

### defining germanium constants - Winkler SOC book ###

gamma1 = 13.38
gamma2 = 4.24
gamma3 = 5.69
#kappa = 3.41
#q = 0.06

'''
Note: The effective heavy hole model has been derived assuming quantum well potential i z direction and  treating the Ez term in z direction 
as  perturbative. Electric field values and z direction terms, if any, should be taken in accordance to the above.
'''

def Em(g1 = gamma1, g2 = gamma2, g3 = gamma3, m=1, Lz=1, E_strain=0):

    '''
    returns Energy of mth heavy hole state on which we are projecting

    Arg:    
        g1, g2, g3 : material constants, defaults to Germanium constants
        m : z subband quantum no. , defaults to 1. 
        m refers to z subband subspace on which Hamiltonian is projected, for effective Heavy hole, we are projecting on one z subband only
        Lz : z confinement length in nm
        E_strain : strain energy in meV. Defaults to 0

    '''

    return (g1-2*g2)*hbm*math.pi**2*m**2/Lz**2 - E_strain

def El_hh(g1 = gamma1, g2 = gamma2, g3 = gamma3, l, Lz=1, E_strain=0):

    '''
    returns Energy of l_th heavy hole state

    Arg:    
        g1, g2, g3 : material constants, defaults to Germanium constants
        l : z subband quantum no.
        Lz : z confinement length in nm
        E_strain : strain energy in meV

    '''

    return (g1-2*g2)*hbm*math.pi**2*l**2/Lz**2 - E_strain

def El_lh(g1 = gamma1, g2 = gamma2, g3 = gamma3, l, Lz=1, E_strain=0):

    '''
    returns Energy of l_th light hole state 

    Arg:    
        g1, g2, g3 : material constants, defaults to Germanium constants
        l : z subband quantum no. 
        Lz : z confinement length in nm
        E_strain : strain energy in meV

    '''

    return (g1+2*g2)*hbm*math.pi**2*l**2/Lz**2 + E_strain