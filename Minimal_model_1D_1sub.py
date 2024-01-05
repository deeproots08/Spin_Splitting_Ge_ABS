import cmath
import math
import numpy as np
from scipy.sparse import csr_matrix, coo_matrix,  bmat, diags, identity, issparse, vstack, hstack
from scipy.sparse.linalg import eigsh,eigs


###defining constants##

hbar = 1.055*10**(-34);   #in J-s
ee = 1.602*10**(-19);     #in C
hbar_mev = hbar/ee*1000   # in meV-s
me = 9.11*10**(-31);      #in kg
hbm = hbar**2*1000/2/me/ee/(1e-18) #in meV (nm)**2

import Operators2D as Op

class Model:

    def __init__(self, n, nsm, nsc, ax, omega, alphat, delta=0, meff=1, mu = 0, px=False, phi=0, kx0=0):

        """
        Initializes a new instance of the Operators_2d class.
        
        Args:
            n (int): The number of horizontal lattice points of the system.
            nsm (int): The number of horizontal lattice points of the semiconductor weak link.
            nsc (int): The number of horizontal lattice points of the superconductor on both sides.
            ax (float): The horizontal (x) lattice spacing.
            omega (float): y subband confinement energy 
            alphat (float): spin orbit coupling constant
            delta (float): superconducting gap (in meV)
            meff (float): Effective mass
            mu (float): chemical potential of the system (in meV)
            px (bool): toggles periodicity on/off
            phi (float): superconduncting phase difference across the Josephson junction
            kx0 (float): value of kx, used only if there is periodicity in x direction( i.e. px==True)
            
        """    

        self.n = n
        self.nsm = nsm
        self.nsc = nsc
        self.ax = ax
        self.omega = omega
        self.alphat = alphat
        self.delta = delta
        self.meff = meff
        self.mu = mu
        self.px = px
        self.phi = phi
        self.kx0 = kx0

        self.iden = Op.Operators_2d(self.n, 1, self.ax, 1, self.px, False).iden()
        self.kx = Op.Operators_2d(self.n, 1, self.ax, 1, self.px, False).Kx(self.kx0)
        self.kxsq = Op.Operators_2d(self.n, 1, self.ax, 1, self.px, False).KxSq(self.kx0)
        self.kxcu = Op.Operators_2d(self.n, 1, self.ax, 1, self.px, False).KxCu(self.kx0)

        self.omega_minus = self.omega/2
        self.omega_plus = self.omega
        self.kappa = math.sqrt(self.omega*self.meff/2/hbm)
        self.total = math.sqrt(self.omega_minus**2+(self.alphat*self.kappa)**2)
        self.cos = self.omega_minus/self.total
        self.sin = self.alphat*self.kappa/self.total
        self.mu_lattice = self.mu*self.iden

    def HA(self):

        return hbm/self.meff*self.kxsq + self.omega_plus*self.iden - self.mu_lattice -self.total*self.iden - (self.alphat*self.sin)**2*self.kxsq/2/self.total

    def HB(self):

        return -self.alphat*self.kx*self.cos + self.alphat**3*self.sin**2*self.cos*self.kxcu/2/self.total**2
    
    def Delta_matrix(self):


        del_nsc1 = np.full(self.nsc, self.delta)
        del_nsc2 = np.full(self.nsc, self.delta*cmath.exp(1j*self.phi))
        del_nsm = np.full(self.nsm, 0 )

        del_array = np.concatenate((del_nsc1, del_nsm, del_nsc2))

        delta_lattice = diags(del_array)

        return delta_lattice
    
    def H_p(self):

        ## particle part of Hamiltonian

        H_particle = bmat([
                [self.HA() + self.HB(), None],
                [None, self.HA()-self.HB()]
            ])
        
        return H_particle

    def H_sys(self):

        H = bmat([
                [self.HA() + self.HB(), None, None, self.Delta_matrix()],
                [None, self.HA()-self.HB(), -self.Delta_matrix(), None],
                [None,-coo_matrix.conjugate(self.Delta_matrix()), -self.HA()+self.HB(), None],
                [coo_matrix.conjugate(self.Delta_matrix()), None, None, -self.HA()-self.HB()]
            ])

        return H
    



            