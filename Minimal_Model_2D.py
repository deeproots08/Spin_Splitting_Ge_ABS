from config import *

import Operators2D as Op


class twoD_Minimal_Model:

    def __init__(self, n, nsm, nsc, m, ax, ay, omega=0, alphat=0, delta=0, meff_sc=1, meff_sm = 1 , mu_sc = 0, mu_sm=0, dis_mean=0, dis_var=1, px=False, py=False,  epsilon=0):

        """
        Initializes a new instance of the oneD_Minimal_Model class.
        
        Args:
            n (int): The number of horizontal lattice points of the system.
            nsm (int): The number of horizontal lattice points of the semiconductor weak link.
            nsc (int): The number of horizontal lattice points of the superconductor on both sides.
            m (int): The number of vertical lattive points
            ax (float): The horizontal (x) lattice spacing.
            ay (float): The vertical (y) lattice spacing.
            meff (float): Effective mass 
            alphat (float): spin orbit coupling constant
            mu_s (float): chemical potential of the superconducting region of SNS
            mu_sm (float): chemical potential of the semiconductingconducting region of SNS
            px (bool, optional): Flag indicating periodicity in x direction. Defaults to False.
            epsilon (float): Kx^4 prefactor to remove higher k states from dipping

        Notes:
            - The ordering of lattice points for kx, kxsq, kxcu(imported for operator class file) is from left to right, bottom to top, similar to x and y axis in 2D g

        """

        self.n = n
        self.nsm = nsm
        self.nsc = nsc
        self.m = m
        self.ax = ax
        self.ay = ay
        self.omega = omega
        self.meff_sc = meff_sc
        self.meff_sm = meff_sm
        self.alphat = alphat
        self.delta = delta
        self.mu_sc = mu_sc
        self.mu_sm  = mu_sm
        self.dis_mean = dis_mean
        self.dis_var = dis_var
        self.px = px
        self.py = py
        self.epsilon = epsilon
        self.twodOP = Op.Operators_2d(self.n, self.m, self.ax, self.ay, self.px, self.py)

    def ysq_matrix(self):

        row = np.empty(0)
        col = np.empty(0)
        data = np.empty(0)

        shift = (((self.m/2)+1)*self.ay)

        for j in range(self.m):

            for i in range(self.n):
            
                row = np.append(row, j*self.n + i)
                col = np.append(col, j*self.n + i)
                data = np.append(data, ((j+1)*self.ay-shift)**2)

        M = csr_matrix((data,(row,col)),shape=(self.m*self.n, self.m*self.n))

        return M

    def chem_pot(self):

        """
        Returns a sparse 2D array diagonal array of chemical potential of size (nmxnm) defining the chemical potential landscape in two dimension.
        """

        nsc_data = np.full(self.nsc, self.mu_sc)
        nsm_data = np.full(self.nsm, self.mu_sm)
        remaining_data = np.full(self.n - self.nsc - self.nsm, self.mu_sc)

        data_one_row = np.concatenate((nsc_data, nsm_data, remaining_data))

        data_lattice = np.tile(data_one_row, self.m)

        M = diags(data_lattice)

        return M

    def dis_pot_matrix(self, dis_mag=0):

        """
        Generate a sparse diagonal matrix of size mnXmn 
        gives onsite disorder potential sampled from gaussian for the semiconducting region only. 
        onsite disorder potential for SC region is taken as zero

        Returns:
        dis_pot_matrix (csr_matrix): Sparse diagonal matrix mnXmn

        """
        final_array = []

        for _ in range(self.m):
            # Generate random values for the non-zero entries
            values = [random.gauss(self.dis_mean, self.dis_var) for _ in range(self.nsm)]
            values = np.array(values)

            # Create the 1D array for this iteration
            array = np.concatenate((np.zeros(self.nsc), values, np.zeros(self.nsc)))

            # Append the array to the final array
            final_array = np.append(final_array, array)

        # Convert final_array to a NumPy array
        final_array = np.array(final_array, dtype=float)

        # Multiply by magnitude
        final_array *= dis_mag

        # Create a sparse diagonal matrix from the final array
        dis_pot_matrix = diags(final_array, 0, format='csr')

        return dis_pot_matrix

    def Disorder_full(self, uncorr_var, lmbda=0):

        """
        returns  mnXmn sparse matrix with diagonal terms giving disorder potential

        The disorder on each site is sampled from a normal distribution with variance uncorr_var and mean zero 

        this def assumes disorder is through out the josephson junction including superconductors unlike dis_pot_mat def
        
        Args: 
            uncorr_var (float): variance of the standard normal distribution from which disorder is sampled
            lmbda (float): to modify the magnitude of disorder  
        
        """

        gaussian_array = lmbda* np.random.normal(loc=0, scale=math.sqrt(uncorr_var), size=self.m*self.n)

        dis_array = diags(gaussian_array)

        return dis_array

    def spin_orbit(self, s, so='cu', kx=0, ky=0):

        """
        returns  <s|H_SO|j> i.e a matris in xy lattice basis i.e. a nmXnm matrix
        
        Args:
            s, sp (int){0 or 1}: refers to spin index up and down resp. for spin orbit, sp is opposite of s for spin orbit to be non zero
            so : spin orbit coupling, 'cu' if cubic rashba, 'li' if linear rashba
        
        """

        if so == 'cu':
            #print('cu')
            
            if s==0:

                return self.twodOP.H_cubic_so(self.alphat, kx, ky)
            
            else:

                return self.twodOP.H_cubic_so(self.alphat, kx, ky).transpose().conjugate()

        elif so == 'li':
            #print('li')
            # Linear SOC #
            ### Note in later definitions, soc coupling is evoked only when s!=s' ###
            ### Hence not coding that statement here ###
            if s ==0 : 

                return self.twodOP.H_linear_so(self.alphat, kx, ky)
            
            else:

                return self.twodOP.H_linear_so(self.alphat, kx, ky).transpose().conjugate()
                   
        else:

            return self.twodOP.zer()

    def H_same_spin(self, dis_mag=0, uncorr_var=1, lmbda=0, kx=0, ky=0):

        return (hbm/self.meff_sm* self.twodOP.KxSq(kx) + hbm/self.meff_sm* self.twodOP.KySq(ky) - self.chem_pot() + 1/4/hbm*self.meff_sm*self.omega**2*self.twodOP.ysq() +self.dis_pot_matrix(dis_mag) +self.Disorder_full(uncorr_var, lmbda))
 
    def H_s_sp(self, s, sp, so='cu', dis_mag=0, kx=0, ky=0):

        """
        returns  <s|H|sp>
        
        Args: 
            s, sp (int){0 or 1}: refers to spin index up and down resp.
            so (str): spin orbit coupling, cu if cubic rashba, li if linear rashba
        
        """

        if s==sp:

            return self.H_same_spin(dis_mag, kx, ky)
        
        else:

            return self.spin_orbit(s, so, kx, ky)

    def H_particle(self, so = 'cu', dis_mag=0, kx=0, ky=0):

        """
        returns  Particle space Hamiltonian
        
        """

        V = csr_matrix((0,2*self.m*self.n))

        for s in range(2):

            H = csr_matrix((self.n*self.m,0))

            for sp in range(2):

                H = hstack( ( H, self.H_s_sp(s, sp, so, dis_mag, kx, ky) ) )

            V = vstack((V,H))

        return V

    def Del_array(self, phi=0):

        """
        Returns a sparse 2D array of delta of size (nmxnm) defining the superducting gap landscape in two dimension.
        """

        nsc_data = np.full(self.nsc, self.delta)
        nsm_data = np.full(self.nsm, 0)
        remaining_data = np.full(self.n - self.nsc - self.nsm, self.delta*cmath.exp(1j*phi))

        data = np.concatenate((nsc_data, nsm_data, remaining_data))

        data_full_lattice = np.tile(data, self.m)

        M = diags(data_full_lattice)

        return M                

    def Delta_s_sp(self,s,sp, phi=0):

        """
        returns  <s|Delta|sp> mnXmn sparse matrix
        
        Args: 
            s, sp (int){0 or 1}: refers to spin index up and down resp.
        
        """
        if s!=sp:

            if s==0:

                return self.Del_array(phi)
            
            else:

                return -self.Del_array(phi)
            
        else:

            return self.twodOP.zer()

    def Delta_mat(self, phi=0):

        """
        returns <particle| Delta | hole> 
        
        """

        V = csr_matrix((0,2*self.m*self.n))

        for s in range(2):

            H = csr_matrix((self.n*self.m,0))

            for sp in range(2):
                #print(self.H_is_jsp(i,j,s,sp))    
                H = hstack( ( H, self.Delta_s_sp(s,sp,phi) ) )
                
            V = vstack((V,H))

        return V
    
    def H_sys(self,phi=0, so = 'cu', dis_mag=0, kx=0, ky=0):

        ### -kx, -ky of hoel space matrix only matters in case of periodic 

        M = bmat([
           [self.H_particle( so, dis_mag, kx, ky), self.Delta_mat(phi)],
           [-self.Delta_mat(phi).conjugate(), -self.H_particle( so, dis_mag, -kx, -ky).conjugate()]
       ]) 

        return M