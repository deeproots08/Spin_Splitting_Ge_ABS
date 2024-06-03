    def Correlation_2D(self, A, i, ip, j, jp, dl=1):

        """
            Returns correlation function between point labelled by i, j on a n X m lattice
            
            Args:
                A (float) : normalizing ocnstant for correlation function
                i, j (int) : index of two poins between which correlation function needs to be calculated
                dl (float) : decay length of correlation function 
            
        """

        #x_dis = ((i%n) - (j%n))*ax   ## x distance between two points on lattice with index i,j labelled in graph way 

        #y_dis = ((i//n) - (j//n))*ay   ## y distance between two points on lattice with index i,j labelled in graph way 

        x_dis = (i -ip)*self.ax   ## x distance between two points on lattice with index i,j labelled in graph way 

        y_dis = (j -jp)*self.ay   ## y distance between two points on lattice with index i,j labelled in graph way 

        dis_sq = x_dis**2 + y_dis**2

        lmbda_sq = dl**2*(self.ax**2+self.ay**2)  #### check if you can do without ax ay, my suspicion is you can ###

        return A**2 * math.exp(-dis_sq/lmbda_sq)
       
    def variance_of_equivalent_uncorrelated_disorder(self, WV, A, dl=1):

        """
            Returns equivalent <VV> or variance for uncorrelated disorder in 2D given correlated disorder magnitude A   
            
            Args:
                WV (array) : Array of dim (nXm, 1) representing wavefunction 
                A (float) : normalizing ocnstant for correlation function
                dl (float) : decay length of correlation function 
            
        """
    
        WVb = WV.reshape((4,self.m,self.n))           ### reshaping one D wavevector for ease of summing with the correlation function which relies on i,j being latice index

        var_corr = np.empty((4*self.m*self.n, 4*self.m*self.n))                        ### intializing variance of uncorrelated disorder

        for k in range(4):

            for kp in range(4):

                for i in range(self.n):

                    for j in range(self.m):

                        for ip in range(self.n):

                            for jp in range(self.m):

                                ### summing square of given wavefunction with correlation function to find corr_var

                                var_corr[i+j*self.n+k*self.m*self.n, ip+jp*self.n+kp*self.m*self.n] = abs(WVb[k,j,i])**2*abs(WVb[kp,jp,ip])**2*self.Correlation_2D(A, i, ip, j, jp, dl)
                                

        ### Please note, take the ax**2*ay**2 (area of single element of integral), 
        ### include it either in normalization of wavefunction, or in variance reimann sum calculation. 
        ### Do not include it at both the places.

        var_corr_sum = np.sum(var_corr)*self.ax**2*self.ay**2

        var_uncorr = np.empty((self.m*self.n, 4, 4))

        for k in range(4):

                for kp in range(4):

                    for j in range(self.m):

                        for i in range(self.n):

                            ### summing square of given wavefunction with correlation function to find corr_var

                            var_uncorr[i+j*self.n, k, kp] = abs(WVb[k,j,i])**2*abs(WVb[kp,j,i])**2

        factor = np.sum(var_uncorr)*self.ax*self.ay

        eqv_uncorr_sum = var_corr_sum/ factor   

        ### returns variance of equivalent uncorrelated disorder
        ### to be used in disorder def 

        return eqv_uncorr_sum
