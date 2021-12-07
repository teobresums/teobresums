#!/usr/bin/python3

"""
Script to auto generate parfiles and run them

Given a template parfile with basic setup, 
generates parfiles for all combination of a given subset of parameters
(e.g. to vary binary masses and spins). Then it run the code.

SB 10/2018
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import time

from EOBUtils import *

if __name__ == "__main__": 


    # Setup -----------------------------------------
    
    # Base dir & parfile
    based = "./"
    basep = "test_HM.par"

    # Set new values/ranges for parameters (Use lists)
    N = 10
    nu = 0.25*np.random.rand(N)
    q    = (1. - 2.*nu + np.sqrt(1. - 4.*nu))/(2.*nu) 
    chi1 = 1-2*np.random.rand(N)
    chi2 = 1-2*np.random.rand(N)
    
    # Pack them into a dictionary
    # NOTE: keys must match those in parfile otherwise ignored
    n = {'q': q,
         'chi1': chi1,
         'chi2': chi2}
    
    # ------------------------------------------
    # DO NOT CHANGE BELOW HERE
    # ------------------------------------------

    # Generate parfiles ----------------------------
    
    # Read the base parfile
    d = read_parfile_dict(basep)

    # Generate combinations
    #x, keys = combine_parameters(n)

    x, keys = [], []
    for key in n.keys():
        keys.append(key)

    for i in range(N):
        tmp = (q[i],chi1[i],chi2[i])
        x.append(tmp)

    # Write parfiles
    parfile = []
    basen, ext = os.path.splitext(basep)
    for s in range(len(x)):
        for i in range(len(keys)):
            d[keys[i]] = str(x[s][i])
        # Output to file
        #print(d)
        parfile.append(based+"/"+basen+"_"+str(s)+ext)
        write_parfile_dict(parfile[-1], d)
        #print("Written {}".format(s))

    # Run  ----------------------------

    j = 0
    
    for p in parfile:
        run(p)
        os.remove(p)

        lvec = [2,3,4,5]
        mvec = [2,3,4,5]
        
        fig, axs = plt.subplots(3, 1)

        for i in range(len(lvec)):
            l = lvec[i];
            m = mvec[i];
            
            HM = np.genfromtxt("highermodes/hlm_interp_l"+str(l)+"_m"+str(m)+".txt",delimiter=" ")

            t   = HM[:,0]
            A   = HM[:,1]
            phi = HM[:,2]

            omg = np.zeros_like(phi)
            omg[1:] = np.diff(phi)/np.diff(t)

            axs[0].plot(t,abs(A),label=r'$'+str(l)+str(m)+'$')
            axs[1].plot(t,phi,label=r'$\phi_{'+str(l)+','+str(m)+'}$')
            axs[2].plot(t,omg,label=r'$\omega_{'+str(l)+','+str(m)+'}$')

        axs[0].legend()
        axs[0].set_ylabel('$ A $')
        axs[1].set_ylabel('$ \phi $')
        axs[2].set_ylabel('$ \omega $')
        axs[2].set_xlabel('$ t $')
        axs[1].text(0, 0.75*phi[-1], '(%.3g,%.3g,%.3g)' % (q[j],chi1[j],chi2[j]), fontsize = 12)
        plt.show()

        j=j+1
    


        
