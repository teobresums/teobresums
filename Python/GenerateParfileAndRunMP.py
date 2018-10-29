#!/usr/bin/python3

"""
Script to auto generate parfiles and run them in parallel with multiprocess

Given a template parfile with basic setup, 
generates parfiles for all combination of a given subset of parameters
(e.g. to vary binary masses and spins). Then it runs the code with multiprocess.

SB 10/2018
"""

import multiprocessing as mp
#from functools import partial
#import numpy as np
import matplotlib.pyplot as plt
import os
import sys
              
if __name__ == "__main__": 


    # Setup -----------------------------------------
    
    # Base dir & parfile
    based = "./"
    basep = "test_NQCIteration.par"

    # Set new values/ranges for parameters (Use lists)
    q = [1., 1.2, 1.4]
    chi1 = [0., 0.8]
    chi2 = [0., 0.8, 0.9]
    
    # Pack them into a dictionary
    # NOTE: keys must match those in parfile otherwise ignored
    n = {'q': q,
         'chi1': chi1,
         'chi2': chi2}

    # Set the number of processes
    nproc = 4
    
    # ------------------------------------------
    # DO NOT CHANGE BELOW HERE
    # ------------------------------------------

    # Generate parfiles ----------------------------
    
    # Read the base parfile
    d = read_parfile_dict(basep)

    # Generate combinations
    x, keys = combine_parameters(n)
    
    # Write parfiles
    parfile = []
    basen, ext = os.path.splitext(basep)
    for s in range(len(x)):
        for i in range(len(keys)):
            d[keys[i]] = str(x[s][i])
        # Output to file
        print(d)
        parfile.append(based+"/"+basen+"_"+str(s)+ext)
        write_parfile_dict(parfile[-1], d)
        print("Written {}".format(s))

    # Run  ----------------------------

    # Launch tasks
    pool = mp.Pool(processes=nproc)
    result_list = pool.map(run, parfile)

    for m in result_list:
        print(m)
    print("done")

    #TODO:
    # - check errors
    # - remove all parfiles that did not fail
    
    #for p in parfile:
    #    run(p)
    #os.remove(p)


