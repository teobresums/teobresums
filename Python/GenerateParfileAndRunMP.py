#!/usr/bin/python3

"""
Script to auto generate parfiles and run them in parallel with multiprocess

Given a template parfile with basic setup, 
generates parfiles for all combination of a given subset of parameters
(e.g. to vary binary masses and spins). 
Then it runs the code with multiprocess.

SB 10/2018
"""

import multiprocessing as mp
from functools import partial
#import numpy as np
import matplotlib.pyplot as plt
import os
import sys

from EOBUtils import *

if __name__ == "__main__": 

    # Setup -----------------------------------------
    
    # Base dir & parfile
    based = "./"
    basep = "test_NQCIteration.par"
    ##basep = "tmp.par"

    # Set new values/ranges for parameters (Use lists)
    q = [1., 1.2, 1.4]
    chi1 = [0., 0.2]
    chi2 = [0., 0.4]
    
    # Pack them into a dictionary
    # NOTE: keys must match those in parfile otherwise ignored
    n = {'q': q,
         'chi1': chi1,
         'chi2': chi2}

    # Set the number of processes
    nproc = 1
    
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
        d['output_dir'] = basen+"_{:04d}".format(s) # outputdir
        # Output to file
        ##print(d)
        parfile.append(based+"/"+basen+"_{:04d}".format(s)+ext)
        write_parfile_dict(parfile[-1], d)
        print("Written {:04d}".format(s))

    # Run  ----------------------------
    
    # Launch tasks
    pool = mp.Pool(processes=nproc)
    task = partial(run_exception, rm_file=1)
    result_list = pool.map(task, parfile)
    pool.close() 
    pool.join()



