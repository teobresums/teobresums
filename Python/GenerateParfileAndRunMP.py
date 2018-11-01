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
import time , datetime

from EOBUtils import *

def split_out_time_info(s):
    """
    Split output and time info from run_exception() output
    TODO: improve, need to handle failed runs with return code...
    """
    T = {}    
    s = s.split("\n");
    T["run"] = s[0].split(" ")[1] # run no and parfile
    i = s[1].split(" ")
    T["user"] = float(i[0].replace("user", ""))
    T["system"] = float(i[1].replace("system", ""))
    ms = i[2].replace("elapsed", "").split(":")
    T["elapsed"] = 60*float(ms[0]) + float(ms[1])
    T["CPU"] = int(i[3].replace("%CPU", ""))
    return T

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
        ##print("# Written {:04d}".format(s))
    print("# Written {:04d} parfiles".format(s))

    # Run  ----------------------------
    
    # Launch tasks
    print("# Running ...")
    pool = mp.Pool(processes=nproc)
    task = partial(run_exception, rm_file=1)
    results = pool.map(task, parfile)
    pool.close() 
    pool.join()
    print("# ... done")

    # Compute timing info
    print("# Computing timing info ...")
    ##TODO: currently no check for failed runs
    T = [split_out_time_info(r) for r in results]
    maxT = max(T, key=lambda x:x['elapsed'])
    minT = min(T, key=lambda x:x['elapsed'])
    avg_elapsed = sum(t['elapsed'] for t in T) / len(T)
    avg_CPU = sum(t['CPU'] for t in T) / len(T)
    print("Average elapsed time = "+str(avg_elapsed))
    print("Average CPU usage = "+str(avg_CPU))
    print("Max: {}".format(maxT))
    print("Min: {}".format(minT))

    # PLOTS (timing and wf)
    # TODO ... 
    
