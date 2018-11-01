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
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os, sys
import time, datetime

from EOBUtils import *

def split_out_time_info(s):
    """
    Split output and time info from run_exception() output
    TODO: improve, need to handle failed runs with return code...
    """
    T = {}    
    ss = s.split("\n");
    if "FAILED:" in s:
        T["run"] = ss[0] #"failed"
        T["parfile"] = ss[1].split(" ")[1] # run no and parfile
        T["error"] = ss[2]
        T["output"] = s
        return T
    T["run"] = "ok"
    T["error"] = []
    T["parfile"] = ss[0].split(" ")[1] # run no and parfile
    i = ss[1].split(" ")
    T["user"] = float(i[0].replace("user", ""))
    T["system"] = float(i[1].replace("system", ""))
    ms = i[2].replace("elapsed", "").split(":")
    T["elapsed"] = 60*float(ms[0]) + float(ms[1])
    T["CPU"] = int(i[3].replace("%CPU", ""))
    return T

def extract_basename(p):
    """
    Return basename and extension from parfilename
    """
    return os.path.splitext( os.path.basename(p) )
    

if __name__ == "__main__": 

    # Setup -----------------------------------------
    
    # Base dir & parfile
    based = "./"
    basep = "test_NQCIteration.par"

    # Set new values/ranges for parameters (Use lists)
    q = [1., 1.2, 1.4]
    chi1 = [0., 1.2] ## try: chi1 = [0., 1.2] to check failures
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

    # Compute timing info
    print("# Computing timing info ...")
    A = [split_out_time_info(r) for r in results] # All runs
    S = [t for t in A if t["run"]=="ok"]
    print("Sucessful_runs = {}/{}".format(len(S),len(x)))
    maxT = max(S, key=lambda x:x['elapsed'])
    minT = min(S, key=lambda x:x['elapsed'])
    avg_elapsed = sum(t['elapsed'] for t in S) / len(S)
    avg_CPU = sum(t['CPU'] for t in S) / len(S)
    print("Average_elapsed_time = "+str(avg_elapsed))
    print("Average_CPU_usage = "+str(avg_CPU))
    print("Max: {}".format(maxT))
    print("Min: {}".format(minT))
    
    # Plot wf
    fig, ax = plt.subplots()
    for p in S:
        b, e = extract_basename(p["parfile"])
        t, re, im = np.loadtxt(b+"/waveform.txt",unpack=True)
        ax.plot(t, re, '--', t, np.fabs(re**2+im**2), '-')## check/improve
    ax.set(xlabel='time', ylabel='wafeform', title=b)
    ax.grid()
    ##fig.savefig(b+".png")
    plt.show()

    # Failed runs
    print("# Info about failed runs ...")
    F = [t for t in A if t["run"]!="ok"] # failed runs    
    for f in F:
        print(f["parfile"]+" ===> "+f["error"])
        #print(f["output"]) ## uncomment for full error message

    #TODO ...
    # one could load the failed parfiles and search for intersection of parameter values ...

