#!/usr/bin/python3

"""
Script to auto generate parfiles 

Given a template parfile with basic setup, 
generates parfiles for all combination of a given subset of parameters
(e.g. to vary binary masses and spins)

SB 10/2018
"""

import os
import itertools

from EOBUtils import *

if __name__ == "__main__": 

    # Setup -----------------------------------------
    
    # Base dir & parfile
    based = "./"
    basep = "test_NQCIteration.par"

    # Set new values/ranges for parameters (Use lists)
    q = [1., 1.2, 1.4, 2., 10.]
    chi1 = [0., 0.8]
    chi2 = [0., 0.8, 0.9]
    
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
    x, keys = combine_parameters(n)

    # For testing:
    # for s in range(len(x)):
    #     for i in range(len(keys)):
    #         d[keys[i]] = str(x[s][i])
    #         print(keys[i],d[keys[i]])
    #     print("----")
    
    # Write parfiles
    basen, ext =  os.path.splitext(basep)
    for s in range(len(x)):
        for i in range(len(keys)):
            d[keys[i]] = str(x[s][i])
        ##print(d)
        # Output to file
        write_parfile_dict(based+"/"+basen+"_"+str(s)+ext, d)
        print("Written {}".format(s))

    
