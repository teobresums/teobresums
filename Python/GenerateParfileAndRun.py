#!/usr/bin/python3

"""
Script to auto generate parfiles and run them

Given a template parfile with basic setup, 
generates parfiles for all combination of a given subset of parameters
(e.g. to vary binary masses and spins). Then it run the code.

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

    # Create a dict of cycling iterators for each key
    xv = {}
    for key in n.keys():
        xv[key] = itertools.cycle(n[key])

    # Count how many parfile
    N = 1
    for key, val in n.items():
        N *= len(val)
    print("Creating {} parfiles".format(N))
        
    # Write parfiles
    parfile = []
    basen, ext =  os.path.splitext(basep)
    for s in range(N):
        for key in set(d.keys()) & set(n.keys()):
            print('%s: %s ' % (key, xv[key].__next__() ))
            d[key] = str(xv[key].__next__())
            print(d)
        # Output to file
        parfile.append(based+"/"+basen+"_"+str(s)+ext)
        write_parfile_dict(parfile[::-1], d)
        print("Written {}".format(s))

    # Run  ----------------------------

    for p in parfile:
        run(p)
        os.remove(p)


