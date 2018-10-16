#!/usr/bin/python3

"""
Script to auto generate parfiles 

Given a template parfile with basic setup, 
generates parfiles for all combination of a given subset of parameters
(e.g. to vary binary masses and spins)

SB 10/2018
"""

import os
from EOBUtils import *

class arrayiter:
    """
    Iterator for array 
    Keep on returning last element
    """
    def __init__(self, a):
        self.a = a
        self.nmax = len(a)
        self.n = 0
    def __next__(self):
        if self.n < self.nmax-1:
            self.n += 1
            return self.a[self.n]
        else:
            return self.a[self.nmax-1]
            #raise StopIteration

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
    # NOTE: keys must match those in parfile
    n = {'q': q,
         'chi1': chi1,
         'chi2': chi2}
    
    # ------------------------------------------
    # DO NOT CHANGE BELOW HERE
    # ------------------------------------------
    
    # Generate parfiles ----------------------------
    
    # Read the base parfile
    d = read_parfile_dict(basep)

    # Create a dict of iterators for each key

    xv = {}
    for key in n.keys():
        xv[key] = arrayiter(n[key])

    # Count how many parfile
    N = 1
    for key, val in n.items():
        N *= len(val)
    print("Creating {} parfiles".format(N))
        
    # Write parfiles
    basen, ext =  os.path.splitext(basep)
    for s in range(N):
        for key in set(d.keys()) & set(n.keys()):
            #print('%s: %s ' % (key, xv[key].__next__() ))
            d[key] = xv[key].__next__()    
        #print(d)
        #write_parfile_dict(based+"/"+"basen_"+str(s)+ext, d)
        #print("Written {}".format(s))

    
    
