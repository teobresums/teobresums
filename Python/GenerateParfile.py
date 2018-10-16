#!/usr/bin/python3

"""
Script to auto generate parfiles 

Given a template parfile with basic setup, 
generates parfiles for all combination of a given subset of parameters
(e.g. to vary binary masses and spins)

SB 10/2018
"""

from EOBUtils import *

def sub(d,n):
    """
    """
    for key in set(d.keys()) & set(n.keys()):
        # Loop on common keys
        for val in n[key]:
            print('%s: %s ' % (key, val))


if __name__ == "__main__": 

    # Setup -----------------------------------------
    
    # Base parfile
    basep = "test_NQCIteration.par"

    # Set new values/ranges for parameters
    q = [1., 1.2, 1.4, 2., 10.]
    chi1 = [0., 0.8]
    chi2 = [0., 0.8, 0.9]
    
    # Pack them into a dictionary
    # NOTE: keys must match those in parfile
    n = {'q': q,
         'chi1': chi1,
         'chi2': chi2}
    

    # DO NOT CHANGE BELOW HERE

    
    # Generate parfiles ----------------------------
    
    # Read the base parfile
    d = read_parfile_dict(basep)

    # Create the others
    sub(d,n)
    
    ##write_parfile_dict("tmp",d)
