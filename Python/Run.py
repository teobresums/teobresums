#!/usr/bin/python3

"""
Script to run TEOBResumS C code with multiple parfiles

Examples

$ python3 Run.py test_NQCIteration.par
$ python3 Run.py *.par
$ python3 Run.py ../C/par/*.par

SB 10/2018
"""

import sys, glob

from EOBUtils import *

if __name__ == "__main__": 
    
    print(sys.argv)
    for parfile in glob.glob(sys.arg):
        run(parfile)
        
