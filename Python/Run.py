#!/usr/bin/python3

"""
Script to run TEOBResumS C code with multiple parfiles

Examples

$ python3 Run.py test_NQCIteration.par
$ python3 Run.py *.par
$ python3 Run.py ../C/par/*.par

NOTE

 * It assumes you have compiled the EOB C code and the executable is somewehere accessible following 
   $TEOBRESUMS/TEOBResumS.x

SB 10/2018
"""

import sys, glob
import subprocess

def run(parfile):
    """
    Run TEOBResumS C code using subprocess call
    """
    x = "$TEOBRESUMS/TEOBResumS.x " + parfile
    return subprocess.call(x, shell=True)

if __name__ == "__main__": 
    
    print(sys.argv)
    for parfile in glob.glob(sys.arg):
        run(parfile)
        
