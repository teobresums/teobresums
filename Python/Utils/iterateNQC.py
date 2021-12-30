#!/usr/bin/python3

"""
Script to run TEOBResumS C code for iterating NQC

Please, see the wiki page
```
    https://bitbucket.org/eob_ihes/teobresums/wiki/NQC.md
```
and understand how to do an iteration by yourself before using this script.

Example

$ python3 NQCIteration.py -i test_NQCIteration.par -n 2

Run the EOB binary set up in `test_NQCIteration.par` for 2 iterations.

NOTES

 * It assumes you have compiled the EOB C code and the exe is
   $TEOBRESUMS/TEOBResumS.x
 * It can be called with several different parfiles so to span the parameter space
 * It automatically changes some options for NQC 

SB 10/2018
RG 12/2021
"""

import os, shutil
import argparse 

from eobutils import *
from parfile import *

def generate_fnamei(fname, i):
    """
    Generates filename with absolute path and number "_{i}"
    Returns and base name
    """
    s = os.path.abspath(fname) # full path
    e = os.path.splitext(s)[1] # extension
    b = os.path.splitext(s)[0] # basename
    p = os.path.split(b)[0]    # path
    b = os.path.split(b)[1]    # name
    b = b + "_{0:02d}".format(i)
    return p+"/"+b+e, b

if __name__ == "__main__": 
    
    # Parse input
    parser = argparse.ArgumentParser(description="Run NQC iteration for an arbitrary number of binary configurations")
    parser.add_argument("-i", dest="filenames", 
                        nargs='+', required=True,
                        help="parfiles", metavar="FILE")
    parser.add_argument("-n", dest="iters", type=int, 
                        nargs='?', default=1,
                        help="Number of iterations")

    args = parser.parse_args()

    # Find/select only existing parfiles
    files_found = []
    for f in args.filenames:
        if os.path.isfile(f):
            files_found.append(f)
            print(f)
    if (args.iters<1):
        raise ValueError("I refuse to do less than 1 iteration.")

    for f in files_found:
        par = Parfile(par=f)
        
        # Make sure you have the options properly set in the parfile  
        d = {'nqc': "manual",
             'output_nqc_coefs': "yes"}              
        d0= {'nqc_coefs_flx' : "none",
             'nqc_coefs_hlm' :"compute"}
        dn= {'nqc_coefs_flx' :"fromfile",
             'nqc_coefs_hlm' : "compute"}

        for n in range(args.iters):
            print(n)
            # Parfile name for nth - iter
            parfile, basename = generate_fnamei(f, n)
            shutil.copy(f, parfile)

            # Set output directory
            d['output_dir']= '"'+ basename +'"'
            par.update_fromdict(d)

            # Set correct option for iteration
            if n == 0:
                # First iteration: do not apply NQC to flux, compute them from wf
                par.update_fromdict(d0)
            else:
                dn['nqc_coefs_flx_file']= basename_prev +'/nqc_coefs.txt'
                par.update_fromdict(dn)

            par.write_parfile_par(fname=parfile)
            # Run the code
            run(parfile)

            # Save previous iteration name
            basename_prev = basename
            
            # Remove parfile (they are copied in the simulation dir)
            os.remove(parfile)
                
            
