#!/usr/bin/python3

"""
Various python utilities to run TEOBResumS.x 
"""

import os
import itertools
import subprocess
import parfile

def run(parfile):
    """
    Run TEOBResumS C code using subprocess call

    NOTE

    * It assumes you have compiled the EOB C code and the exe is
      $TEOBRESUMS/TEOBResumS.x
    """
    x = "$TEOBRESUMS/TEOBResumS.x -p " + parfile
    return subprocess.call(x, shell=True)

def run_exception(parfile, rm_file=0):
    """
    Run TEOBResumS C code using subprocess call

    NOTE

    * It assumes you have compiled the EOB C code and the exe is
      $TEOBRESUMS/TEOBResumS.x
    * Please compile with no debug/verbose options
    * Return timing info
    * Optionally delete the parfile of run if successful 
    """
    x = "echo $TEOBRESUMS/TEOBResumS.x -p '"+parfile+"'; time $TEOBRESUMS/TEOBResumS.x -p " + parfile
    try:
        p = subprocess.check_output(x, stderr=subprocess.STDOUT, shell=True, universal_newlines=True)
    except subprocess.CalledProcessError as e:
        ##print(e.output)
        return("FAILED:\n"+e.output)
        ##return(e.returncode)
    #print("Output: \n{}\n".format(p))
    if (rm_file): os.remove(parfile)
    ##if (rm_dir): os.remove(dirname)
    return p 

def combine_parameters(pars):
    """
    Given a dictionary of parameters and their values in lists, 
    Generate all the combinations (return the key list for the ordering)
    """
    x, k = [], []
    for key in pars.keys(): 
        x.append(pars[key])
        k.append(key)
    x = list(itertools.product(*x))
    return x, k

def generate_parfiles(n, based, basep):
    """
    Generate parfiles from a baseline with combinations of the specifiedd parameters.

    Returns a list with the generated parfile paths.
    """

    par = parfile.Parfile(path=based, par=basep)

    # Generate combinations
    x, keys = combine_parameters(n)

    # Write parfiles
    parfiles = []
    basen, ext =  os.path.splitext(basep)
    for s in range(len(x)):
        d = {}
        for i in range(len(keys)):
            d[keys[i]] = str(x[s][i])
        new_out = basen+"_"+str(s)
        d['output_dir'] = new_out
        print(d)
        par.update_fromdict(d)
        # Output to file
        parfiles.append(based+"/"+new_out+ext)
        par.write_parfile_par(based+"/"+basen+"_"+str(s)+ext)
        print("# Written {}".format(s))

    return parfiles

def run_batch_parfiles(parfiles):
    print("# Running ...")
    for p in parfiles:
        run(p)
        os.remove(p)
    print("done")
    return 0

def run_batch_parfiles_MP(parfiles, nproc):
    import multiprocessing as mp
    from functools import partial
    print("# Running with "+ str(nproc)+" processes ...")
    pool = mp.Pool(processes=nproc)
    task = partial(run_exception, rm_file=1)
    results = pool.map(task, parfiles)
    pool.close() 
    pool.join()   
    print("done")
    return results 
