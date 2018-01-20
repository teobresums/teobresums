#!/usr/bin/env python

""" Python warpper to test waveforms generation """

__author__ = "Sebastiano Bernuzzi"
__email__ = "sebastiano.bernuzzi@gmail.com"

import multiprocessing
import subprocess
import shlex
import sys
import numpy as np
import math

from multiprocessing.pool import ThreadPool

def dummy_run(cmd):
    """ Launch the cmd in a subprocess (thread). Dummy: echo the cmd. """
    cmd = "echo " + cmd
    p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return (out, err)
    
def run(cmd):
    """ Launch the cmd in a subprocess (thread). """
    try:
        p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        if err:
            out = None
            err = cmd
        else:
            out = cmd
            err = None
    except:
        dump_error(cmd)
        err = sys.exc_info()[0]
        out = None
        print("Sys Error: ",err)
    return (out, err)

if __name__ == '__main__':

    """ Build parameters space """

    nm = 20
    ns = 20
    mass_range = np.linspace(10., 100., num=nm)
    spin_range = np.linspace(-0.8, .8, num=ns)

    mass1 = []
    mass2 = []
    spin1 = []
    spin2 = []

    for m1 in mass_range:
        for s1 in spin_range:
            for m2 in mass_range[:nm//2]:
                for s2 in spin_range[:ns//2]:
                    # fix/exclude cases
                    # equal spins up/down
                    #if math.isclose(s1,-s2, rel_tol=1e-12):
                    #    s2 = s2 + 1e-8
                    mass1.append(m1)
                    mass2.append(m2)
                    spin1.append(s1)
                    spin2.append(s2)

    """ Run on available CPUs """

    n = multiprocessing.cpu_count()
    pool = ThreadPool(n)
    results = []
    for m1,m2,s1,s2 in zip(mass1,mass2,spin1,spin2):
        cmd = "./TEOBResumS.x -m1 %f -m2 %f -chi1 %f -chi2 %f"%(m1,m2,s1,s2)
        #results.append(pool.apply_async(dummy_run, (cmd,))) # test
        results.append(pool.apply_async(run, (cmd,)))
    pool.close()
    pool.join()

    failed = 0
    for result in results:
        out, err = result.get()
        #print("out: {} err: {}".format(out, err))
        if err is not None:
            failed += 1
            print("err: {}".format(err))

    print("#CPU = ",n)
    print("#wfs = ",len(mass1))
    print("#failed = ",failed)
    print("#mass_range = ",min(mass_range),max(mass_range))
    print("#spin_range = ",min(spin_range),max(spin_range))
    print("#done")
