#!/usr/bin/python3

"""
Various python utilities to run TEOBResumS
"""

import os, fileinput, glob, shutil
import argparse, re
import itertools
import subprocess

def run(parfile):
    """
    Run TEOBResumS C code using subprocess call

    NOTE

    * It assumes you have compiled the EOB C code and the exe is
      $TEOBRESUMS/TEOBResumS.x
    """
    x = "$TEOBRESUMS/TEOBResumS.x " + parfile
    return subprocess.call(x, shell=True)

def run_exception(parfile, rm_file=0):
    """
    Run TEOBResumS C code using subprocess call

    NOTE

    * It assumes you have compiled the EOB C code and the exe is
      $TEOBRESUMS/TEOBResumS.x
    * 
    * optionally delete the parfile of run if successful 
    """
    x = "$TEOBRESUMS/TEOBResumS.x " + parfile
    try:
        subprocess.check_output(x, shell=True)
    except subprocess.CalledProcessError as e:
        print(e.output)
        return(e.returncode)
    if (rm_file): os.remove(parfile)
    return "ok"
        
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

#
# Files, I/O, etc
#

def read_parfile_dict(fname):
    """
    Read a parfile into a py dictionary
    """
    d = {}
    with open(fname) as f:
        lines = [l for l in f.readlines() if l.strip()] # rm empty
        for line in lines:
            if line[0]=="#": continue
            line = line.rstrip("\n")
            line = line.split("#", 1)[0]
            key, val = line.split("=")
            d[key.strip()] = val.strip().replace('"','')
    return d

def write_parfile_dict(fname, d):
    """
    Write a parfile from a py dictionary
    """
    with open(fname, 'w') as f:
        for key, val in d.items():
            if re.match('^[a-z]+', val, re.IGNORECASE):
                f.write('%s = "%s"\n' % (key, val))
            else:
                f.write('%s = %s\n' % (key, val))
    return

def substitute_refline(fname, s1, s2):
    """
    Substitute file-line matching with regular expression in s1 with string in s2
    """
    r = re.compile(s1)
    for line in fileinput.input(fname, inplace=True):
        print( r.sub( s2, line.rstrip() ) )
    return 

def search_refline(fname, s):
    """
    Search regular expression in s in file
    """
    with open(fname) as f:
        m = re.search(s,f.read(), re.MULTILINE)
    return m 

def add_fline(fname, s):
    """
    Add the line in s to a file
    """
    with open(fname, "a") as f:
        f.write(s+"\n")
    return 

