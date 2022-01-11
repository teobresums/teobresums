#!/usr/bin/python3

"""
Script to run the benchmarks. Uses the parfiles.json under teobresums/Python/Utils

RG 01/2022
"""

import eobutils as utils; import json; import parfile as p

# read the parfiles
with open("./parfiles.json", "r") as f:
            data = json.load(f)

# write, run and delete
parlist = []
for s in data['examples']:
    print("Writing {}.par ".format(s['name']))
    par = p.Parfile(par=s)
    parname = fname=s['name']+".par"
    par.write_parfile_par(parname)
    parlist.append(parname)

utils.run_batch_parfiles(parlist)


