#!/usr/bin/python3

"""
Script to run timing tests WIP

SB 08/2019
"""

import os
import timeit
import subprocess

PARFILE="""
# template parfile 
# Post-adiabatic + ODE run

#
# use geometric units for I/O
# everything is given geometric units and binary mass-rescaled 
#

use_geometric_units = 0

#
# Binary properties
#

M = 2.7
q = 1. # mass ratio

use_spins = 0
chi1 = 0.0 # dimensionless spin 1 z-component 
chi2 = 0.0 # dimensionless spin 2 z-component 

initial_frequency = {f0}

tides = "TEOBRESUM3"
tides_gravitomagnetic = "PN"
LambdaAl2 = 400.     # SLy mA = mB = 1.35*Msun
LambdaBl2 = 400.
LambdaAl3 = 709.518721258149
LambdaBl3 = 709.518721258149
LambdaAl4 = 1245.540991372862
LambdaBl4 = 1245.540991372862
use_Yagi_fits = 1

centrifugal_radius = "NNLO" 
use_flm = "SSNLO"

#
# Output dir, output modes 
#

output_dir = "{outdir}"
output_hpc = 0 # switch this off to test performances!!!
output_multipoles = 0
output_lm = [ 0 ]
output_dynamics = 0

#
# Setup
#

size = 200 # size of the arrays (chunks, dynamically extended)
ode_timestep = "adaptive"
ode_stop_at_radius = 4.5 # this should be safe, but depends on Lambda
compute_LR = 0 # calculate LR ?
compute_LSO = 0 # calculate LSO ?

# interpolate on uniform grid?
#interp_uniform_grid = 0 # no 
interp_uniform_grid = 1 # yes, only h+,hx
#interp_uniform_grid = 2 # yes, all modes

srate_interp = 4096.

#
# Setup for PA
#

postadiabatic_dynamics = "yes"
postadiabatic_dynamics_N = 8 # post-adiabatic order
postadiabatic_dynamics_size = {gridsize} # grid size 
postadiabatic_dynamics_rmin = 14. # stop post-adiabatic dynamics at this radius
postadiabatic_dynamics_stop = "no" # stop after post-adiabatic dynamics
"""


if __name__ == "__main__": 

    # gridsize is actually useless
    runs = [
        {'gridsize': 800, 'f0': 10., 'outdir': 'ttest10'},
        {'gridsize': 800, 'f0': 15., 'outdir': 'ttest15'},
        {'gridsize': 800, 'f0': 20., 'outdir': 'ttest20'},
        {'gridsize': 800, 'f0': 25., 'outdir': 'ttest25'},
        {'gridsize': 800, 'f0': 35., 'outdir': 'ttest35'},
        {'gridsize': 800, 'f0': 50., 'outdir': 'ttest50'}
    ]
    
    for d in runs:
        open(os.path.join("./ttest.par"), "w").write(PARFILE.format(gridsize=d['gridsize'],f0=d['f0'],outdir=d['outdir']))
        t = timeit.timeit("os.system('$TEOBRESUMS/TEOBResumS.x ttest.par')", "import os", number=20)
        print(t/20.)
  
        # we should have a better measurement... WIP
        #open(os.path.join("./"+d['outdir']+".par"), "w").write(PARFILE.format(gridsize=d['gridsize'],f0=d['f0'],outdir=d['outdir']))
        #os.system('time $TEOBRESUMS/TEOBResumS.x ttest.par 2>'+d['outdir']+'.time')
        #r = subprocess.check_output("time $TEOBRESUMS/TEOBResumS.x "+d['outdir']+".par", shell=True, stderr=subprocess.STDOUT)
        #r = subprocess.check_output("time --output="+d['outdir']+".time $TEOBRESUMS/TEOBResumS.x ttest.par", shell=True, stderr=subprocess.STDOUT)
        #r = subprocess.call('time $TEOBRESUMS/TEOBResumS.x ttest.par', shell=True)
        
        # and add a fft estimate ... TODO
