# Run a test

import os
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(1,'../')
import EOBRun_module

def modes_to_k(modes):
    """
    Map multipolar (l,m) -> linear index k
    """
    return [int(x[0]*(x[0]-1)/2 + x[1]-2) for x in modes]

# Routines to map mem usage

import psutil
def memory_usage_psutil():
    process = psutil.Process(os.getpid())
    mem = process.memory_info()[0] / float(2 ** 20)
    return mem

import resource
def memory_usage_resource():
    rusage_denom = 1024.
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem


if __name__ == "__main__":
    
    # Example run

    # Use 22 and 31
    modes = [[2,2], [3,1]]
    k = modes_to_k(modes)
    print(k)

    # Define a dictionary with input parameters
    pars = {
        'M'                  : 50.,
        'q'                  : 1.2,
        'chi1'               : 0.,
        'chi2'               : 0.,
        'LambdaAl2'          : 0.,
        'LambdaBl2'          : 0.,     
        'domain'             : 0,      #Set 1 for FD. Default = 0
        'arg_out'            : "yes",      #Output hlm/hflm. Default = 0
        'use_mode_lm'        : k,      #List of modes to use/output through EOBRunPy
        'output_lm'          : k,      #List of modes to print on file
        'srate_interp'       : 4096.,  #srate at which to interpolate. Default = 4096.
        'use_geometric_units': "no",   #output quantities in geometric units. Default = "yes"
        'initial_frequency'  : 20.,    #in Hz if use_geometric_units = 0, else in geometric units
        'interp_uniform_grid': "yes"   #interpolate mode by mode on a uniform grid. Default = "no" (no interpolation)
    }

    # Run the WF generator (TD)
    t, hp, hcm, hlm, dyn = EOBRun_module.EOBRunPy(pars)

    # Plot h+
    plt.plot(t, hp)
    plt.show()

    # Plot Re[h_22]
    Ah22   = hlm['1'][0]
    Phih22 = hlm['1'][1]
    plt.plot(t, Ah22*np.cos(Phih22)) 
    plt.show()

    # To test for leaks, uncomment below
    # print('Memory test')
    # for i in range(100):
    #     t, hp, hcm, hlm, dyn = EOBRun_module.EOBRunPy(pars)
    #     print(memory_usage_psutil(),memory_usage_resource())

