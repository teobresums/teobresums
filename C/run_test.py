# Run a test
import EOBRun_module
import matplotlib.pyplot as plt
import numpy as np
import os
import psutil
import resource

def memory_usage_psutil():
    process = psutil.Process(os.getpid())
    mem = process.memory_info()[0] / float(2 ** 20)
    return mem

def memory_usage_resource():
    rusage_denom = 1024.
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem

def modes_to_k(modes):
    return [int(x[0]*(x[0]-1)/2 + x[1]-2) for x in modes]
# Example run

# use 22 and 31
modes = [[2,2], [3,1]]
k = modes_to_k(modes)
print(k)

# Define a dictionary with input parameters
pars = {
    'M'                  : 50.,
    'q'                  : 1.,
    'chi1'               : 0.,
    'chi2'               : 0.,
    'Lambda1'            : 0.,
    'Lambda2'            : 0.,     
    'domain'             : 0,      #Set 1 for FD. Default = 0
    'arg_out'            : 1,      #Output hlm/hflm. Default = 0
    'use_mode_lm'        : k,      #List of modes to use/output through EOBRunPy
    'output_lm'          : k,      #List of modes to print on file
    'srate_interp'       : 4096.,  #srate at which to interpolate. Default = 4096.
    'use_geometric_units': 0,      #output quantities in geometric units. Default = 1
    'df'                 : 0.01,   #df for FD interpolation
    'initial_frequency'  : 20.,    #in Hz if use_geometric_units = 0, else in geometric units
    'interp_uniform_grid': 2,      #interpolate mode by mode on a uniform grid. Default = 0 (no interpolation)
    'ecc'                : 0.1,    #Eccentricity. Default = 0.
    'j_hyp'              : 0.,     #J_hyp. Default = 0.
    'r_hyp'              : 0.,     #r_hyp. Default = 0.
    'H_hyp'              : 0.,     #H_hyp. Default = 0.
}

#Run the WF generator
t, hp, hcm, hlm = EOBRun_module.EOBRunPy(pars)

# Plot h+
plt.plot(t, hp)
plt.show()

# Plot Re[h_22]
Ah22   = hlm['1'][0]
Phih22 = hlm['1'][1]
plt.plot(t, Ah22*np.cos(Phih22)) 
plt.show()

# If arg_out = 0
#t, hp, hcm = EOBRun_module.EOBRunPy(pars)          

# If FD
#f, hpr, hpi, hcr, hci = EOBRun_module.EOBRunPy(pars) 

#to test for leaks, uncomment below
#print(memory_usage_psutil(),memory_usage_resource())

