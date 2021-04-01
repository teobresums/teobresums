# Run a test
import EOBRun_module
import matplotlib.pyplot as plt
import time

# SXS:BBH:1355
f0  = 18.132161296614207
e0  = 0.0890000000
# SXS:BBH:1356
f0  = 12.328512556553836
e0  = 0.1503800000
# SXS:BBH:1358
f0  = 13.725058954557340
e0  = 0.1807800000
# SXS:BBH:1359

# SXS:BBH:324
f0 = 12.239976297170919
e0 = 0.2948000000

pars = {
    'M'                  : 50.,
    'q'                  : 1.,
    'chi1'               : 0.,
    'chi2'               : 0.,
    'Lambda1'            : 0.,
    'Lambda2'            : 0.,     
    'domain'             : 0,      #Set 1 for FD. Default = 0
    'arg_out'            : 1,      #Output hlm/hflm. Default = 0
    'use_mode_lm'        : [1],    #List of modes to use/output through EOBRunPy
    'output_lm'          : [1],    #List of modes to print on file
    'srate_interp'       : 4096.,  #srate at which to interpolate. Default = 4096.
    'use_geometric_units': 0,      #output quantities in geometric units. Default = 1
    'df'                 : 0.01,   #df for FD interpolation
    'initial_frequency'  : f0,     #in Hz if use_geometric_units = 0, else in geometric units
    'interp_uniform_grid': 1,      #interpolate mode by mode on a uniform grid. Default = 0 (no interpolation)
    'ecc'                : e0,     #Eccentricity. Default = 0.
    'ecc_freq'           : 2,      #Use periastron (0), average (1) or apastron (2) frequency for initial condition computation. Default = 1
}

#Run the WF generator
start = time.time()
t, hp, hcm, hlm, dyn = EOBRun_module.EOBRunPy(pars)
end = time.time()
DeltaT = end-start
print("Full time=%s"%DeltaT);

plt.plot(t, hp)
plt.show()
