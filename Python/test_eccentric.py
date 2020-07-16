# Run a test
import EOBRun_module
import matplotlib.pyplot as plt

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
    'initial_frequency'  : 20.,    #in Hz if use_geometric_units = 0, else in geometric units
    'interp_uniform_grid': 2,      #interpolate mode by mode on a uniform grid. Default = 0 (no interpolation)
    'ecc'                : 0.1,    #Eccentricity. Default = 0.
}

#Run the WF generator
t, hp, hcm, hlm = EOBRun_module.EOBRunPy(pars)

plt.plot(t, hp)
plt.show()
