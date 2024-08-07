""" 
Example running for and eccentric, precessing BBH system.
This plot generates the EOB waveform shown in Fig. 12 of Gamba:2024cvy.

If you use this model, please cite:
@article{Gamba:2024cvy,
    author = {Gamba, Rossella and Chiaramello, Danilo and Neogi, Sayan},
    title =  {Towards efficient Effective One Body models for generic, non-planar orbits},
    eprint = {2404.15408},
    archivePrefix = {arXiv},
    primaryClass = {gr-qc},
    year = {2024},
}

The initial frequency is twice the orbit-averaged frequecy expressed in geometric units.
Eccentricity and true anomaly are specify uniquely the evolution of the system.
"""

import EOBRun_module
import matplotlib.pyplot as plt; import numpy as np; import time

if __name__ == '__main__':

    q      = 1.
    chi1   = [0.31708349, 0.62449976, 0.] # chi1, chi2 rotated around L; |chi1| = |chi2| = 07
    chi2   = [0.31708347, 0.62449972, 0.]
    f0     = 0.00371968856863
    e0     = 0.3168941
    z0     = 3.51120126

    pars = {
        # System parametes    
        'M'                  : 1.,    
        'q'                  : q,       # Mass ratio m1/m2 > 1
        'chi1x'              : chi1[0], # x component of chi1
        'chi1y'              : chi1[1], # y component of chi1   
        'chi1z'              : chi1[2], # z component of chi1
        'chi2x'              : chi2[0], # x component of chi2
        'chi2y'              : chi2[1], # y component of chi2
        'chi2z'              : chi2[2], # z component of chi2
        'chi1'               : chi1[2],
        'chi2'               : chi2[2],
        'LambdaAl2'          : 0.,
        'LambdaBl2'          : 0.,
        'ecc'                : e0,      # Eccentricity. Default = 0.
        'anomaly'            : z0,      # True anomaly. Default = Pi
        'distance'           : 1.,

        # Initial conditions and output time grid
        'domain'             : 0,       # Time domain. EOBSPA is not available for eccentric waveforms!
        'use_geometric_units': "yes",   # output quantities in geometric units. Default = 1
        'initial_frequency'  : f0,      # in Hz if use_geometric_units = 0, else in geometric units
        'interp_uniform_grid': "yes",   # interpolate mode by mode on a uniform grid. Default = 0 (no interpolation)
        'srate_interp'       : 4096.,

        # Modes
        'use_mode_lm'        : [0,1],  # List of modes to use/output through EOBRunPy [(2,|2|), (2,|1|)]

        # Output parameters (Python)
        'arg_out'            : "yes",   # Output hlm/hflm. Default = 0
        'output_hpc'         : "no",
        
        # Precession settings (do not touch!)
        'spin_flx'           : 'EOB',   # Use the EOB frequency evolution to drive the PN spins dynamics
        'spin_interp_domain' : 0        # interpolate in time rather than frequency
    }

    print("Run the WF generator...")
    start = time.time()
    t, _, _, hlm, dyn = EOBRun_module.EOBRunPy(pars)
    end = time.time()
    DeltaT = end-start
    print("Full time=%s"%DeltaT)

    print("Plot...")
    _, (ax1, ax2)  = plt.subplots(2, 1, sharex=True)
    for k, ax in zip(['1', '0'], [ax1, ax2]):
        ax.plot(t, hlm[k][0]*np.cos(hlm[k][1]), linestyle='-', color='k')
        ax.plot(t, hlm[k][0], color='gray', linestyle='--')
        ax.grid()
    ax2.set_xlabel("t/M")
    ax1.set_ylabel(r'$h_{22}/\nu$')
    ax2.set_ylabel(r'$h_{21}/\nu$')
    plt.show()
