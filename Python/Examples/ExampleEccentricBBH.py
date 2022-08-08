""" 
Example running for some eccentric inspirals. The script is set so to reproduce EOB runs in Table III of arXiv:2101.08624

@article{Nagar:2021gss,
author = "Nagar, Alessandro and Bonino, Alice and Rettegno, Piero",
title = "{All in one: effective one body multipolar waveform model for spin-aligned, quasi-circular, eccentric, hyperbolic black hole binaries}",
eprint = "2101.08624",
archivePrefix = "arXiv",
primaryClass = "gr-qc",
month = "1",
year = "2021"}

The initial apastron frequency is expressed in Hz and the code runs by default in physical units.
The total reference mass is fixed to M=50*Msun. The dimensionless apastron frequency omg_a is transformed into
the physical apastron frequency f0 according to the formula

f = omg_a/(2*pi*M)

where Msun = 4.925490947e-6 sec.

Running times for these initial data are of the order of 0.1 sec on  a 3.5GHz IntelCore i7, 16G RAM,
except for 1149 and 1169 when one gets arund 0.2sec. Initial physical frequency is always below 20Hz.
"""

import EOBRun_module
import matplotlib.pyplot as plt; import numpy as np; import time

# SXS waveforms and corresponding initial frequency and eccentricity
# ID : [f0 at 50 Msun (Hz), e0, q, s1z, s2z]
SXS_params = {
    'SXS:BBH:1355' : [18.132161296614207, 0.0890000000,     1., 0., 0.],      # q=1, omg_a = 0.02805750
    'SXS:BBH:1356' : [12.328512556553836, 0.1503800000,     1., 0., 0.],      # q=1, omg_a = 0.019077
    'SXS:BBH:1357' : [12.666501342373286, 0.19201,          1., 0., 0.],      # q=1, omg_a = 0.01960
    'SXS:BBH:1358' : [13.725058954557340, 0.1807800000,     1., 0., 0.],      # q=1, omg_a = 0.021238
    'SXS:BBH:1359' : [13.823288964967581, 0.18240,          1., 0., 0.],      # q=1, omg_a = 0.02139
    'SXS:BBH:1360' : [12.6600388, 0.2429,                   1., 0., 0.],      # q=1, omg_a = 0.01959
    'SXS:BBH:1361' : [13.571251438257093, 0.23557,          1., 0., 0.],      # q=1, omg_a = 0.0210
    'SXS:BBH:1362' : [12.369226310868608, 0.3019,           1., 0., 0.],      # q=1, omg_a = 0.01914
    'SXS:BBH:1363' : [12.330451306759301, 0.30479,          1., 0., 0.],      # q=1, omg_a = 0.01908
    'SXS:BBH:1364' : [16.305535478031654, 0.08464,          2., 0., 0.],      # q=2, omg_a = 0.025231
    'SXS:BBH:1365' : [15.501600392832042, 0.11015,          2., 0., 0.],      # q=2, omg_a = 0.023987
    'SXS:BBH:1366' : [16.673251767001570, 0.1496,           2., 0., 0.],      # q=2, omg_a = 0.02580
    'SXS:BBH:1367' : [16.673251767001570, 0.15065,          2., 0., 0.],      # q=2, omg_a = 0.02580
    'SXS:BBH:1368' : [16.330739230702701, 0.14951,          2., 0., 0.],      # q=2, omg_a = 0.02527
    'SXS:BBH:1369' : [11.212438688274313, 0.3134,           2., 0., 0.],      # q=2, omg_a = 0.01735
    'SXS:BBH:1370' : [10.931319908481843, 0.31445,          2., 0., 0.],      # q=2, omg_a = 0.016915
    'SXS:BBH:1371' : [18.778734490136888, 0.0912,           3., 0., 0.],      # q=3, omg_a = 0.029058
    'SXS:BBH:1372' : [16.847739285493446, 0.14915,          3., 0., 0.],      # q=3, omg_a = 0.026070
    'SXS:BBH:1373' : [16.350126732757353, 0.15035,          3., 0., 0.],      # q=3, omg_a = 0.0253
    'SXS:BBH:1374' : [10.951353660604985, 0.31388,          3., 0., 0.],      # q=3, omg_a = 0.016946
    'SXS:BBH:89'   : [11.496788718409222, 0.07201,          1.,-0.5, 0.],     # (q=1,-0.50,0), omg_a = 0.01779
    'SXS:BBH:324'  : [12.239976297170919, 0.2948000000,     1.22, 0.33,-0.44],# (q=1.22,+0.33,-0.44), omg_a = 0.01894
    'SXS:BBH:1149' : [17.216101824531854, 0.0623,           3.,0.70,0.60],    # (q=3,+0.70,+0.60), omg_a = 0.02664
    'SXS:BBH:1169' : [15.694182913241596, 0.04895,          3,-0.70,-0.60]    # (q=3,-0.70,-0.60), omg_a = 0.024285
}

if __name__ == '__main__':

    ID     = 'SXS:BBH:1149'
    M      = 50.
    f0_ref, e0, q, s1, s2 = SXS_params[ID]

    # rescale the initial frequency to the desired mass
    f0     = f0_ref*(50/M) 

    print("Compute the EOB waveform for ", ID, '...')

    pars = {
        # System parametes, assuming aligned spins        
        'M'                  : M,      # Total mass
        'q'                  : q,      # Mass ratio m1/m2 > 1
        'chi1'               : s1,     # Z component of chi_1
        'chi2'               : s2,     # Z component of chi_2
        'LambdaAl2'          : 0.,     # Quadrupolar tidal parameter of body 1 (A)
        'LambdaBl2'          : 0.,     # Quadrupolar tidal parameter of body 2 (B)
        'ecc'                : e0,     # Eccentricity. Default = 0.
        'ecc_freq'           : 2,      # Use periastron (0), average (1) or apastron (2) frequency for initial condition computation. Default = 1

        # Initial conditions and output time grid
        'domain'             : 0,      # Time domain. EOBSPA is not available for eccentric waveforms!
        'srate_interp'       : 4096.,  # srate at which to interpolate. Default = 4096.
        'use_geometric_units': "no",   # output quantities in geometric units. Default = 1
        'initial_frequency'  : f0,     # in Hz if use_geometric_units = 0, else in geometric units
        'interp_uniform_grid': "yes",  # interpolate mode by mode on a uniform grid. Default = 0 (no interpolation)

        # Modes
        'use_mode_lm'        : [1],    # List of modes to use/output through EOBRunPy

        # Output parameters (Python)
        'arg_out'            : "yes",  # Output hlm/hflm. Default = 0
    }

    print("Run the WF generator...")
    start = time.time()
    t, hp, hc, hlm, dyn = EOBRun_module.EOBRunPy(pars)
    end = time.time()
    DeltaT = end-start
    print("Full time=%s"%DeltaT);

    A       = np.abs(hp -1j*hc)
    p       = -np.unwrap(np.angle(hp -1j*hc))
    omg     = np.zeros_like(p)
    omg[1:] = np.diff(p)/np.diff(t)

    print("Plot...")
    _, ax1  = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(t, hp, color='b', linestyle='-')
    ax1.plot(t, A, color='b')
    ax2.plot(t, omg, color='r')
    ax1.set_xlabel("t [s]")
    ax1.set_ylabel(r'$h_+$', color='b')
    ax2.set_ylabel('$\omega$', color='r')
    plt.title(ID)
    plt.show()
