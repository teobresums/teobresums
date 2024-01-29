# Produce a waveform compatible with GW150914

import sys
sys.path.insert(1,'../')
import EOBRun_module

import numpy as np
import matplotlib.pyplot as plt

GMo = 4.925490947*1e-6 # s
# Map HMs (l,m) indices to internal TEOBResumS indices convention
def modes_to_k(modes):
    return [int(x[0]*(x[0]-1)/2 + x[1]-2) for x in modes]

"""

This file contains a minimal aligned-spins BBH example to generate a GW150914-like waveform using TEOBResumS,
both in the time and frequency domain.

"""

# Parameters from: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.061102
M1   = 36 # Solar masses
M2   = 29 # Solar masses
Deff = 410 # Mpc
iota = 150./360.*np.pi # radiants

# Use 22 mode only
k = modes_to_k([[2,2]])

TEOBResumS_domain = {'TD':0,  'FD':1}
YESNO2INT         = {'yes':1, 'no':0}
domain = 'TD'
seglen = 8.0 #s

# Define a dictionary with input parameters
pars = {
    # System parametes, assuming aligned spins
    'M'                  : M1+M2,
    'q'                  : M1/M2,  
    'chi1'               : 0.,
    'chi2'               : 0.,
    'distance'           : Deff,
    'inclination'        : iota,
    'coalescence_angle'  : 0.0,

    # Initial conditions and output time grid
    'use_geometric_units': "no",   # Output quantities in geometric units
    'initial_frequency'  : 35.,    # in Hz if use_geometric_units = "no", else in geometric units
    'domain'             : TEOBResumS_domain[domain],
    'interp_uniform_grid': "yes",  # Interpolate the output waveform a uniform grid, mode by mode
    'srate_interp'       : 4096.,  # srate at which to interpolate, fixes f_max in 'FD' too
    # Modes
    'use_mode_lm'        : k,      # List of modes to use/output through EOBRunPy
    'model'              : "Giotto",
    # Output parameters
    'arg_out'            : "yes",      # Request multipoles and dynamics as output of the function call. Default is "no". Allowed values: ["no","yes"].
}

# Run the WF generator
if YESNO2INT[pars['arg_out']]:
    t, hp, hc, hlm, dyn = EOBRun_module.EOBRunPy(pars)
else:
    t, hp, hc           = EOBRun_module.EOBRunPy(pars)

# Approximate merger time
h    = hp - 1j*hc
tmrg = t[np.argmax(np.abs(h))]
t    = t - tmrg # shift to mrg

# Plot the result
plt.figure()
plt.plot(t, hp, label=r'$h_{+}$')
plt.plot(t, hc, label=r'$h_{\times}$')
plt.xlim([-0.4,0.05])
plt.xlabel(r'$t-t_{\rm mrg}$ $[s]$')
plt.legend()
plt.savefig('h_TD.png')

if YESNO2INT[pars['arg_out']]:

    # t/M, geometric units
    t *= 1./((M1+M2)*GMo) 

    # hlm are a dict containing amplitude and phase, Amp, as unpacked below.
    k     = modes_to_k([[2,2]])[0]
    A22   = hlm['{}'.format(k)][0]
    Phi22 = hlm['{}'.format(k)][1]
    # GW frequency
    omg22     = np.zeros_like(Phi22)
    omg22[1:] = np.diff(Phi22)/np.diff(t)

    # Plot 22 mode in geometric units
    plt.figure()
    plt.plot(t, A22*np.cos(Phi22), '--', label=r'$\Re{(Rh_{22}/(M\nu))}$')
    plt.plot(t, A22, label=r'$|Rh_{22}/(\nu M)|$')
    plt.plot(t, omg22, label=r'$M\omega_{22}$')
    plt.xlim([-500.,100])
    plt.xlabel(r'$t/M$')
    plt.legend()
    plt.savefig('h22_TD.png')

# FD
pars['domain'] = TEOBResumS_domain['FD']
pars['df']     = 1.0/seglen

if YESNO2INT[pars['arg_out']]:
    f, hp_re, hp_im, hc_re, hc_im, hflm, htlm, dyn = EOBRun_module.EOBRunPy(pars)
else:
    f, hp_re, hp_im, hc_re, hc_im                  = EOBRun_module.EOBRunPy(pars)

# Spell out waveform conventions
hp, hc = hp_re-1j*hp_im, hc_re-1j*hc_im
h      = hp - 1j*hc

# Plot the result
plt.figure()
plt.plot(f, hp_re, label=r'$Re(h_{+})$')
plt.plot(f, hc_re, label=r'$Re(h_{\times})$')
plt.xlabel(r'$f$ $[Hz]$')
plt.xlim([pars['initial_frequency'],750])
plt.legend()
plt.savefig('h_re_FD.png')

plt.figure()
plt.plot(f, hp_im, label=r'$Im(h_{+})$')
plt.plot(f, hc_im, label=r'$Im(h_{\times})$')
plt.xlabel(r'$f$ $[Hz]$')
plt.xlim([pars['initial_frequency'],750])
plt.legend()
plt.savefig('h_im_FD.png')