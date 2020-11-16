# Produce a waveform bhns test

import sys
sys.path.insert(1,'../')
import EOBRun_module

import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import minimize
from scipy import optimize
import cmath
import math
import array

# Source parameters
M1 =  8.4 # Mo
M2 = 1.4 # Mo
Deff = 56 # Mpc
iota = 0.
#iota = 150./360.*np.pi

# Use 22 mode only
def modes_to_k(modes):
    return [int(x[0]*(x[0]-1)/2 + x[1]-2) for x in modes]


k = modes_to_k([[2,2]])

# Define a dictionary with input parameters
pars = {
    'M'                  : M1+M2,
    'q'                  : M1/M2,
    'Lambda1'            : 0.,
    'Lambda2'            : 0.,  
    'chi1'               : 0.,
    'chi2'               : 0,
    'domain'             : 0,      # TD
    'arg_out'            : 1,      # Output hlm/hflm. Default = 0
    'use_mode_lm'        : k,      # List of modes to use/output through EOBRunPy
    'output_lm'          : k,      # List of modes to print on file
    'srate_interp'       : 4096.,  # srate at which to interpolate. Default = 4096.
    'use_geometric_units': 1,      # output quantities in geometric units. Default = 1
    'initial_frequency'  : 0.005,    # in Hz if use_geometric_units = 0, else in geometric units 
    'interp_uniform_grid': 2,      # Interpolate mode by mode on a uniform grid. Default = 0 (no interpolation)
    'distance'           : Deff,
    'inclination'        : iota,
    'output_dynamics'    : 1,
    'output_hpc'         : 1
}

print("Read pars..")

# Run the WF generator
t, hp, hc, hlm = EOBRun_module.EOBRunPy(pars)

print("EOBRunPy done..")

# Approximate merger time
h = hp - 1j*hp
tmrg = t[np.argmax(np.abs(h))]
print('tmrg = {}'.format(tmrg))
#t = t - tmrg # shift to mrg
q = M1/M2
nu = q/((1.+q)*(1.+q))

# Plot h
plt.plot(t, hp, label=r'$h_{+}$')
plt.plot(t, hc, label=r'$h_{\times}$')
#plt.xlim([-0.4,0.05])
#plt.xlim([0.,5500])
plt.xlabel(r'$t-t_{\rm mrg}$ $[s]$')
plt.legend()
#plt.savefig('plots/h.png')
plt.show()

# 22 mode
A22   = hlm['1'][0]
Phi22 = hlm['1'][1]

# t/M, geometric units
GMo = 4.925490947*1e-6 # s
#t *= 1./((M1+M2)*GMo)

# GW frequency
omg22 = np.zeros_like(Phi22)
omg22[1:] = np.diff(Phi22)/np.diff(t)

# Plot 22 mode in geometric units
plt.plot(t, A22*np.cos(Phi22)*nu, '--', label=r'$\Re{(Rh_{22}/(M))}$') 
plt.plot(t, A22*nu, label=r'$|Rh_{22}/(\nu M)|$')
plt.plot(t, omg22, label=r'$M\omega_{22}$')
#plt.xlim([-500.,100])
#plt.xlim([0.,3000])
plt.xlabel(r'$t/M$')
plt.legend()
#plt.savefig('plots/amp_phs_frq.png')
plt.show()

