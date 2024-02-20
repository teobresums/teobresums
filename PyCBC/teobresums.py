""" TEOBResumS pycbc waveform plugin 
"""

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../', 'Python'))
import EOBRun_module

import numpy as np

from pycbc.types import TimeSeries
from pycbc.types import FrequencySeries

def modes_to_k(modes):
    """ Map (l, m) to linear index """
    return sorted([int(x[0]*(x[0]-1)/2 + x[1]-2) for x in modes])

DOMAIN   = {'TD':0,'FD':1}
BIT      = {'no':0, 'yes':1}
USETIDAL = {'none':0, 'TIDES_NNLO':1, 'TIDES_TEOBRESUM':2, 'TIDES_TEOBRESUM3':3}

def get_par(k,d):
    if k not in d.keys():
        raise ValueError("Need parameter {}.".format(k))
    if d[k] is None:
        raise ValueError("Need value for parameter {}.".format(k))
    return d[k] 

def teobresums_pars_update(par, domain):
    """ Update TEOBResumS default parameters with pars """

    # Initial frequency and sampling
    flow = par['f_lower'] 
    if DOMAIN['TD'] == domain:
        srate     = 1./par['delta_t']
        interp    = "yes"   # interpolate on grid with given dt
    elif DOMAIN['FD'] == domain:
        df        = par['delta_f']
        par['df'] = df
        interp    = "no"    # no interpolation needed
        srate     = 4096. # default
        if par['delta_t'] is not None:
            srate = 1./par['delta_t']
    else:
        raise ValueError("Choose TD or FD waveform.")
    
    # TEOBResumS requires some parameters 
    m1          = get_par('mass1',par)
    m2          = get_par('mass2', par)
    spin1z      = get_par('spin1z', par)
    spin2z      = get_par('spin2z', par)
    lambda1     = get_par('lambda1', par)
    lambda2     = get_par('lambda2', par)
    distance    = get_par('distance', par)
    inclination = get_par('inclination', par)
    coa_phase   = get_par('coa_phase', par)
    ecc         = get_par('eccentricity', par)
    # NOTE: In teob `anomaly` is actually not the mean anomaly, but the true one
    #       since the parameters are linked with one another by Kepler Equation, 
    #       and have similar physical meanings we use one in place of the other 
    anomaly     = get_par('mean_per_ano', par)
    
    spin1x = spin1y = spin2x = spin2y = 0.
    if par['spin1x'] is not None:
        spin1x = par['spin1x']
    if par['spin1y'] is not None:
        spin1y = par['spin1y']
    if par['spin2x'] is not None:
        spin2x = par['spin2x']
    if par['spin2y'] is not None:
        spin2y = par['spin2y']

    if par['mode_array'] is None:
        k = modes_to_k([(2,2)])
    else:
        k = modes_to_k(par['mode_array'])
        if k[0] < 0 or k[-1] > 34:
            raise ValueError("Invalid mode list.")
    
    # Check masses convention
    q = m1/m2
    if q < 1.0 :
        m1,m2 = m2,m1
        spin1z,spin2z = spin2z,spin1z
        lambda1,lambda2 = lambda2,lambda1
        q = 1./q

    # Always use physical units
    par['use_geometric_units'] = "no"
    
    # Set TEOBResumS parameters
    # Below we list all possible parameters as a reference for an 
    # advanced used, though they are set to None and removed before
    # updating the dictionary with the use input 'pars'.
    # If those values are not set by user, then internal default are
    # used. 
    default = {
        'domain' : domain, # 0 = TD, 1 = FD 
        'M' : m1+m2,
        'q' : q,
        'LambdaAl2' : lambda1,
        'LambdaBl2' : lambda2,     
        'chi1' : spin1z,
        'chi2' : spin2z,
        #
        'chi1x': spin1x,
        'chi1y': spin1y,
        'chi1z': spin1z,
        'chi1x': spin2x,
        'chi1y': spin2y,
        'chi1z': spin2z,
        #
        'distance'           : distance, 
        'inclination'        : inclination,
        'coalescence_angle'  : coa_phase, # reference angle/phase at coalescence
        'srate_interp'       : srate, # in Hz, srate at which to interpolate output
        'initial_frequency'  : flow,  # in Hz if use_geometric_units = 0, else in geometric units
        'use_mode_lm'        : k, # List of wvf modes to use
        #
        'use_geometric_units': 1,  # I/O units output: 1 = geometric, 0 = physical
        'interp_uniform_grid': interp, # Interpolate mode by mode on a uniform grid. Default = 0 (no interpolation)
        'arg_out'            : "no", # return modes hlm/hflm. Default = 0 (no)
        #
        'ringdown_extend_array' : None,
        'centrifugal_radius'    : None,
        'use_flm'               : None,
        'nqc'                   : None,
        'nqc_coefs_flx'         : None,
        'nqc_coefs_hlm'         : None,
        'use_speedytail'        : None,
        #
        'use_tidal'             : None, # Tidal model, BBH default = TIDES_OFF, BNS default = TIDES_TEOBRESUM3 (see USETIDAL)
        'use_tidal_gravitomagnetic' : None, # Gravitomagnetic tides model, BBH default = TIDES_GM_OFF, BNS default = TIDES_GM_PN
        'pGSF_tidal'            : None, # default = 4.0
        #
        'output_hpc'        : "no",
        'output_lm'         : [-1], 
        'output_multipoles' : "no",
        'output_dynamics'   : "no",
        #
        'compute_LR'        : None,
        'compute_LR_guess'  : None,
        'compute_LSO'       : None,
        'compute_LSO_guess' : None, 
        #
        'postadiabatic_dynamics'      : None,
        'postadiabatic_dynamics_N'    : None,
        'postadiabatic_dynamics_size' : None,
        'postadiabatic_dynamics_stop' : None,
        'postadiabatic_dynamics_rmin' : None,
        #
        'r0' : None,
        'dt' : None,
        'dt_merger_interp' : None,
        'dt_interp' : None,
        'srate' : None,
        'size' : None,
        #
        'ode_timestep' : None,
        'ode_abstol' : None,
        'ode_reltol' : None,
        'ode_tmax' : None,
        'ode_stop_radius' : None,
        'ode_stop_afterNdt' : None,
        #
        'ecc'     : ecc, #this actually affects the evolution *only* if on teobresums-eccentric branch
        'anomaly' : anomaly
    }
    
    # Remove unset parameters (use internal defaults)
    # and return updated parameter list
    rm_none = {k: v for k, v in default.items() if v is not None}
    default.clear()
    default.update(rm_none)
    default.update(par)
    return default

    
def teobresums_td(**par):
    """ Time-domain TEOBResumS """

    # Required parameters pycbc
    flow = par['f_lower'] 
    dt   = par['delta_t']

    # TEOBResumS parameters
    par = teobresums_pars_update(par, DOMAIN['TD'])

    # Run the WF generator
    t, hp, hc  = EOBRun_module.EOBRunPy(par)
    wf = hp - 1j* hc
    
    # Return product should be a pycbc time series in this case for
    # each GW polarization
    #
    # Note that by convention, the time at 0 is a fiducial reference.
    # For CBC waveforms, this would be set to where the merger occurs

    offset = t[0]
    wf     = TimeSeries(wf, delta_t=dt, epoch=offset)
    return wf.real(), wf.imag()


def teobresums_fd(**par):
    """ Frequency-domain TEOBResumS """

    # Required parameters pycbc
    flow = par['f_lower'] 
    df   = par['delta_f']

    # TEOBResumS parameters
    par = teobresums_pars_update(par, DOMAIN['FD'])

    # Run the WF generator
    f, rhp, ihp, rhc, ihc  = EOBRun_module.EOBRunPy(par)
    hp = rhp -1j*ihp
    hc = rhc -1j*ihc

    offset = 0 # GPS time, fixed at 0

    # we need to extend the TEOBResumS output, which automatically 
    # starts at f_low instead of 0    
    N = int(flow/df)
    hp = np.pad(hp, (N,0),  mode='constant')
    hc = np.pad(hc, (N,0),  mode='constant')

    hpf = FrequencySeries(hp, delta_f=df, epoch=offset)
    hcf = FrequencySeries(hc, delta_f=df, epoch=offset)
    return hpf, hcf
