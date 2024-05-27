import numpy as np; from math import factorial as fact
import matplotlib.pyplot as plt

import EOBRun_module as EOB

def CreateDict( M=1.0, q=1.0,
                chi1=[0,0,0], chi2=[0,0,0],
                lambda1=0, lambda2=0,
                ecc=0., anomaly=np.pi,
                iota=0, distance=1.,
                f0=20., srate=4096*2., df=0.05,
                interp='yes', domain=0,
                modes=[1], coa=0, argout="no"
                ):
    """
    Create the dictionary of parameters for EOBRunPy
    ---
    Input:
        - M        : total mass of the system [Msun]
        - q        : mass ratio, m1/m2 > 1
        - chi1     : spin vector of the primary
        - chi2     : spin vector of the secondary
        - lambda1  : quadrupolar tidal parameter of the primary
        - lambda2  : quadrupolar tidal parameter of the secondary
        - iota     : inclination angle between the orbital angular momentum and the line of sight [rad]
        - distance : luminosity distance to the source [Mpc]
        - f0       : initial frequency of the waveform [Hz]
        - srate    : sampling rate [Hz]
        - df       : frequency step [Hz]
        - interp   : interpolate the waveform to a uniform grid [yes/no]
        - domain   : domain of the waveform [0/1] for time/frequency
        - modes    : list of co-precessing modes to employ
        - coa      : coalescence angle [rad]
        - argout   : output additional arguments (hlm, dynamics) [yes/no]
    
    Output:
        - pardic   : dictionary of EOB parameters, to be fed to the wf generator
    """

    chi1x, chi1y, chi1z = chi1
    chi2x, chi2y, chi2z = chi2

    pardic = {
        'M'                  : M,
        'q'                  : q,
        'chi1x'              : chi1x,
        'chi1y'              : chi1y,
        'chi1z'              : chi1z,
        'chi2x'              : chi2x,
        'chi2y'              : chi2y,
        'chi2z'              : chi2z,
        'chi1'               : chi1z,
        'chi2'               : chi2z,
        'LambdaAl2'          : lambda1,
        'LambdaBl2'          : lambda2,
        'distance'           : distance,
        'initial_frequency'  : f0,
        'use_geometric_units': "no",
        'coalescence_angle'  : coa,
        'interp_uniform_grid': interp,
        'df'                 : df,
        'domain'             : domain,
        'srate_interp'       : srate,
        'inclination'        : iota,
        'use_mode_lm'        : modes,      # List of modes to use/output through EOBRunPy
        'output_lm'          : modes,
        'arg_out'            : argout,
        'output_hpc'         : "no",
        'output_multipoles'  : "no", 
        'output_dynamics'    : "no",
        'time_shift_FD'      : "no",
        'time_shift_TD'      : "no",
        'ecc'                : ecc,
        'anomaly'            : anomaly
    }
    return pardic

def gen_wf(m1, m2, s1z, s2z, lam1, lam2, additional_pars={}, return_zero=True):
    """
    Wrapper to EOBRunPy
    ---
    Input:
        - m1, m2, s1z, s2z, lam1, lam2: masses, spins and tidal parameters of the binary
        - additional_pars: dictionary with additional options
    Output:
        Depending on the value of par['argout'], par['domain'] returns:
            - (no, 0)  : t, hp, hc
            - (yes, 0) : t, hp, hc, hlm, dyn
            - (no, 1)  : f, rhp, ihp, rhc, ihc
            - (yes, 1) : f, rhp, ihp, rhc, ihc, hlm, dyn
    """
    par = CreateDict(M=m1+m2, q=m1/m2, chi1=[0,0,s1z], chi2=[0,0,s2z], lambda1=lam1, lambda2=lam2)
    pp = {**par, **additional_pars}
    result = EOB.EOBRunPy(pp)
    if(return_zero):
        del(result)
        return 0
    else:
        return result

def spinsphericalharm(s, l, m, phi, i):
    """
    Compute spin-weighted spherical harmonics
    """
    c = pow(-1.,-s) * np.sqrt( (2.*l+1.)/(4.*np.pi) )
    dWigner = c * wigner_d_function(l,m,-s,i)
    rY = np.cos(m*phi) * dWigner
    iY = np.sin(m*phi) * dWigner
    return rY + 1j*iY

def wigner_d_function(l,m,s,i):
    """
    Compute wigner d functions, following Ref
    TODO: add reference
    """
    costheta = np.cos(i*0.5)
    sintheta = np.sin(i*0.5)
    norm = np.sqrt( (fact(l+m) * fact(l-m) * fact(l+s) * fact(l-s)) )
    ki = max( 0  , m-s )
    kf = min( l+m, l-s )
    dWig = 0.
    for k in range(int(ki), int(kf)+1):
        div = 1.0/( fact(k) * fact(l+m-k) * fact(l-s-k) * fact(s-m+k) )
        dWig = dWig+div*( pow(-1.,k) * pow(costheta,2*l+m-s-2*k) * pow(sintheta,2*k+s-m) )
    return (norm * dWig)

def compute_hphc(hlm, phi=0, i=0, modes=['1']):
    """
    For aligned spins, assuming usual symmetry between hlm and hl-m
    """
    h = 0+1j*0
    for k in modes:
        ki  = int(k)
        ell = k_to_ell(ki)
        emm = k_to_emm(ki)
        Alm = hlm[str(k)][0]
        plm = hlm[str(k)][1]
        Hp  = Alm*np.exp(-1j*plm)
        Hn  = (-1)**ell*Alm*np.exp( 1j*plm)

        Ylmp = spinsphericalharm(-2, ell, emm, np.pi/2-phi, i)
        Ylmn = spinsphericalharm(-2, ell,-emm, np.pi/2-phi, i)
        h   += Ylmp*Hp + Ylmn*Hn

    hp =  np.real(h)
    hc = -np.imag(h)
    return hp, hc

def k_to_ell(k):
    LINDEX = [\
    2,2,\
    3,3,3,\
    4,4,4,4,\
    5,5,5,5,5,\
    6,6,6,6,6,6,\
    7,7,7,7,7,7,7,\
    8,8,8,8,8,8,8,8]
    return LINDEX[k]

def k_to_emm(k):
    MINDEX = [\
    1,2,\
    1,2,3,\
    1,2,3,4,\
    1,2,3,4,5,\
    1,2,3,4,5,6,\
    1,2,3,4,5,6,7,\
    1,2,3,4,5,6,7,8]
    return MINDEX[k]


def modes_to_k(modes):
    """
    Map multipolar (l,m) -> linear index k
    """
    return [int(x[0]*(x[0]-1)/2 + x[1]-2) for x in modes]

if __name__ == '__main__':

    m1, m2, s1z, s2z, l1, l2 = 40., 20., 0.8, 0.5, 0., 0.
    # reference_distance = 100.
    # test_distance      = [200., 1e3, 1e5, 1e9]

    # t,hp_ref,hc_ref = gen_wf(m1, m2, s1z, s2z, l1, l2, additional_pars={'distance':reference_distance}) 
    # A_ref = np.sqrt(hp_ref**2+hc_ref**2)
    # for dis in test_distance:
    #     tt,hp_tst,hc_tst = gen_wf(m1, m2, s1z, s2z, l1, l2, additional_pars={'distance':dis})
    #     A_tst = np.sqrt(hp_tst**2+hc_tst**2)
    #     scale_factor = dis/reference_distance
    #     print(abs(A_ref-A_tst*scale_factor))
    #     #plt.semilogy(tt, abs(A_ref-A_tst*scale_factor))
    #     plt.plot(tt, A_tst*scale_factor)
    # plt.plot(t, A_ref, linestyle='--')
    # plt.show()