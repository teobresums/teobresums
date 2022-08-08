"""
--------------------------------------------------
Test the priors for hyperbolic encounters
---------------------------------------------------
For a given r0, nu we need:
 * j > pphi_lso = EOBRun_module.pph_lso_orbital_py(nu)
 * Emin < E < Emax;
"""

import EOBRun_module
import numpy as np; import matplotlib.pyplot as plt
import time

def D4(f, dx):
    """
    4th order drvt 
    """
    n  = len(f)
    f  = np.array(f)
    df = np.zeros_like(f)
    i  = np.arange(2, n-2)
    oo12dx  = 1./(12*dx)

    df[i] = (8.*(f[i+1]-f[i-1]) - f[i+2] + f[i-2])*oo12dx
    i = 0
    df[i] = (-25.*f[i] + 48.*f[i+1] - 36.*f[i+2] + 16.*f[i+3] - 3.*f[i+4])*oo12dx
    i = 1
    df[i] = (-3.*f[i-1] - 10.*f[i] + 18.*f[i+1] - 6.*f[i+2] + f[i+3])*oo12dx
    i = n-2
    df[i] = - (-3.*f[i+1] - 10.*f[i] + 18.*f[i-1] - 6.*f[i-2] + f[i-3])*oo12dx
    i = n-1
    df[i] = - (-25.*f[i] + 48.*f[i-1] - 36.*f[i-2] + 16.*f[i-3] - 3.*f[i-4])*oo12dx
    
    return df

def E(r, pph, nu):
    """
    E, nonspinning case
    """
    A, dA, d2A  = EOBRun_module.eob_metric_A5PNlog_py(r, nu)
    Heff0       = np.sqrt(A*(1+(pph/r)**2))
    E0          = np.sqrt(1 + 2*nu*(Heff0-1))
    return E0

def Espin(r, pph, q, chi1, chi2):
    """
    E, spinning case
    """
    # New energy potential energy function with the full spin dependence.
    hatH = EOBRun_module.eob_ham_s_py(r, q, pph, 0., chi1, chi2)
    nu   = q/(1+q)**2
    E0   = nu*hatH[0]
    return E0

def EnergyLimits(rmx, nu, q, pph_hyp, chi1, chi2, do_plots=False):
    """
    Determine min and max energy allowed
    """
    if chi1!=0 or chi2!=0:
        # important: rmin should be smaller with spin
        rmin = 1.1
    else:
        rmin = 1.3
        
    x    = np.linspace(rmin,rmx+10, 100000)
    dx   = x[1]-x[0]
    E0   = [Espin(xi, pph_hyp, q, chi1, chi2) for xi in x]
    Emin = Espin(rmx, pph_hyp, q, chi1, chi2)
    
    dE0  = D4(E0,  dx) 
    d2E0 = D4(dE0, dx)

    #-----------------------------------
    #  determine the max energy allowed
    #-----------------------------------

    Emx = np.max(E0)

    #---------------------------------------
    # determine the inflection point of the
    # potential energy E0.
    #--------------------------------------

    jflex = np.where(d2E0 >= 0)[0][0] #FIXME: this needs fixing
    Einfl = E0[jflex]

    if do_plots:
        plt.plot(x, E0)
        plt.xlabel('$r$')
        plt.ylabel('$\hat{E}$')
        plt.hlines(Einfl, x[0], x[-1], colors = 'b',   label='Inflection point')
        plt.hlines(Emx,  x[0], x[-1],  linestyle='--', label=r'$E_{\rm max}$')
        plt.hlines(Emin, x[0], x[-1],  linestyle=':',  label=r'$E_{\rm min}$')
        plt.legend()
        plt.show()

    return Emin, Emx, Einfl

def EobAdm(flag, q_vec, p_vec, nu):
    """
    functions to convert EOB q and p to ADM.
    flag can be 'eob2adm' or 'adm2eob'
    """

    # shorthands
    q2      = np.dot(q_vec, q_vec); # x, y
    q       = np.sqrt(q2);
    q3      = q*q2;
    q4      = q*q3;
    p2      = np.dot(p_vec, p_vec);
    p       = np.sqrt(p2);
    p3      = p*p2;
    p4      = p*p3; 
    qdotp   = np.dot(q_vec, p_vec);
    qdotp2  = qdotp*qdotp;
    nu2     = nu*nu;

    if flag == 'eob2adm':
        # coefficients for coordinates
        cQ_1PN_q = nu*p2/2 - (1 + nu/2)/q;
        cQ_1PN_p = nu*qdotp; 
        cQ_2PN_q = -nu/8*(1 + nu)*p4 + 3/4*nu*(nu/2 - 1)*p2/q - nu*(2 + 5/8*nu)*qdotp2/q3 + (-nu2 + 7*nu - 1)/4/q2;
        cQ_2PN_p = qdotp*(nu*(nu - 1)/2*p2 + nu/2*(-5 + nu/2)/q);
        # coefficients for momenta
        cP_1PN_q = -(1 + nu/2)*qdotp/q3;
        cP_1PN_p = -nu/2*p2 + (1 + nu/2)/q;
        cP_2PN_q = qdotp/q3*(3/4*nu*(nu/2 - 1)*p2 + 3/8*nu2*qdotp2/q2 + (-3/2 + 5/2*nu - 3/4*nu2)/q);
        cP_2PN_p = nu*(1 + 3*nu)/8*p4 - nu/4*(1 + 7/2*nu)*p2/q + nu*(1 + nu/8)*qdotp2/q3 + (5/4 - 3/4*nu + nu2/2)/q2;

    elif flag == 'adm2eob':
        # coefficients for coordinates
        cQ_1PN_q = -nu/2*p2 + 1/q*(1 + nu/2);
        cQ_1PN_p = -qdotp*nu; 
        cQ_2PN_q = nu/8*(1 - nu)*p4 + nu/4*(5 - nu/2)*p2/q + nu*(1 + nu/8)*qdotp2/q3 + 1/4*(1 - 7*nu + nu2)/q2;
        cQ_2PN_p = qdotp*(nu/2*(1 + nu)*p2 + 3/2*nu*(1 - nu/2)/q);
        # coefficients for momenta
        cP_1PN_q = qdotp/q3*(1 + nu/2);
        cP_1PN_p = nu/2*p2 - 1/q*(1 + nu/2);
        cP_2PN_q = qdotp/q3*(nu/8*(10 - nu)*p2 + 3/8*nu*(8 + 3*nu)*qdotp2/q2 + 1/4*(-2 - 18*nu + nu2)/q);
        cP_2PN_p = nu/8*(-1 + 3*nu)*p4 - 3/4*nu*(3 + nu/2)*p2/q  - nu/8*(16 + 5*nu)*qdotp2/q3 + 1/4*(3 + 11*nu)/q2;
    
    # Put all together 
    Q_vec = q_vec + cQ_1PN_q*q_vec + cQ_1PN_p*p_vec + cQ_2PN_q*q_vec + cQ_2PN_p*p_vec; 
    P_vec = p_vec + cP_1PN_q*q_vec + cP_1PN_p*p_vec + cP_2PN_q*q_vec + cP_2PN_p*p_vec;

    return Q_vec, P_vec

if __name__ == "__main__":

    ############################
    # Parameters of the system
    ############################
    q    = 4
    chi1 = 0.0        # spin larger object
    chi2 = 0.0        # spin smaller object
    r    = 50         # initial separation, mass rescaled. r=10000 to study scattering, r=1500 for parspace
    j0   = 4
    E0   = 1.0027     # > 1 for a capture
    nu   = q/(1+q)**2

    do_plots = True   # generate plots
    admeob   = False; # Compute the EOB to ADM transformation to 2PN
    
    ############################
    # Compute energy limits
    ############################
    pphi_lso        = EOBRun_module.pph_lso_orbital_py(nu);
    Emn, Emx, Einfl = EnergyLimits(r, nu, q, j0, chi1, chi2, do_plots=do_plots)
    print("Emin = %s, Emax = %s" %(Emn, Emx))
    print("q    = %s" %q)
    print("r0   = %s" %r)
    print("pph0 = %s" %j0)
    print("E0   = %s" %E0)

    ############################
    # Generate the waveform
    ############################
    pars = {
    'M'                  : 1,
    'q'                  : q,
    'chi1'               : chi1,
    'chi2'               : chi2,
    'LambdaAl2'          : 0.,
    'LambdaBl2'          : 0.,
    'j_hyp'              : j0,                # J_hyp
    'r_hyp'              : r,                 # r_hyp
    'H_hyp'              : E0,                # H_hyp

    'domain'             : 0,                 # Set 1 for FD. Default = 0
    'arg_out'            : "yes",             # Output hlm/hflm. Default = 0
    'nqc'                : "manual",
    'nqc_coefs_hlm'      : "none",
    'nqc_coefs_flx'      : "none",
    'use_mode_lm'        : [1],               # List of modes to use/output through EOBRunPy
    'output_lm'          : [1],               # List of modes to print on file
    'output_dynamics'    : "no",              # Output of the dynamics
    'ode_tstep_opt'      : "adaptive",        # Fixing uniform or adaptive. Default = 1 
    'srate_interp'       : 300000,            # Srate at which to interpolate. Default = 4096.
    'use_geometric_units': "no",              # Output quantities in geometric units. Default = 1
    'interp_uniform_grid': "yes",             # Interpolate mode by mode on a uniform grid. Default = 0 (no interpolation)
    'ode_tmax'           : 20e4,
    'output_hpc'         : "no",              # Output waveform. Default = 1.
    }

    # Run the wf generator
    start = time.time()
    t, hp, hc, hlm, dyn = EOBRun_module.EOBRunPy(pars)    
    end = time.time()
    DeltaT = end-start
    print("ODE time=%s"%DeltaT);
   
    # Unpack the dynamics
    Momg_o = dyn['MOmega_orb']
    T      = dyn['t']
    r      = dyn['r']
    phi    = dyn['phi']
    Pphi   = dyn['Pphi']
    Prstar = dyn['Prstar']

    ############################
    # Plots
    ############################
    if do_plots:

        # plot
        plt.title(r"Polarizations $h_{+}$ and $h_{\times}$")
        plt.plot(t, hp, label=r'$h_+$')
        plt.plot(t, hc, label=r'$h_{\times}$')
        plt.xlabel('$t$ [s]')
        plt.legend()
        plt.show()

        # plot r
        plt.title("Relative EOB distance, $r(t/M)$")
        plt.plot(T,r)
        plt.xlabel("$t/M$")
        plt.ylabel("$r$")
        plt.show()
   
        # plot dynamics
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.set_title("EOB Trajectory")
        plt.polar(phi, r)
        plt.show()

    ############################
    # EOB to ADM
    ############################
    if admeob:

        dashes    = '----------------------------'
        
        # initial coordinate and momenta
        r0      = r[0]
        phi0    = phi[0]
        prstar0 = Prstar[0]
        pph0    = Pphi[0]
        A, B    = EOBRun_module.eob_metricAB_py(r0, q)
        pr0     = prstar0*np.sqrt(B/A)
    
        # cartesian 
        x   = r0
        y   =  0
        py  = pph0/x
        p   = np.sqrt((pph0/r0)**2+pr0**2)
        ph  = np.arcsin(pph0/(x*p))
        px  = -p*np.cos(ph)
        
        qe = np.array([x , y])
        pe = np.array([px, py])

        qa, pa = EobAdm('eob2adm', qe, pe, nu)
        print('\n', dashes, 'EOB', dashes)
        print(' qe  = [',  qe[0], ', ',  qe[1], ']')
        print(' pe  = [',  pe[0], ', ',  pe[1], ']')
        print('\n', dashes, 'ADM', dashes)
        print(' qa  = [',  qa[0], ', ',  qa[1], ']')
        print('|qa| =  ', np.sqrt(np.dot(qa,qa)))
        print(' pa  = [',  pa[0], ', ',  pa[1], ']')
        print('\n', dashes[:-3], 'EOB check', dashes[:-3])
        qe_check, pe_check = EobAdm('adm2eob', qa, pa, nu);
        print(' qe  = [',  qe_check[0], ', ',  qe_check[1], ']')
        print(' pe  = [',  pe_check[0], ', ',  pe_check[1], ']')