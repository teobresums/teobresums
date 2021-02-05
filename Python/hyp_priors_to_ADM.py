import EOBRun_module
import numpy as np
import math as math
import cmath as cmath
import matplotlib.pyplot as plt
import time

# --------------------------------------------------
# file to test the priors for hyperbolic encounters
#---------------------------------------------------
# for a given r0, nu we need:
# * j > pphi_lso = EOBRun_module.pph_lso_orbital_py(nu);
# * Emin < E < Emax;

def D4(f, dx):
    n = len(f)
    df = np.zeros_like(f)
    oo12dx  = 1./(12*dx);
    
    for i in range(2, n-2):
        df[i] = (8.*(f[i+1]-f[i-1]) - f[i+2] + f[i-2])*oo12dx;

    i = 0;
    df[i] = (-25.*f[i] + 48.*f[i+1] - 36.*f[i+2] + 16.*f[i+3] - 3.*f[i+4])*oo12dx;
    i = 1;
    df[i] = (-3.*f[i-1] - 10.*f[i] + 18.*f[i+1] - 6.*f[i+2] + f[i+3])*oo12dx;
    i = n-2;
    df[i] = - (-3.*f[i+1] - 10.*f[i] + 18.*f[i-1] - 6.*f[i-2] + f[i-3])*oo12dx;
    i = n-1;
    df[i] = - (-25.*f[i] + 48.*f[i-1] - 36.*f[i-2] + 16.*f[i-3] - 3.*f[i-4])*oo12dx;
    
    return df;

def E(r, pph, nu):
    A, dA, d2A  = EOBRun_module.eob_metric_A5PNlog_py(r, nu)
    Heff0       = np.sqrt(A*(1+(pph/r)**2))
    E0          = np.sqrt(1 + 2*nu*(Heff0-1))
    return E0

def Espin(r, pph, q, chi1, chi2):
    # New energy potential energy function with the full spin dependence.
    hatH = EOBRun_module.eob_ham_s_py(r, q, pph, 0., chi1, chi2)
    nu   = q/(1+q)**2
    E0   = nu*hatH[0]
    return E0


def EnergyLimits(rmx, nu, q, pph_hyp, chi1, chi2):

    if chi1!=0 or chi2!=0:
        # important: rmin should be smaller with spin
        rmin = 1.1
    else:
        rmin = 1.3
        
    x    = np.linspace(rmin,rmx+10, 100000)
    dx   = x[1]-x[0]
    E0   = [Espin(xi, pph_hyp, q, chi1, chi2) for xi in x]
    Emin = Espin(rmx, pph_hyp, q, chi1, chi2)
    
    #E0   = [E(xi,  pph_hyp, nu) for xi in x]
    #Emin = E(rmx, pph_hyp, nu)
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

    plt.plot(x, E0)
    plt.xlabel('r')
    plt.ylabel('E')
    plt.hlines(Einfl, x[0], x[-1], colors = 'b')
    plt.hlines(Emx,  x[0], x[-1])
    plt.hlines(Emin, x[0], x[-1])
    plt.show()

    return Emin, Emx, Einfl

def EobAdm(flag, q_vec, p_vec, nu):
# functions to convert EOB q and p to ADM.
# flag can be 'eob2adm' or 'adm2eob'

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

    # main
    # r = 10000 for scattering angle. For parspace survey, r=1500 is ok
    r    = 1500. 
    chi1 =  0.0       # spin larger object
    chi2 =  0.0       # spin smaller object
    #============================
    # FIX E0 > 1
    #============================= 
    q    =  1
    E0   =  1.0055    # initial energy
    j    =  3.97      # initial angular momentum
    q    = 1.498
    E0   = 1.0076
    j    = 3.95
    q    = 1.5
    j    = 4.19
    E0   = 1.012
    j    = 4.22
    q    = 1.6264
    # Glitch-like signal
    #q    = 1
    #E0   = 1.010
    #j    = 4.2
    #q   =  1.0
    #E0  = 1.01
    #j   = 4.0
    nu  = q/(1+q)**2
    
    # nonspinning LSO computation
    pphi_lso = EOBRun_module.pph_lso_orbital_py(nu);
    #j = 1.15*pphi_lso  # initial angular momentumn as multiple of LSO

    # energy limits
    Emn, Emx, Einfl = EnergyLimits(r, nu, q, j, chi1, chi2)
    print("j = %s" %j)
    print("Emin = %s, Emax = %s" %(Emn, Emx))
    print("q    = %s" %q)
    print("r0   = %s" %r)
    print("pph0 = %s" %j)
    print("E0   = %s" %E0)
    
    # define the input parameters
    pars = {
    'M'                  : 1,
    'q'                  : q,
    'chi1'               : chi1,
    'chi2'               : chi2,
    'Lambda1'            : 0.,
    'Lambda2'            : 0.,
    'dt'                 : 0.5,
    'dt_interp'          : 0.5,
    'domain'             : 0,                 #Set 1 for FD. Default = 0
    'arg_out'            : 1,                 #Output hlm/hflm. Default = 0
    'nqc'                : 2,
    'nqc_coefs_hlm'      : 0,
    'nqc_coefs_flx'      : 0,
    'use_mode_lm'        : [1],               #List of modes to use/output through EOBRunPy
    'output_lm'          : [1],               #List of modes to print on file
    'output_dynamics'    : 0,                 #output of the dynamics
    'ode_tstep_opt'      : 1,                 #fixing uniform or adaptive. Default = 1 
    'srate_interp'       : 300000.,           #srate at which to interpolate. Default = 4096.
    'use_geometric_units': 0,                 #output quantities in geometric units. Default = 1
    'r0'                 : r,
    'interp_uniform_grid': 1,                 #interpolate mode by mode on a uniform grid. Default = 0 (no interpolation)
    'ecc'                : 0.18,              #Eccentricity. Default = 0.
    'j_hyp'              : j,                 #J_hyp. Default = 0.
    'r_hyp'              : r,                 #r_hyp. Default = 0.
    'H_hyp'              : E0,                #H_hyp. Default = 0.
    'ode_tmax'           : 20e4,
    'output_hpc'         : 0,                 #output waveform. Default = 1.
    }

    #run the wf generator
    start = time.time()
    t, hp, hc, hlm, dyn = EOBRun_module.EOBRunPy(pars)    
    Momg_o = dyn['MOmega_orb']
    T      = dyn['t']
    r      = dyn['r']
    phi    = dyn['phi']
    Pphi   = dyn['Pphi']
    Prstar = dyn['Prstar']
    #t, hp, hc = EOBRun_module.EOBRunPy(pars)
    end = time.time()
    DeltaT = end-start
    print("ODE time=%s"%DeltaT);
   
    #---------------------------------------------------------------------------------
    do_plots = True;
    debug    = False;

    if do_plots:
        # plot
        plt.plot(t, hp, label=r'h_+')
        plt.plot(t, hc, label=r'h_x')
        plt.legend()
        plt.show()

        # plot Omg
        plt.plot(T,r)
        plt.show()
    
        # plot amplitude and phase
        A22   = hlm['1'][0]
        Phi22 = hlm['1'][1]
        Reh22   = -A22*np.cos(-Phi22)
        Imh22   = -A22*np.sin(-Phi22)

        # plot amplitude and phase
        #A21   = hlm['0'][0]
        #Phi21 = hlm['0'][1]
        #Reh21   = A21*np.cos(-Phi21)
        #Imh21   = A21*np.sin(-Phi21)

        #print("A22 = %s" %A22)
        #print("phi22 = %s" %Phi22)
        #print("A21 = %s" %A21)
        #print("phi21 = %s" %Phi21)

        plt.plot(t,-Reh22, label=r'h22')
        #plt.plot(t,Imh22)
        plt.legend()
        plt.show()

        plt.plot(t,A22)
        plt.show()

        plt.plot(Phi22)
        plt.show()


    dashes    = '----------------------------';
    # initial coordinate and momenta
    r0      = r[0];
    phi0    = phi[0];
    prstar0 = Prstar[0];
    pph0    = Pphi[0];
    A, B    = EOBRun_module.eob_metricAB_py(r0, q);
    pr0     = prstar0*np.sqrt(B/A);
   
    # cartesian 
    x   = r0;
    y   =  0;
    py  = pph0/x;
    p   = np.sqrt((pph0/r0)**2+pr0**2);
    ph  = np.arcsin(pph0/(x*p));
    px  = -p*np.cos(ph);
    
    qe = np.array([x , y]);
    pe = np.array([px, py]);

    qa, pa = EobAdm('eob2adm', qe, pe, nu);
    print('\n', dashes, 'EOB', dashes)
    print('qe = [',  qe[0], ', ',  qe[1], ']')
    print('pe = [',  pe[0], ', ',  pe[1], ']')
    print('\n', dashes, 'ADM', dashes)
    print('qa = [',  qa[0], ', ',  qa[1], ']')
    print('pa = [',  pa[0], ', ',  pa[1], ']')
    print('\n', dashes, 'EOB2', dashes)
    qe, pe = EobAdm('adm2eob', qa, pa, nu);
    print('qe = [',  qe[0], ', ',  qe[1], ']')
    print('pe = [',  pe[0], ', ',  pe[1], ']')

    if debug:
        A_test, dA_test, d2A_test = EOBRun_module.eob_metric_A5PNlog_py(r[0], nu);
        u0 = 1/r0;
        u02 = u0*u0;
        u03 = u0*u02;
        Dp      = 1.0 + 6.*nu*u02 - 2.*(3.0*nu-26.0)*nu*u03; # Pade resummation of D
        D       = 1./Dp;
        B_test  = D/A;
        print('----------A------------')
        print(A)
        print(A_test)
        print('----------B------------')
        print(B)
        print(B_test)
       

