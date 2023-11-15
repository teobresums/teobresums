#
# Functions for generating a large number of TEOBResumS waveforms
# over the BBH and BNS parameter spaces, and plotting the results
#
# RG, 12/21


import numpy as np; import EOBRun_module as EOB; 
import matplotlib.pyplot as plt 

def compute_lambda_tilde(m1, m2 ,l1 , l2):
    """ Compute Lambda Tilde from masses and tides components
        --------
        m1 = primary mass component [solar masses]
        m2 = secondary mass component [solar masses]
        l1 = primary tidal component [dimensionless]
        l2 = secondary tidal component [dimensionless]
        """
    M       = m1 + m2
    m1_4    = m1**4.
    m2_4    = m2**4.
    M5      = M**5.
    comb1   = m1 + 12.*m2
    comb2   = m2 + 12.*m1
    return (16./13.)*(comb1*m1_4*l1 + comb2*m2_4*l2)/M5

def compute_delta_lambda(m1, m2 ,l1 , l2):
    """ Compute delta Lambda Tilde from masses and tides components
        --------
        m1 = primary mass component [solar masses]
        m2 = secondary mass component [solar masses]
        l1 = primary tidal component [dimensionless]
        l2 = secondary tidal component [dimensionless]
        """
    M       = m1+m2
    q       = m1/m2
    eta     = q/((1.+q)*(1.+q))
    X       = np.sqrt(1.-4.*eta)
    m1_4    = m1**4.
    m2_4    = m2**4.
    M4      = M**4.
    comb1   = (1690.*eta/1319. - 4843./1319.)*(m1_4*l1 - m2_4*l2)/M4
    comb2   = (6162.*X/1319.)*(m1_4*l1 + m2_4*l2)/M4
    return comb1 + comb2

def compute_chi_prec(m1,m2,s1x,s1y,s2x,s2y):
    """ Compute chi precessing spin parameter (for given 3-dim spin vectors)
        --------
        m1 = primary mass component [solar masses]
        m2 = secondary mass component [solar masses]
        s1 = primary spin megnitude [dimensionless]
        s2 = secondary spin megnitude [dimensionless]
        tilt1 = primary spin tilt [rad]
        tilt2 = secondary spin tilt [rad]
    """

    s1_perp = np.sqrt(s1x**2+s1y**2)
    s2_perp = np.sqrt(s2x**2+s2y**2)
    one_q   = m2/m1

    # check that m1>=m2, otherwise switch
    if one_q > 1. :
        one_q = 1./one_q
        s1_perp, s2_perp = s2_perp, s1_perp

    return np.max([s1_perp , s2_perp*one_q*(4.*one_q+3.)/(3.*one_q+4.)])


def Plot2D(x, y, z=None, labels=[None,None], savef=0):
    """
    Make a 2D plot of x vs y
    """

    print("...plot "+labels[0]+'-'+labels[1])
    fig = plt.figure()
    plt.scatter(x, y, c=z, marker='.')
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])

    if(savef):
        fig.savefig(labels[0]+'_'+labels[1]+'.png')

def PlotParspace(f, tides="no",precessing="no"):
    """
    Plot the output of TestParBBH or TestParBNS
    """
    pars = np.genfromtxt(f, names=True)
    q  = pars['q']
    M  = pars['M']
    
    nu = q/(1+q)**2 
    m1 = M*q/(1.+q)
    m2 = M*1/(1.+q)
    Plot2D(nu, M, labels=[r'$\nu$', r'$M$'], savef=1)

    if(tides=='yes'):
        try:
            lam1 = pars['LambdaAl2']
        except ValueError:
            lam1 = 0.
        lam2 = pars['LambdaBl2']
        lamt = compute_lambda_tilde(m1, m2, lam1, lam2)
        dlam = compute_delta_lambda(m1, m2, lam1, lam2)
        Plot2D(nu, lamt, labels=[r"$\nu$", r"$\tilde\Lambda$"], savef=1)

    if(precessing=='yes'):
        c1x = pars['chi1x'];  c2x = pars['chi2x']
        c1y = pars['chi1y'];  c2y = pars['chi2y']
        cp  = list(map(compute_chi_prec, m1, m2, c1x,c1y, c2x, c2y))
        Plot2D(nu, cp, labels=[r"$\nu$", r"$\chi_{\rm prec}$"], savef=1)

def write_dict_to_txt(file, dic, initialize):
    """
    Write a dictionary to a .txt
    """
    with open(file, "a") as f:
        if(initialize):
            f.write("# ")
            for k in dic.keys():
                f.write(k)
                f.write(" ")
            f.write("\n")
        for k in dic.keys():
            f.write(str(dic[k]))
            f.write(" ")
        f.write("\n")

def modes_to_k(modes):
    """
    Map multipolar (l,m) -> linear index k
    """
    return [int(x[0]*(x[0]-1)/2 + x[1]-2) for x in modes]

def gen_random_pars(Mint, qint, chi_int=None, lambda_int=None, precessing=0, N=100):
    """
    Generate randomly N combinations of binary parameters
    """
    M    = np.random.uniform(Mint[0], Mint[1], N)
    q    = np.random.uniform(qint[0], qint[1], N)
    
    params = {'M': M, 'q':q}

    if(chi_int is not None):
        chi1 = np.random.uniform(chi_int[0], chi_int[1], N)
        chi2 = np.random.uniform(chi_int[0], chi_int[1], N)

        params['chi1']=chi1; params['chi2'] = chi2

    if(lambda_int is not None):
        lambda1 =  np.random.uniform(lambda_int[0], lambda_int[1], N)
        lambda2 =  np.random.uniform(lambda_int[0], lambda_int[1], N)

        params['lambda1'] = lambda1; params['lambda2'] = lambda2

    if(precessing):
        cos_theta1 = np.random.uniform(-1., 1., N)
        cos_theta2 = np.random.uniform(-1., 1., N)
        phi1 = np.random.uniform(0., 2*np.pi, N)
        phi2 = np.random.uniform(0., 2*np.pi, N)

        sin_theta1 = np.sqrt(1-cos_theta1**2)
        sin_theta2 = np.sqrt(1-cos_theta2**2)

        chi1z = cos_theta1*chi1
        chi1x = sin_theta1*np.cos(phi1)*chi1
        chi1y = sin_theta1*np.sin(phi1)*chi1

        chi2z = cos_theta2*chi2
        chi2x = sin_theta2*np.cos(phi2)*chi2
        chi2y = sin_theta2*np.sin(phi2)*chi2

        params['chi1x'] = chi1x; params['chi2x'] = chi2x
        params['chi1y'] = chi1y; params['chi2y'] = chi2y
        params['chi1z'] = chi1z; params['chi2z'] = chi2z
    else:
        params['chi1x'] = np.zeros(N); params['chi2x'] = np.zeros(N);
        params['chi1y'] = np.zeros(N); params['chi2y'] = np.zeros(N);
        params['chi1z'] = chi1; params['chi2z'] = chi2

    return params

def CreateDict(M, q, chi1, chi2, l1, l2, iota, f0, srate, df, interp, domain, modes, coa=0, argout="no"):
    """
    Create the dictionary of parameters for EOBRunPy
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
    'LambdaAl2'          : l1,
    'LambdaBl2'          : l2,
    'distance'           : 1.,
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
    'output_hpc'         : "no",
    'output_multipoles'  : "no", 
    'arg_out'            : argout,
    'output_dynamics'    : "no",
    'project_spins'      : "no",
    'time_shift_FD'      : "no",
    }
    return pardic

def TestParspaceBBH(precessing):
    """
    Generate 1e4 precessing waveforms with parameters within standard BBH bounds
    """
    modes = [[2,1], [2,2], [3,1], [3,2], [3,3], [4,1], [4,2], [4,3], [4,4]]

    k = modes_to_k(modes)

    Mmin, Mmax = 20, 300
    qmin, qmax = 1, 10

    chi_min    = 1e-3
    chi_max    = 0.99

    N = 10000 # generate 1e4 waveforms

    f0    = 11. # initial frequency
    srate = 4096.

    ### DO NOT CHANGE BELOW UNLESS YOU KNOW WHAT YOU ARE DOING ###

    print("...generate the parspace")
    rand_pars = gen_random_pars([Mmin, Mmax], [qmin,qmax], chi_int=[chi_min, chi_max], lambda_int=None, precessing=precessing, N=N)

    # start generating waveforms

    init_ok  = 1
    init_err = 1
    print("...run")
    for m, qi, c1x, c1y, c1z, c2x, c2y, c2z in zip(rand_pars['M'], rand_pars['q'], rand_pars['chi1x'], rand_pars['chi1y'], rand_pars['chi1y'], rand_pars['chi2x'], rand_pars['chi2y'], rand_pars['chi2z']):
        
        par = CreateDict(m, qi, [c1x, c1y, c1z], [c2x, c2y, c2z], 0., .0, np.pi/3, f0, srate, 0., "yes", 0, k)
        try:
            t, hp, hc = EOB.EOBRunPy(par)
            # remove keys that must not be written (this is a dumb way..)
            for lbl in ['distance','initial_frequency','use_geometric_units','interp_uniform_grid','domain','srate_interp','inclination','use_mode_lm','output_lm','output_hpc','output_multipoles','arg_out','output_dynamics','use_spins']:
                par.pop(lbl)
            write_dict_to_txt("ParBBH.txt", par, init_ok)
            if(init_ok):
                init_ok=0
    
        except Exception:
            write_dict_to_txt("ErrorsParBBH.txt", par, init_err)
            if(init_err):
                init_err=0

def TestParspaceBNS(precessing):
    """
    Generate 1e4 precessing waveforms with parameters within standard BNS bounds
    """

    modes = [[2,2]]

    k = modes_to_k(modes)

    Mmin, Mmax = 1.5, 3.4
    qmin, qmax = 1, 2

    lambda_min = 3
    lambda_max = 5000

    chi_min    = 1e-3
    chi_max    = 0.6

    N = 1000 # generate 1e4 waveforms

    f0    = 20. # initial frequency
    srate = 2048

    ### DO NOT CHANGE BELOW UNLESS YOU KNOW WHAT YOU ARE DOING ###

    print("...generate the parspace")
    rand_pars = gen_random_pars([Mmin, Mmax], [qmin,qmax], chi_int=[chi_min, chi_max], lambda_int=[lambda_min, lambda_max], precessing=precessing, N=N)

    # start generating waveforms

    init_ok  = 1
    init_err = 1
    print("...run")
    for m, qi, c1x, c1y, c1z, c2x, c2y, c2z, l1, l2 in zip(rand_pars['M'], rand_pars['q'], rand_pars['chi1x'], rand_pars['chi1y'], rand_pars['chi1y'], rand_pars['chi2x'], rand_pars['chi2y'], rand_pars['chi2z'], rand_pars['lambda1'], rand_pars['lambda2']):
        
        par = CreateDict(m, qi, [c1x, c1y, c1z], [c2x, c2y, c2z], l1, l2, np.pi/3, f0, srate, 0., "yes", 0, k)
        try:
            t, hp, hc = EOB.EOBRunPy(par)

            # remove keys that must not be written (this is a dumb way..)
            for lbl in ['distance','initial_frequency','use_geometric_units','interp_uniform_grid','domain','srate_interp','inclination','use_mode_lm','output_lm','output_hpc','output_multipoles','arg_out','output_dynamics','use_spins']:
                par.pop(lbl)
            write_dict_to_txt("ParBNS.txt", par, init_ok)
            if(init_ok):
                init_ok=0
    
        except Exception:
            write_dict_to_txt("ErrorsParBNS.txt", par, init_err)
            if(init_err):
                init_err=0

def TestParspaceBHNS(precessing):
    """
    Generate 1e4 precessing waveforms with parameters within standard BHNS bounds
    """

    modes = [[2,1], [2,2], [3,2], [3,3], [4,4]] #[[2,2]]

    k = modes_to_k(modes)

    Mmin, Mmax = 3, 40
    qmin, qmax = 1, 10

    lambda_min = 3
    lambda_max = 5000

    chi_min    = 1e-3
    chi_max    = 0.99

    N = 10000 # generate 1e4 waveforms

    f0    = 20. # initial frequency
    srate = 2048

    ### DO NOT CHANGE BELOW UNLESS YOU KNOW WHAT YOU ARE DOING ###

    print("...generate the parspace")
    rand_pars = gen_random_pars([Mmin, Mmax], [qmin,qmax], chi_int=[chi_min, chi_max], lambda_int=[lambda_min, lambda_max], precessing=precessing, N=N)

    # start generating waveforms

    init_ok  = 1
    init_err = 1
    print("...run")
    for m, qi, c1x, c1y, c1z, c2x, c2y, c2z, l2 in zip(rand_pars['M'], rand_pars['q'], rand_pars['chi1x'], rand_pars['chi1y'], rand_pars['chi1y'], rand_pars['chi2x'], rand_pars['chi2y'], rand_pars['chi2z'], rand_pars['lambda2']):
        
        par = CreateDict(m, qi, [c1x, c1y, c1z], [c2x, c2y, c2z], 0, l2, np.pi/3, f0, srate, 0., "yes", 0, k)
        
        try:
            t, hp, hc = EOB.EOBRunPy(par)

            # remove keys that must not be written (this is a dumb way..)
            for lbl in ['LambdaAl2','distance','initial_frequency','use_geometric_units','interp_uniform_grid','domain','srate_interp','inclination','use_mode_lm','output_lm','output_hpc','output_multipoles','arg_out','output_dynamics','use_spins']:
                par.pop(lbl)
            write_dict_to_txt("ParPrecBHNS_HM.txt", par, init_ok)
            if(init_ok):
                init_ok=0
    
        except Exception:
            write_dict_to_txt("ErrorsParBHNS.txt", par, init_err)
            if(init_err):
                init_err=0

def PlotPrecWF():
    """
    Plot one example BBH precessing waveform
    """
    M = 250
    q = 1.1
    chiA = [0.8, 0., 0.]
    chiB = [0., -0.7, -0.1]
    par = CreateDict(M, q, chiA, chiB, 0, 0, np.pi/3, 11., 4096., 0., "yes", 0, [0,1,2,3,4,5,6,7,8])

    t, hp, hc = EOB.EOBRunPy(par)
    
    fig = plt.figure()
    plt.plot(t, hp)
    plt.plot(t, hc)
    plt.xlabel('t')
    fig.savefig('h_test.png')
    
def PlotAlignedSpinLimit():
    """
    Plot a precessing wf in the aligned spins limit and compare it with an aligned spins one
    """
    M = 60
    q = 6
    
    # above the 1e-4 threshold 
    chiA = [0.0002, 0.,   0.]
    chiB = [0.000,  0., -0.1]

    # blow the 1e-4 threshold
    chiAt = [0.000001, 0., 0.0]
    chiBt = [0.000,    0.,-0.1]

    # aligned
    chiAa = [0., 0.,  0.0]
    chiBa = [0., 0., -0.1]

    par_p = CreateDict(M, q, chiA,  chiB,  0, 0, np.pi/3, 20., 4096., 0., "yes", 0, [0,1,2,3,4,5,6,7,8,13], argout="no")
    par_t = CreateDict(M, q, chiAt, chiBt, 0, 0, np.pi/3, 20., 4096., 0., "yes", 0, [0,1,2,3,4,5,6,7,8,13], argout="no")    
    par_a = CreateDict(M, q, chiAa, chiBa, 0, 0, np.pi/3, 20., 4096., 0., "yes", 0, [0,1,2,3,4,5,6,7,8,13], argout="no")

    fig = plt.figure()
    for pp in [par_a, par_p, par_t]:
        t,hp,hc = EOB.EOBRunPy(pp)
        plt.plot(t, hc)

    plt.xlabel('t')
    plt.show()

def PlotSPAWF():
    """
    Plot one example BNS+SPA precessing waveform, compare with the TD
    """
    from pycbc.types.frequencyseries import FrequencySeries

    M = 2.7
    q = 1
    nu= q/(1+q)**2
    chiA  = [0.0, 0., 0.5]
    chiB  = [0.0, 0., 0.1]
    f0    = 20.
    df    = 1./(128*2)
    # parameter dictionaries
    part     = CreateDict(M, q, chiA, chiB, 400, 400, np.pi/2.5, f0, 4096., 0., "yes",  0, [1], coa=np.pi/4, argout="yes")
    parf     = CreateDict(M, q, chiA, chiB, 400, 400, np.pi/2.5, f0, 4096., df, "yes" , 1, [1], coa=np.pi/4, argout="yes")
    
    # Gen TD wf
    t,hpt,hct,hlm,_= EOB.EOBRunPy(part)
    # Gen FD wf
    f_spa, rhpf, ihpf, rhcf,ihcf,_,_,_= EOB.EOBRunPy(parf)
    hf = rhpf - 1j*ihpf
    hfc= rhcf - 1j*ihcf
    
    # ifft
    N = int(f0/(f_spa[1]-f_spa[0]))
    hf_spa_tmp  = np.pad(hf,  (N,0), mode='constant')
    hfc_spa_tmp = np.pad(hfc, (N,0), mode='constant')
    hf_spa_FD   = FrequencySeries(hf_spa_tmp, delta_f=f_spa[1]-f_spa[0])
    hfc_spa_FD  = FrequencySeries(hfc_spa_tmp, delta_f=f_spa[1]-f_spa[0])
    hf_spa_TD   = hf_spa_FD.to_timeseries()
    hfc_spa_TD  = hfc_spa_FD.to_timeseries()
    # plot
    plt.plot(t, np.sqrt(hpt**2+hct**2), color='k')
    plt.plot(hf_spa_TD.sample_times, np.sqrt(hf_spa_TD**2+hfc_spa_TD**2), linestyle='--', color='r', label=r'TEOBResumSPA')    
    plt.xlabel('t [s]')
    plt.ylabel('$A$')
    plt.legend()
    plt.show()

def SwapPrecWF():
    """
    Test the swapping of components when q<1.
    --------------------------------------------------
    According to LAL conventions, if m1<->m2 the x axis is flipped.
    Therefore, the waveform has to remain identical when:
        - label_1<->label_2;
        - coalescence_angle -> coalescence_angle + Pi;
        - In plane spins are rotated by Pi
    """
    M    = 100
    q    = 2.
    chiA     = [ 0.8,  0.,  0. ]
    chiA_rot = [-0.8,  0.,  0. ]
    chiB     = [  0.,-0.7, -0.1]
    chiB_rot = [  0., 0.7, -0.1]
    phi_ref  = np.pi/4

    par     = CreateDict(M, q,    chiA,     chiB,     0, 0, np.pi/3, 20., 4096., 0, "yes", 0, [0,1], coa=phi_ref, argout="yes")
    par_rot = CreateDict(M, 1./q, chiB_rot, chiA_rot, 0, 0, np.pi/3, 20., 4096., 0, "yes", 0, [0,1], coa=phi_ref+np.pi, argout="yes")

    t1, hp1,_,_,_ = EOB.EOBRunPy(par)
    t2, hp2,_,_,_ = EOB.EOBRunPy(par_rot)

    plt.plot(t1, hp1, color='r', label=r'$m_1 > m_2$', linewidth=2)
    plt.plot(t2, hp2, color='royalblue', linestyle='--', linewidth=2, label=r'$m_2 > m_1$')
    plt.grid()
    plt.xlabel(r'$t$ [s]')
    plt.ylabel(r'$h_+$')
    plt.legend()
    plt.show()

if __name__ == "__main__":

    # run some tests

    if 0:
        print("##### Generate the BBH precessing parameter space #####")
        TestParspaceBBH(1)
        print("...done")

    if 0:
        print("##### Generate the BNS precessing parameter space #####")
        TestParspaceBNS(1)
        print("...done")

    if 0:
        print("##### Generate the BHNS precessing parameter space #####")
        TestParspaceBHNS(1)
        print("...done")

    if 0:
        print("##### Plot the BBH precessing parameter space #####")
        PlotParspace("ParPrecBBH.txt", precessing="yes")
        print("...done")

    if 0:
        print("##### Plot the BNS precessing parameter space #####")
        PlotParspace("ParPrecBNS.txt", tides="yes", precessing="yes")
        print("...done")
    
    if 0:
        print("##### Plot the BHNS precessing parameter space #####")
        PlotParspace("ParPrecBHNS.txt", tides="yes", precessing="yes")
        print("...done")

    if 0:
        print('##### Plot Prec WF #####')
        PlotPrecWF()
        print("...done")

    if 1:
        print('##### Test the spin aligned limit')
        PlotAlignedSpinLimit()
        print('...done')

    if 1:
        print('##### Plot a FD SPA WF, compare to TD #####')
        PlotSPAWF()
        print("...done")

    if 1:
        print("##### Test Swap #####")
        SwapPrecWF()
        print("...done")