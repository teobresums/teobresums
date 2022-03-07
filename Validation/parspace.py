#
# Functions for generating a large number of TEOBResumS waveforms
# over the BBH and BNS parameter spaces, and plotting the results
#
# RG, 12/21


import numpy as np; import EOBRun_module as EOB; 

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
    import matplotlib.pyplot as plt

    print("...plot "+labels[0]+'-'+labels[1])
    fig = plt.figure()
    plt.scatter(x, y, c=z)
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
        lam1 = pars['Lambda1']
        lam2 = pars['Lambda2']
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

def CreateDict(M, q, chi1, chi2, l1, l2, iota, f0, srate, df, interp, domain, modes):
    """
    Create the dictionary of parameters for EOBRunPy
    """

    chi1x, chi1y, chi1z = chi1
    chi2x, chi2y, chi2z = chi2

    pardic = {
    'M'                  : M,
    'q'                  : q,
    'chi1'               : chi1z,
    'chi2'               : chi2z,
    'chi1x'              : chi1x,
    'chi1y'              : chi1y,
    'chi1z'              : chi1z,
    'chi2x'              : chi2x,
    'chi2y'              : chi2y,
    'chi2z'              : chi2z,
    'Lambda1'            : l1,
    'Lambda2'            : l2,
    'distance'           : 1.,
    'initial_frequency'  : f0,
    'use_geometric_units': "no",
    'interp_uniform_grid': interp,
    'domain'             : domain,
    'srate_interp'       : srate,
    'inclination'        : iota,
    'use_mode_lm'        : modes,      # List of modes to use/output through EOBRunPy
    'output_lm'          : modes,
    'output_hpc'         : "no",
    'output_multipoles'  : "no", 
    'arg_out'            : 0,
    'output_dynamics'    : "no",
    'use_spins'          : 2,
    'project_spins'      : "no",
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
    
    import matplotlib.pyplot as plt 
    fig = plt.figure()
    plt.plot(t, hp)
    plt.plot(t, hc)
    plt.xlabel('t')
    fig.savefig('h_test.png')
    

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
        print("##### Plot the BBH precessing parameter space #####")
        PlotParspace("ParPrecBBH.txt", precessing="yes")
        print("...done")

    if 0:
        print("##### Plot the BNS precessing parameter space #####")
        PlotParspace("ParPrecBNS.txt", tides="yes", precessing="yes")
        print("...done")
    if 0:
        print('##### Plot Prec WF #####')
        PlotPrecWF()
        print("...done")
