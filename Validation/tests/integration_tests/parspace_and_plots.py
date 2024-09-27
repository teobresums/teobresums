#
# Functions for generating a large number of TEOBResumS waveforms
# over the BBH and BNS parameter spaces, and plotting the results
#
# RG, 12/21

import numpy as np; import EOBRun_module as EOB; 
import matplotlib.pyplot as plt 
import utilities as utils

def eob_dyn_omg_from_omgbar(omgbar, ecc, zeta):
    """
    Compute the instanteous omega from omega_bar (orb avg)
    """
    fact1 = (ecc*np.cos(zeta)+1.)
    den   = pow((1.-ecc*ecc), 3./2)
    return omgbar*fact1*fact1/den

def Espin(r, pph, pr, q, chi1, chi2):
    """
    E, spinning case
    """
    # New energy potential energy function with the full spin dependence.
    hatH = EOB.eob_ham_s_py(r, q, pph, pr, chi1, chi2)
    return hatH

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

def PlotParspace(f, tides="no",precessing="no", eccentric="no", hyperbolic="no"):
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

    if(hyperbolic=='yes'):
        E0 = pars['E0']
        L0 = pars['L0']
        Plot2D(L0, E0, labels=[r'$L_0$', r'$E_0$'], savef=1)

    if(eccentric=='yes'):
        ecc = pars['ecc']
        Plot2D(nu, ecc, labels=[r'$\nu$', r'$e_0$'], savef=1)

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
        else:
            for k in dic.keys():
                f.write(str(dic[k]))
                f.write(" ")
            f.write("\n")

def gen_random_pars(Mint, qint, chi_int=None, lambda_int=None, ecc_int=None, omg_int=None, E0_int=None, L0_int=None, anomaly_int=None, precessing=0, N=100):
    """
    Generate randomly N combinations of binary parameters
    """
    M    = np.random.uniform(Mint[0], Mint[1], N)
    q    = np.random.uniform(qint[0], qint[1], N)
    m1   = q*M/(1. + q)
    m2   =   M/(1.+ q)
    params = {'M': M, 'q':q, 'm1':m1, 'm2':m2 }

    if(chi_int is not None):
        chi1 = np.random.uniform(chi_int[0], chi_int[1], N)
        chi2 = np.random.uniform(chi_int[0], chi_int[1], N)

        params['chi1']=chi1; params['chi2'] = chi2

    if(lambda_int is not None):
        lambda1 =  np.random.uniform(lambda_int[0], lambda_int[1], N)
        lambda2 =  np.random.uniform(lambda_int[0], lambda_int[1], N)

        params['LambdaAl2'] = lambda1; params['LambdaBl2'] = lambda2

    if(ecc_int is not None):
        ecc = np.random.uniform(ecc_int[0], ecc_int[1], N)

        params['ecc'] = ecc

    if(omg_int is not None):
        omg = np.random.uniform(omg_int[0], omg_int[1], N)

        params['omg'] = omg
    
    if(E0_int is not None):
        E0 = np.random.uniform(E0_int[0], E0_int[1], N)

        params['E0'] = E0

    if(L0_int is not None):
        L0 = np.random.uniform(L0_int[0], L0_int[1], N)

        params['L0'] = L0

    if(anomaly_int is not None):
        anomaly = np.random.uniform(anomaly_int[0], anomaly_int[1], N)

        params['anomaly'] = anomaly

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

def TestParspaceBBH(precessing, eccentric):
    """
    Generate 1e4 precessing waveforms with parameters within standard BBH bounds
    """
    import multiprocessing; import tqdm

    modes = [[2,1], [2,2], [3,1], [3,2], [3,3], [4,1], [4,2], [4,3], [4,4]]
    k = utils.modes_to_k(modes)

    Mmin, Mmax = 20, 150
    qmin, qmax = 1, 10
    chi_min    = 1e-3
    chi_max    = 0.99

    N = 10000 # generate 1e4 waveforms

    if eccentric:
        if precessing:
            outfile = "ParBBH_eccprec.txt"
        else:
            outfile = "ParBBH_ecc.txt"
    elif precessing:
        outfile = "ParBBH_prec.txt"
    else:
        outfile = "ParBBH.txt"

    if(eccentric):
        ecc_int = [1.e-5, 0.7]
        anomaly_int = [0., 2*np.pi]
    else:
        ecc_int = None
        anomaly_int = None

    print("...generate the parspace")
    pp    = gen_random_pars([Mmin, Mmax], [qmin,qmax], chi_int=[chi_min, chi_max], lambda_int=None, ecc_int=ecc_int, anomaly_int=anomaly_int, precessing=precessing, N=N)
    write_dict_to_txt(outfile, pp, 1) # write keys to file

    print("...run")
    if(eccentric):
        with multiprocessing.Pool() as p:
            for j in tqdm.tqdm(
                                [[i, p.apply_async(utils.gen_wf, (pp['m1'][i],pp['m2'][i], pp['chi1z'][i],pp['chi2z'][i], 0., 0.,{'ecc':pp['ecc'][i], 'use_mode_lm':k, 'anomaly': pp['anomaly'][i]}))]
                                               for i in range(N)]
                                ):
                j[1].get()
                this_par = {}
                for key in pp.keys():
                    this_par[key] = pp[key][j[0]]
                write_dict_to_txt(outfile, this_par, 0)
    else:
        with multiprocessing.Pool() as p:
            for j in tqdm.tqdm(
                                [[i, p.apply_async(utils.gen_wf, (pp['m1'][i],pp['m2'][i], pp['chi1z'][i],pp['chi2z'][i], 0., 0.,{'use_mode_lm':k}))]
                                               for i in range(N)]
                                ):
                j[1].get()
                this_par = {}
                for key in pp.keys():
                    this_par[key] = pp[key][j[0]]
                write_dict_to_txt(outfile, this_par, 0)

def TestParspaceHyp(precessing):
    """
    Generate 1e4 (precessing) hyperbolic waveforms with parameters within standard BBH bounds
    """
    import multiprocessing; import tqdm

    modes = [[2,1], [2,2], [3,1], [3,2], [3,3], [4,1], [4,2], [4,3], [4,4]]
    k = utils.modes_to_k(modes)

    Mmin, Mmax = 20, 100
    qmin, qmax = 1, 10
    chi_min    = -0.99
    chi_max    = 0.99
    r0         = 3000.

    N = 10000 # generate 1e4 waveforms

    E0_int = [1.001, 1.1]
    L0_int = [2., 10.]

    print("...generate the parspace")
    pp    = gen_random_pars([Mmin, Mmax], [qmin,qmax], chi_int=[chi_min, chi_max], lambda_int=None, E0_int=E0_int, L0_int=L0_int, precessing=precessing, N=N)
    write_dict_to_txt("ParHyp.txt", pp, 1) # write keys to file

    print("...run")
    with multiprocessing.Pool() as p:
        for j in tqdm.tqdm(
                            [[i, p.apply_async(utils.gen_wf, (pp['m1'][i],pp['m2'][i], pp['chi1z'][i],pp['chi2z'][i], 0., 0.,{'H_hyp':pp['E0'][i], 'j_hyp':pp['L0'][i], 'r_hyp':r0, 'use_mode_lm':k}))]
                                           for i in range(N)]
                            ):
            j[1].get()
            this_par = {}
            for key in pp.keys():
                this_par[key] = pp[key][j[0]]
            write_dict_to_txt("ParHyp.txt", this_par, 0)

def TestParspaceBNS(precessing):
    """
    Generate 1e4 precessing waveforms with parameters within standard BNS bounds
    """
    import multiprocessing; import tqdm
    import utilities as utils

    modes = [[2,2]]
    k = utils.modes_to_k(modes)

    Mmin, Mmax = 1.5, 3.4
    qmin, qmax = 1, 2
    lambda_min = 3
    lambda_max = 5000
    chi_min    = 1e-3
    chi_max    = 0.6

    N          = 10000 # generate 1e4 waveforms

    print("...generate the parspace")
    pp = gen_random_pars([Mmin, Mmax], [qmin,qmax], chi_int=[chi_min, chi_max], lambda_int=[lambda_min, lambda_max], precessing=precessing, N=N)
    write_dict_to_txt("ParBNS.txt", pp, 1) # write keys to file

    # start generating waveforms
    print("...run")
    with multiprocessing.Pool() as p:
        for j in tqdm.tqdm(
                            [[i, p.apply_async(utils.gen_wf, (pp['m1'][i],pp['m2'][i], pp['chi1z'][i],pp['chi2z'][i], pp['LambdaAl2'][i], pp['LambdaBl2'][i],{'ecc':1e-8, 'use_mode_lm':k}, True))]
                                           for i in range(N)]
                            ):
            j[1].get()
            this_par = {}
            for key in pp.keys():
                this_par[key] = pp[key][j[0]]
            write_dict_to_txt("ParBNS.txt", this_par, 0)

def TestParspaceBHNS(precessing):
    """
    Generate 1e4 precessing waveforms with parameters within standard BHNS bounds
    """
    import multiprocessing; import tqdm
    import utilities as utils

    modes = [[2,1], [2,2], [3,2], [3,3], [4,4]]
    k = utils.modes_to_k(modes)

    Mmin, Mmax = 3, 40
    qmin, qmax = 1, 10
    lambda_min = 3
    lambda_max = 5000
    chi_min    = 1e-3
    chi_max    = 0.99

    N = 10000 # generate 1e4 waveforms

    print("...generate the parspace")
    pp = gen_random_pars([Mmin, Mmax], [qmin,qmax], chi_int=[chi_min, chi_max], lambda_int=[lambda_min, lambda_max], precessing=precessing, N=N)
    write_dict_to_txt("ParBHNS.txt", pp, 1) # write keys to file

    # start generating waveforms
    print("...run")
    with multiprocessing.Pool() as p:
        for j in tqdm.tqdm(
                            [[i, p.apply_async(utils.gen_wf, (pp['m1'][i],pp['m2'][i], pp['chi1z'][i],pp['chi2z'][i], 0., pp['LambdaBl2'][i],{'ecc':1e-8, 'use_mode_lm':k}))]
                                           for i in range(N)]
                            ):
            j[1].get()
            this_par = {}
            for key in pp.keys():
                this_par[key] = pp[key][j[0]]
            write_dict_to_txt("ParBHNS.txt", this_par, 0)

def PlotPrecWF():
    """
    Plot one example BBH precessing waveform
    """
    M = 250
    q = 1.1
    chiA = [0.8, 0., 0.]
    chiB = [0., -0.7, -0.1]
    par = utils.CreateDict(M, q, chiA, chiB, lambda1=0, lambda2=0, iota=np.pi/3, f0=11., srate=4096., coa=0., interp="yes", domain=0, modes=[0,1,2,3,4,5,6,7,8])
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

    par_p = utils.CreateDict(M, q, chiA,  chiB,  lambda1=0, lambda2=0, iota=np.pi/3, f0=20., srate=4096., coa=0., interp="yes", domain=0, modes=[0,1,2,3,4,5,6,7,8,13], argout="no")
    par_t = utils.CreateDict(M, q, chiAt, chiBt, lambda1=0, lambda2=0, iota=np.pi/3, f0=20., srate=4096., coa=0., interp="yes", domain=0, modes=[0,1,2,3,4,5,6,7,8,13], argout="no")    
    par_a = utils.CreateDict(M, q, chiAa, chiBa, lambda1=0, lambda2=0, iota=np.pi/3, f0=20., srate=4096., coa=0., interp="yes", domain=0, modes=[0,1,2,3,4,5,6,7,8,13], argout="no")

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
    part     = utils.CreateDict(M, q, chiA, chiB, lambda1=400, lambda2=400, iota=np.pi/2.5, f0=f0, srate=4096., df=0., interp="yes",  domain=0,  modes=[1], coa=np.pi/4, argout="yes")
    parf     = utils.CreateDict(M, q, chiA, chiB, lambda1=400, lambda2=400, iota=np.pi/2.5, f0=f0, srate=4096., df=df, interp="yes" ,  domain=1, modes=[1], coa=np.pi/4, argout="yes")
    
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

    par     = utils.CreateDict(M, q,    chiA,     chiB,     lambda1=0, lambda2=0, iota=np.pi/3, f0=20., srate=4096., df=0, interp="yes", domain=0, modes=[0,1], coa=phi_ref, argout="yes")
    par_rot = utils.CreateDict(M, 1./q, chiB_rot, chiA_rot, lambda1=0, lambda2=0, iota=np.pi/3, f0=20., srate=4096., df=0, interp="yes", domain=0, modes=[0,1], coa=phi_ref+np.pi, argout="yes")

    t1, hp1,_,_,_ = EOB.EOBRunPy(par)
    t2, hp2,_,_,_ = EOB.EOBRunPy(par_rot)

    plt.plot(t1, hp1, color='r', label=r'$m_1 > m_2$', linewidth=2)
    plt.plot(t2, hp2, color='royalblue', linestyle='--', linewidth=2, label=r'$m_2 > m_1$')
    plt.grid()
    plt.xlabel(r'$t$ [s]')
    plt.ylabel(r'$h_+$')
    plt.legend()
    plt.show()

def TestTetradConventions():
    """
    Test tetrad conventions by plotting phase differences between 22 and higher
    modes during the inspiral
    """
    import utilities as ut

    modes = [0, 1, 3, 4, 7, 8, 13]
    par   = utils.CreateDict(100, 1.5, [0.,0., 0.], [0., 0., 0.], lambda1=0, lambda2=0, iota=0., f0=20., srate=4096., df=0, interp="yes", domain=0, modes=modes, coa=0, argout="yes")

    t,_,_,hlm,_ = EOB.EOBRunPy(par)

    # Extract h22 and remove from modes list
    phi22 = -1*hlm['1'][1]
    modes.remove(1)

    # plot
    fig, ax = plt.subplots(len(modes), 1, sharex=True)
    for i, k in enumerate(modes):
        emm   = ut.k_to_emm(k)
        ell   = ut.k_to_ell(k)
        philm = -1*hlm[str(k)][1] 
        diff  = np.mod((2.*philm-float(emm)*phi22),2.*np.pi)
        ax[i].plot(t, diff)
        ax[i].set_ylabel(f"({ell},{emm})")
    
    ax[-1].set_xlabel('t [s]')
    fig.suptitle(r'$2\phi_{lm}-m\phi_{22}$')
    plt.show()

def TestAnomaly():
    """
    Check the variation of the waveform with the anomaly
    """
    import utilities as ut
    modes = [1]
    par   = utils.CreateDict(60, 1.5, [0., 0., 0.], [0., 0., 0.], lambda1=0, lambda2=0, iota=0., 
                             f0=20., srate=4096., df=0, interp="yes", domain=0, modes=modes, coa=0, argout="yes",
                             ecc=0.3)

    n         = 30
    anomalies = np.linspace(0., 2*np.pi, n)
    colors    = plt.cm.magma(np.linspace(0,1,n))

    fig, ax = plt.subplots()
    for z,c in zip(anomalies, colors):
        par['anomaly'] = z
        
        t,_,_,hlm,dyn  = EOB.EOBRunPy(par)
        # plot waveform
        ax.plot(t, hlm['1'][0], c=c)
    
    for z in [0., np.pi]:
        par['anomaly'] = z
        t,_,_,hlm,dyn  = EOB.EOBRunPy(par)
        ax.plot(t, hlm['1'][0], c='k', linewidth=2, 
                linestyle='--',
                label=f'$z = {z}')
        
    ax.set_xlabel('$t [s]$')
    ax.set_ylabel(r'$h_{22}/\nu$')
    plt.show()

def TestEccInitialConditions(plot=False):
    """
    Find the eccentric parameter space corresponding to physical initial conditions.
    Determine whether the desired eccentricity and orbit-averaged frequency imply:
        - a valid initial semilatus rectum p0
        - a valid initial pr*
    This comes from checking that:
        - there exist a value of p0 such that Omg_{EOB}(r(p0,e), j0(p0,e)) - Omg_{target} = 0
        - the value of p0 thus found implies that the periastron is located outside of the EOB horizon
    """
    import json; import tqdm
    
    Msuns       = 4.9254910255435e-6
    Mmin, Mmax  = 20, 400
    qmin, qmax  = 1, 10
    chi_min     = -0.99
    chi_max     = 0.99
    ecc_int     = [0., 0.95]
    anomaly_int = [0., 2*np.pi]
    f0          = 10.
    N = 100000 # look at 1e5 systems
    print("...generate the parspace")
    pp    = gen_random_pars([Mmin, Mmax], [qmin,qmax], chi_int=[chi_min, chi_max], lambda_int=None, ecc_int=ecc_int, anomaly_int=anomaly_int, precessing=0, N=N)

    print("...run")
    errors = []
    
    for n in tqdm.tqdm(range(N)):

        failed      = False
        par         = {k: pp[k][n] for k in pp.keys()}
        nu          = par['q']/(1+par['q'])**2
        omg_orb_avg = f0*np.pi*(par['M']*Msuns)
        par['omega_orb_avg'] = omg_orb_avg
        p0_min = 2.

        # average omega from orbit-average omega
        omg_a = eob_dyn_omg_from_omgbar(omg_orb_avg, par['ecc'], np.pi)
        omg_p = eob_dyn_omg_from_omgbar(omg_orb_avg, par['ecc'], 0.)
        omg_avg_target = 0.5*(omg_a+omg_p)

        # Compute the initial conditions (grid search, for speed)
        p0      = np.linspace(p0_min, 100., 5000)
        omg_avg = []
        for this_p0 in p0:
            # compute j0(p0, ecc)
            j0  = EOB.eob_dyn_j0_py(this_p0, par['q'], par['chi1z'], par['chi2z'], par['ecc'])
            if np.isnan(j0):
                omg_avg.append(0)
                continue
            # compute omega_orb average (between p and a)
            this_omg = 0
            for rr in [this_p0/(1-par['ecc']), this_p0/(1+par['ecc'])]:
                hatH = Espin(rr, j0, 0., par['q'], par['chi1z'], par['chi2z'])
                this_omg += hatH[5]/(nu*hatH[0])
            omg_avg.append(0.5*this_omg)
        omg_avg = np.array(omg_avg)

        # find sample closest to 0, J0 and E0
        idx = np.argmin(np.abs(omg_avg-omg_avg_target))
        J0  = EOB.eob_dyn_j0_py(p0[idx], par['q'], par['chi1z'], par['chi2z'], par['ecc'])
        E0  = nu*Espin(p0[idx]/(1+par['ecc']), J0, 0., par['q'], par['chi1z'], par['chi2z'])[0]
        rs  = np.linspace(2., 30., 100)
        prs = np.linspace(-0.7, 0.7, 3000)

        # check zero: is the energy < 1?
        if E0 > 1.:
            par['notes'] = "E0 > 1"
            par['error_key'] = 0
            errors.append(par)
            failed = True

        # check one: is there a solution?
        # we need some points > 0 and some < 0 for this to happen
        plus  = omg_avg[np.where(omg_avg - omg_avg_target > 0)]
        minus = omg_avg[np.where(omg_avg - omg_avg_target < 0)]
        if len(plus) == 0 or len(minus) == 0:
            par['notes'] = "No solution found for omega_orb_avg"
            par['error_key'] = 1
            if not failed: errors.append(par)
            if not failed: failed = True    

        # check two: can we find pr?
        rma   = p0[idx]/(1+par['ecc']*np.cos(par['anomaly']))
        hatH  = np.array([Espin(rma, J0, pr, par['q'], par['chi1z'], par['chi2z'])[0] for pr in prs])
        plus  = prs[np.where(hatH*nu - E0 > 0)]
        minus = prs[np.where(hatH*nu - E0 < 0)]

        if len(plus) == 0 or len(minus) == 0:
            # is pr=0 a solution?
            if abs(Espin(rma, J0, 0, par['q'], par['chi1z'], par['chi2z'])[0]*nu - E0) > 1e-5:
                # this is a failed case. Did we already fail?
                if not failed:
                    par['notes'] = "No solution found for pr"
                    par['error_key'] = 2
                    errors.append(par)
                    failed = True

        if (failed and plot):

            _, (ax1, ax2) = plt.subplots(2, 1)
            ax1.plot(p0, omg_avg-omg_avg_target)
            ax1.axhline(0., linestyle='--', color='k')
            ax1.set_xlabel(r'$p0$')
            ax1.set_ylabel(r'$\omega_{\rm orb}-\omega_{\rm orb}^{\rm target}$')
            
            hatH = [Espin(x, J0, 0., par['q'], par['chi1z'], par['chi2z'])[0] for x in rs]
            ax2.plot(rs, np.array(hatH)*nu)
            ax2.axhline(E0, linestyle='--', color='r')
            ax2.axvline(p0[idx]/(1+par['ecc']), linestyle='-', color='gray', linewidth=1.5)
            ax2.axvline(p0[idx]/(1-par['ecc']), linestyle='-', color='gray', linewidth=1.5)
            ax2.set_xlabel(r'$r$')
            ax2.set_ylabel(r'$\hat{W}$, effective radial potential')

            _, ax = plt.subplots()
            rma   = p0[idx]/(1+par['ecc']*np.cos(par['anomaly']))
            hatH  = np.array([Espin(rma, J0, pr, par['q'], par['chi1z'], par['chi2z'])[0] for pr in prs])
            ax.plot(prs, hatH*nu)
            ax.axhline(E0, linestyle='--', color='r')
            ax.set_xlabel(r'$p_r$')
            ax.set_ylabel(r'$\hat{H}$')
            plt.show()

        #t,_,_,hlm,dyn = EOB.EOBRunPy(this_p)
    
    print("Errors found:", len(errors))
    # dump errors to file
    with open("ic_errors.json", "w") as f:
        json.dump(errors, f, indent=2)

def PlotInitialConditions():
    """
    Plot the failed/succesful initial conditions;
    Compare to the excluded region in the (e, omega_orb_avg) plane
    """
    import json; import corner
    import matplotlib.pyplot as plt

    with open("ic_errors.json", "r") as f:
        errors = json.load(f)

    eccs  = [par['ecc'] for par in errors]
    zetas = [par['anomaly'] for par in errors]
    omgs  = [par['omega_orb_avg'] for par in errors]
    chis1 = [par['chi1z'] for par in errors]
    chis2 = [par['chi2z'] for par in errors]
    q     = [par['q'] for par in errors]
    omg0  = [par['omega_orb_avg'] for par in errors]

    # check whether the errors below the excluded region are truly
    # due to the initial conditions
    matrix = [q, chis1, chis2, eccs, zetas, omg0]
    labels = [r'$q$', r'$\chi_1$', r'$\chi_2$', 'e', 
              r'$\zeta$', r'$\omega_{\rm orb avg}$'
            ]

    fig = corner.corner(np.transpose(matrix),
                                    labels          = labels,
                                    label_kwargs    = {'size':12.},
                                    hist2d_kwargs     = {'label':labels},
                                    plot_datapoints = True,
                                    show_titles     = False,
                                    plot_contours   = False,
                                    plot_density    = False,
                                    color           = 'red'
                    )
    
    plt.savefig('ic_errors.png')

    def e_cut(omg, chieff):
        A0, A1, A2, A3 = -1.43957927,0.28010553,6.58872891,-9.66104022
        B0, B1, B2, B3 =  0.07065843,0.01772211,0.00653292, 0.87136482
        A = A3 + A2*chieff + A1*chieff**2 + A0*chieff**3
        B = B3 + B2*chieff + B1*chieff**2 + B0*chieff**3
        return A*omg + B

    # now plot the parspace cut
    _, ax = plt.subplots()
    omg_lsp = np.linspace(min(omgs), max(omgs), 100)
    c = ax.scatter(omgs, eccs, c=chis1, label='Failed', alpha=0.4)
    cbar = plt.colorbar(c)
    cbar.set_label(r'$\chi_1$')
    ax.plot(omg_lsp, e_cut(omg_lsp, -0.9), linestyle='--', color='k', label=r'Cut $\chi =-0.9$')
    ax.set_xlabel(r'$M \Omega_{\rm orb avg}$')
    ax.set_ylabel(r'$e$')
    ax.legend()
    plt.savefig('ic_errors_cut.png')
    plt.close()

    pass

def FindInitialConditionsCut():
    """
    Find the parspace cut for the initial conditions:
        - Find the cut in the (e, omega_orb_avg) plane for fixed values of the effective spin;
        - Perform global fits to find the coefficients of the cut
        - Check that the cut is consistent with the excluded region
    """
    import json; from scipy.optimize import curve_fit
    with open("ic_errors.json", "r") as f:
        errors = json.load(f)

    def line(x, A, B):
        return A*x + B
    def cubic(x, A, B, C, D):
        return A*x**3 + B*x**2 + C*x + D
    
    eccs  = [par['ecc'] for par in errors]
    omgs  = [par['omega_orb_avg'] for par in errors]
    chis1 = [par['chi1z'] for par in errors]
    chis2 = [par['chi2z'] for par in errors]
    q     = [par['q'] for par in errors]
    omg0  = [par['omega_orb_avg'] for par in errors]    

    chi_effs_thrs = np.linspace(-0.95, 0.95, 40)

    colors = plt.cm.plasma(np.linspace(0,1,len(chi_effs_thrs)))

    emms = []
    qs   = []

    for kk, this_chi_eff in enumerate(chi_effs_thrs):
        q = np.array(q)
        m1      = q/(1.+q)
        m2      = 1./(1.+q) 
        chi_eff = (m1*chis1 + m2*chis2)

        # choose the bins for the 2d hist
        eccs_bins = np.linspace(0., 0.95, 40)
        omg_bins  = np.linspace(min(omgs), max(omgs), 40)        
        # mask
        mask = np.where(chi_eff > this_chi_eff)
        these_eccs = np.array(eccs)[mask]
        these_omgs = np.array(omgs)[mask]
        # fill the 2d hist
        hist, _, _ = np.histogram2d(these_eccs, these_omgs, bins=[eccs_bins, omg_bins])
        # find the smallest filled value of e for each omega
        e_cut = []
        for i in range(len(omg_bins)-1):
            mask2 = np.where(hist[:,i] > 0)
            if len(mask2[0]) == 0:
               this_ecc = 0.95
            else:
                this_ecc = min(eccs_bins[mask2])
            e_cut.append(this_ecc)
        
        # fit these values with a line
        popt, pcov = curve_fit(line, omg_bins[:-1], e_cut)
        plt.plot(omg_bins[:-1], e_cut, color=colors[kk], alpha=0.6)
        # error bars
        upper = line(omg_bins[:-1], popt[0] + pcov[0,0]**0.5, popt[1] - pcov[1,1]**0.5)
        lower = line(omg_bins[:-1], popt[0] - pcov[0,0]**0.5, popt[1] + pcov[1,1]**0.5)
        plt.fill_between(omg_bins[:-1], upper, lower, color=colors[kk], alpha=0.6)
        plt.plot(omg_bins[:-1], line(omg_bins[:-1], *popt), color=colors[kk], linestyle='--')
        emms.append(popt[0])
        qs.append(popt[1])

    plt.xlabel(r'$\Omega_{\rm orb avg}$')
    plt.ylabel(r'$e_{\rm cut}$')
    plt.savefig('ic_cut.png')
    plt.close()

    # fit and plot the coefficients
    popt1, pcov1 = curve_fit(cubic, chi_effs_thrs, emms)
    popt2, pcov2 = curve_fit(cubic, chi_effs_thrs, qs)

    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    ax.plot(chi_effs_thrs, emms, label=r'$A$')
    ax.plot(chi_effs_thrs, cubic(chi_effs_thrs, *popt1), linestyle='--', color='dodgerblue')
    ax2.plot(chi_effs_thrs, qs, color='r', label=r'$B$')
    ax2.plot(chi_effs_thrs, cubic(chi_effs_thrs, *popt2), linestyle='--', color='darkred')
    ax.set_xlabel(r'$\chi_{\rm eff}$')
    ax.set_ylabel(r'$A$')
    ax2.set_ylabel(r'$B$')
    plt.savefig('ic_coeffs.png')
    plt.close()

    # final cut:
    def e_cut(omg, chi_eff):
        A = cubic(chi_eff, *popt1)
        B = cubic(chi_eff, *popt2)
        e = A*omg + B
        return e
    
    # plot the parspace cut
    _, ax = plt.subplots()
    for kk, chi_eff in [(0, -0.85), (10, -0.5), (20, 0.), (30, 0.5), (39, 0.85)]:
        q       = np.array(q)
        m1      = q/(1.+q)
        m2      = 1./(1.+q) 
        t_chi_eff = (m1*chis1 + m2*chis2)
        mask    = np.where(t_chi_eff > chi_effs_thrs[kk])
        these_eccs = np.array(eccs)[mask]
        these_omgs = np.array(omgs)[mask]
        ax.scatter(these_omgs, these_eccs, color=colors[kk], alpha=0.4)
        ax.plot(omg_bins[:-1], e_cut(omg_bins[:-1], chi_eff), color=colors[kk], linestyle='--', label=f'$\chi_{{\\rm eff}} = {chi_eff}$')
    plt.legend()
    ax.set_xlabel(r'$\Omega_{\rm orb avg}$')
    ax.set_ylabel(r'$e$')
    plt.savefig('ic_cut_final.png')
    plt.close()

    print("Coefficients:")
    print('*\tA = A3 + A2*chi_eff + A1*chi_eff^2 + A0*chi_eff^3')
    print(" \tA0, A1, A2, A3", popt1)
    print('*\tB = B3 + B2*chi_eff + B1*chi_eff^2 + B0*chi_eff^3')
    print(" \tB0, B1, B2, B3", popt2)
    print("Final relation: e_cut = A*omg + B")
    pass

if __name__ == "__main__":

    # run some tests

    if 0:
        print("##### Generate the BBH parameter space #####")
        TestParspaceBBH(0, 1)
        print("...done")

    if 0:
        print("##### Generate the hyperbolic BBH parameter space #####")
        TestParspaceHyp(0)
        print("...done")

    if 0:
        print("##### Generate the BNS parameter space #####")
        TestParspaceBNS(0)
        print("...done")

    if 0:
        print("##### Generate the BHNS parameter space #####")
        TestParspaceBHNS(0)
        print("...done")

    if 0:
        print("##### Plot the BBH parameter space #####")
        PlotParspace("ParBBH.txt", precessing="no")
        print("...done")

    if 1:
        print("##### Plot the BBH eccentric parameter space #####")
        PlotParspace("ParBBH_ecc.txt", precessing="no", eccentric="yes")
        print("...done")

    if 0:
        print("##### Plot the hyp parameter space #####")
        PlotParspace("ParHyp.txt", precessing="no", hyperbolic="yes")
        print("...done")

    if 0:
        print("##### Plot the BNS parameter space #####")
        PlotParspace("ParBNS.txt", tides="yes", precessing="no")
        print("...done")
    
    if 0:
        print("##### Plot the BHNS parameter space #####")
        PlotParspace("ParBHNS.txt", tides="yes", precessing="no")
        print("...done")

    if 0:
        print('##### Plot Prec WF #####')
        PlotPrecWF()
        print("...done")

    if 0:
        print('##### Test the spin aligned limit')
        PlotAlignedSpinLimit()
        print('...done')

    if 0:
        print('##### Plot a FD SPA WF, compare to TD #####')
        PlotSPAWF()
        print("...done")

    if 0:
        print("##### Test Swap #####")
        SwapPrecWF()
        print("...done")
    
    if 0:
        print("##### Test Tetrad conventions #####")
        TestTetradConventions()
        print("...done")

    if 0:
        print('##### Test True Anomaly #####')
        TestAnomaly()
        print("...done")

    if 0:
        print('##### Test Eccentric ICs Failures #####')
        TestEccInitialConditions(plot=False)
        print("...done")

    if 0:
        print('##### Plot eccentric IC failures #####')
        PlotInitialConditions()
        print("...done")

    if 0:
        print('##### Find ICs Cut #####')
        FindInitialConditionsCut()
        print("...done")