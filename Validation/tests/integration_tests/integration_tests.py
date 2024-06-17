"""
List all tests to be run here
Run with
    pytest -v integration_tests.py

RG, 12/2023
"""

import pytest
import numpy as np
import utilities as utils
params = {

        # test the same waveform is always the same (if not passed can point to leaks)
        'test_same_wf':     [ (40., 20., 0.8, 0.5, 0., 0.),   #BBH
                             (1.7, 1., 0.2, 0.1, 400, 1500.), #BNS
                             (6., 1.4, 0.8, 0.1, 0., 1500.)], #BHNS
        # test the PA approximation (add as required)
        'test_PA_phasing': [ ((2., 1.4, 0.1, 0.1,   1000, 1000), 5e-3),
                             ((2., 1.4, 0.99, 0.99, 10,   5000), 5e-3)],
        # test specific binaries (add as required)
        # format: ((m1, m2, s1z, s2z, l1, l2), expected # of output arrays)
        'test_gen_wf':     [ ((397.75416888160163, 4.45654408641580,  -0.96955312160861, 0.82775370084938, 0., 0.), 3), 
                             ((1.4, 1.4, 0., 0., 1000, 1000), 3),
                             ((2., 1.4, 0.1, 0.1, 1000, 1000), 3)],
        
        # test swap of components for BBH, BNS and BHNS
        # format: (m1, m2, s1z, s2z, l1, l2)
        'test_swap':       [ (40., 20., 0.8, 0.5, 0., 0.),    #BBH
                            (1.7, 1., 0.2, 0.1, 400, 1500.)], #BNS
                            #(6., 1.4, 0.8, 0.1, 0., 1500.)], #BHNS always assumes m1 is the BH
        
        # test distance scaling
        # format: (m1, m2, s1z, s2z, l1, l2), reference distance, list of distances to test
        'test_distance':   [ ((40., 20., 0.8, 0.5, 0., 0.),    100., (200., 1e3, 1e5, 1e9)),  #BBH
                             ((1.7, 1., 0.2, 0.1, 400, 1500.), 100., (200., 1e3, 1e5, 1e9)),  #BNS
                             ((6., 1.4, 0.8, 0.1, 0., 1500.),  100., (200., 1e3, 1e5, 1e9))], #BHNS

        # test behavior under small perturbations of the parameters
        # format: (m1, m2, s1z, s2z, l1, l2), maximum tolerance, Number of draws, max variation of the parameters
        'test_parameters_perturb'    : [ ((40., 20., 0.8, 0.5, 0., 0.),        1e-3, 3, [1e-9, 1e-10, 1e-11]),  #BBH
                                         ((1.7, 1.,  0.2, 0.1, 400, 1500.),    1e-3, 3, [1e-9, 1e-10, 1e-11]),  #BNS
                                         ((6., 1.4,  0.8, 0.1, 0., 1500.),     1e-3, 3, [1e-9, 1e-10, 1e-11])], #BHNS

        # Test that the plus and cross polarizations are reconstructed correctly
        'test_polariz_reconstruction': [ ((40., 20., 0.8, 0.5, 0., 0.),     1e-10),  #BBH
                                         ((1.7, 1.,  0.2, 0.1, 400, 1500.), 1e-10),  #BNS
                                         ((6., 1.4,  0.8, 0.1, 0., 1500.),  1e-10)], #BHNS

        # Test the periodicity of 22-only waveforms with phase shifts multiple of Pi
        'test_phase_shifts_22'       :  [((40., 20., 0.8, 0.5, 0., 0.),     1e-10, [np.pi, 2*np.pi, 3*np.pi]),  #BBH
                                         ((1.7, 1.,  0.2, 0.1, 400, 1500.), 1e-10, [np.pi, 2*np.pi, 3*np.pi]),  #BNS
                                         ((6., 1.4,  0.8, 0.1, 0., 1500.),  1e-10, [np.pi, 2*np.pi, 3*np.pi])], #BHNS

        # Test the periodicity of (ell, emm)-only waveforms with phase shifts such that phi_r = 2 Pi/emm
        'test_phase_shifts_lm'       :  [((40., 20., 0.8, 0.5, 0., 0.),    1e-10, [0, 1, 4, 8, 13]),   #BBH  (21, 22, 33, 44, 55)
                                         ((1.7, 1.,  0.2, 0.1, 400, 1500.), 1e-10, [0, 1, 4, 8, 13]),  #BNS  (21, 22, 33, 44, 55)
                                         ((6., 1.4,  0.8, 0.1, 0., 1500.),  1e-10, [0, 1, 4, 8, 13])], #BHNS (21, 22, 33, 44, 55)       
    }

@pytest.mark.parametrize("pars, out", params['test_PA_phasing'])
def test_PA_phasing(pars, out):
    """
    Test that two waveforms computed with/without PA are
    consistent with one another
    """
    mass1, mass2, s1z, s2z, lam1, lam2 = pars
    _, hp_pa0, hc_pa0 = utils.gen_wf(mass1, mass2, s1z, s2z, lam1, lam2, additional_pars={'use_postadiabatic_dynamics':"no"}, return_zero=False)
    _, hp_pa1, hc_pa1 = utils.gen_wf(mass1, mass2, s1z, s2z, lam1, lam2, additional_pars={'use_postadiabatic_dynamics':"yes"}, return_zero=False)
    phase_pa0 = np.unwrap(np.angle(hp_pa0 - 1j*hc_pa0))
    phase_pa1 = np.unwrap(np.angle(hp_pa1 - 1j*hc_pa1))

    if len(phase_pa0) != len(phase_pa1):
        assert False
    deltap = abs(phase_pa0 - phase_pa1)
    assert (deltap < out).all()

@pytest.mark.parametrize("pars, out", params['test_gen_wf'])
def test_gen_wf(pars, out):
    """
    Test the waveform generation for the binary specified by pars
    """
    m1, m2, s1z, s2z, l1, l2 = pars
    if m1+m2 > 100:
        result = utils.gen_wf(m1, m2, s1z, s2z, l1, l2, additional_pars={'initial_frequency':5},return_zero=False)
    else:
        result = utils.gen_wf(m1, m2, s1z, s2z, l1, l2, return_zero=False)
    assert len(result) == out

@pytest.mark.parametrize("pars", params['test_swap'])
def test_swap(pars):
    """
    Test that waveforms for swapped binaries are consistent.
    According to LAL conventions, if m1<->m2 the x axis is flipped.
    Therefore, the waveform has to remain identical when:
        - label_1<->label_2;
        - coalescence_angle -> coalescence_angle + Pi;
        - In plane spins are rotated by Pi
    """
    m1, m2, s1z, s2z, l1, l2 = pars
    _,_,h1 = utils.gen_wf(m1, m2, s1z, s2z, l1, l2, additional_pars={'use_geometric_units':"yes", 
                                                                     'initial_frequency':0.003}, 
                                                    return_zero=False)
    _,_,h2 = utils.gen_wf(m2, m1, s2z, s1z, l2, l1, additional_pars={'coa': np.pi, 
                                                                     'use_geometric_units':"yes", 
                                                                     'initial_frequency':0.003}, 
                                                    return_zero=False)

    assert np.allclose(h1, h2, atol=1e-15, rtol=1e-15)

@pytest.mark.parametrize("pars, ref, tst", params['test_distance'])
def test_distance_scaling(pars, ref, tst):
    """
    Test that the amplitude scales correctly with the distance.
    """
    m1, m2, s1z, s2z, l1, l2 = pars
    reference_distance = ref
    test_distance      = tst

    _,hp_ref,hc_ref = utils.gen_wf(m1, m2, s1z, s2z, l1, l2, additional_pars={'distance':reference_distance}, return_zero=False) 
    A_ref = np.sqrt(hp_ref**2+hc_ref**2)
    for dis in test_distance:
        _,hp_tst,hc_tst = utils.gen_wf(m1, m2, s1z, s2z, l1, l2, additional_pars={'distance':dis}, return_zero=False)
        A_tst = np.sqrt(hp_tst**2+hc_tst**2)
        scale_factor = dis/reference_distance
        if not np.allclose(A_ref/A_ref[1], scale_factor*A_tst/A_ref[1], atol=1e-15, rtol=1e-15): 
            print("Failing for a test distance of", dis)
            assert False
    assert True

def test_srate_change():
    """
    Test that waveforms are unmodified when changing the sampling rate
    """
    t1,_, h1 = utils.gen_wf(1.4, 1.4, 0, 0, 400, 400,  additional_pars={'srate_interp':8192.}, return_zero=False)
    t2,_, h2 = utils.gen_wf(1.4, 1.4, 0, 0, 400, 400,  additional_pars={'srate_interp':4096.}, return_zero=False)
    h2_int   = np.interp(t2, t1, h1)
    assert t1[1]-t1[0] == 1./8192.
    assert t2[1]-t2[0] == 1./4096.
    assert np.allclose(h2/h2_int, np.ones(len(h2_int)), atol=1e-15, rtol=1e-15)

@pytest.mark.parametrize("pars, ref, N, tst", params['test_parameters_perturb'])
def test_parameters_perturb(pars, ref, N, tst):
    """
    Test that waveforms are approximately unchanged with small variations of the parameters.
    Note: we do NOT modify lambda1, as the BHNS to BNS and BBH to BNS limits are known to be 
          discontinuous (different rad. reaction etc)
    """
    m1, m2, s1z, s2z, l1, l2 = pars
    threshold                = ref
    epsilon                  = tst

    t_ref,hp_ref,hc_ref = utils.gen_wf(m1, m2, s1z, s2z, l1, l2, return_zero=False,
                          additional_pars={'use_geometric_units':"yes",
                                           'initial_frequency':0.001
                                           }
                            )
    A_ref = np.sqrt(hp_ref**2+hc_ref**2)
    for eps in epsilon:
        for _ in range(N):
            m1e, m2e = m1  + np.random.uniform(0,1)*eps,  m2 + np.random.uniform(0,1)*eps
            s1e, s2e = s1z + np.random.uniform(0,1)*eps, s2z + np.random.uniform(0,1)*eps
            l2e      = l2  + np.random.uniform(0,1)*eps
            t_tst,hp_tst,hc_tst = utils.gen_wf(m1e, m2e, s1e, s2e, l1, l2e, return_zero=False,
                                additional_pars={'use_geometric_units':"yes",
                                'initial_frequency':0.001},
                                )
            A_tst = np.sqrt(hp_tst**2+hc_tst**2)
            if not np.allclose(max(A_ref), max(A_tst), atol=threshold, rtol=threshold): 
                print("Amplitude failing for a perturbation of", eps)
                print(max(A_ref), max(A_tst))
                assert False
            if not np.allclose(t_ref[np.argmax(A_ref)], t_tst[np.argmax(A_tst)], atol=threshold, rtol=threshold): 
                print("Merger time failing for a perturbation of", eps)
                print(t_ref[np.argmax(A_ref)], t_tst[np.argmax(A_tst)])
                assert False
    assert True

@pytest.mark.parametrize("pars, eps", params['test_polariz_reconstruction'])
def test_polariz_reconstruction(pars, eps):
    """
    Test that the plus and cross polarizations are reconstructed as expected from the
    modes dictionary hlm
    """
    m1, m2, s1z, s2z, l1, l2 = pars
    iota = np.random.uniform(0, 2*np.pi)
    phi  = np.random.uniform(0, 2*np.pi)
    modes = [0,1,5,8,13]
    t,hp,hc,hlm,_ = utils.gen_wf(m1, m2, s1z, s2z, l1, l2, return_zero=False,
                                additional_pars={'arg_out':"yes",
                                                 'use_mode_lm': modes,
                                                 'use_geometric_units':'yes',
                                                 'initial_frequency':0.001,
                                                 'coalescence_angle':phi,
                                                 'inclination': iota,
                                                 }
                                )
    hp_rec, hc_rec = utils.compute_hphc(hlm, phi=phi, i=iota, modes=modes)

    # because of the `LVK` tetrad convention, we need to multiply by -1 
    # the reconstructed polarizations
    if not np.allclose(hp, -1*hp_rec, atol=eps, rtol=eps):
        assert False
    if not np.allclose(hc, -1*hc_rec, atol=eps, rtol=eps): 
        assert False
    assert True

@pytest.mark.parametrize("pars, eps, phis", params['test_phase_shifts_22'])
def test_phase_shifts_22(pars, eps, phis):
    """
    Check how phase shifts affect the waveform and ensure that (for 22-only waveforms)
    multiples of np.pi leave the waveform unchanged.
    """
    m1, m2, s1z, s2z, l1, l2 = pars
    modes = [1]
    t,hp,hc= utils.gen_wf(m1, m2, s1z, s2z, l1, l2, return_zero=False,
                                additional_pars={
                                                 'use_mode_lm': modes,
                                                 'use_geometric_units':'yes',
                                                 'initial_frequency':0.003
                                                 }
                                )
    for phi in phis:
        _,hp_p,hc_p = utils.gen_wf(m1, m2, s1z, s2z, l1, l2,return_zero=False,
                                additional_pars={
                                                 'use_mode_lm': modes,
                                                 'coalescence_angle':phi,
                                                 'use_geometric_units':'yes',
                                                 'initial_frequency':0.003
                                                 }
                                )
    if not np.allclose(hp, hp_p, atol=eps, rtol=eps):
        assert False
    if not np.allclose(hc, hc_p, atol=eps, rtol=eps):
        assert False
    assert True

@pytest.mark.parametrize("pars, eps, modes", params['test_phase_shifts_lm'])
def test_phase_shifts_lm(pars, eps, modes):
    """
    Check how phase shifts affect the waveform and ensure that (for ell emm-only waveforms)
    multiples of 2*np.pi/emm leave the waveform unchanged.
    """
    m1, m2, s1z, s2z, l1, l2 = pars    
    t,hp,hc,hlm,_= utils.gen_wf(m1, m2, s1z, s2z, l1, l2, return_zero=False,
                                additional_pars={
                                                 'use_mode_lm': modes,
                                                 'arg_out':'yes',
                                                 'use_geometric_units':'yes',
                                                 'initial_frequency':0.003
                                                 }
                                )
    for k in modes:
        emm    = utils.k_to_emm(k)
        phi_lm = 2*np.pi/emm
        print(emm, phi_lm)
        # single reconstruct the single mode wf
        hp_p, hc_p = utils.compute_hphc(hlm, modes=[k], phi=phi_lm)
        hp_0, hc_0 = utils.compute_hphc(hlm, modes=[k], phi=0)
        if not np.allclose(hp_0, hp_p, atol=eps, rtol=eps):
            assert False
        if not np.allclose(hc_0, hc_p, atol=eps, rtol=eps):
            assert False
    assert True

@pytest.mark.parametrize("pars", params['test_same_wf'])
def test_same_wf(pars):
    """
    Test that a waveform generated 10 times is always the same
    (this can flag e.g. leaks or non-allocated memory being 
    accessed)
    """
    mass1, mass2, s1z, s2z, lam1, lam2 = pars
    _, hp0, hc0 = utils.gen_wf(mass1, mass2, s1z, s2z, lam1, lam2, return_zero=False,
                                    additional_pars={
                                                 'use_geometric_units':'yes',
                                                 'initial_frequency':0.003
                                                 }
                                    )
    for _ in range(10):
        _, hp1, hc1 = utils.gen_wf(mass1, mass2, s1z, s2z, lam1, lam2, return_zero=False, 
                                   additional_pars={
                                                 'use_geometric_units':'yes',
                                                 'initial_frequency':0.003
                                                 }
                                    )
        
        assert np.allclose(hp0, hp1, atol=1e-15, rtol=1e-15)
        assert np.allclose(hc0, hc1, atol=1e-15, rtol=1e-15)

def test_phiref_meaning():
    """
    Test that the meaning of the reference phase phi_ref is 
    correct
    """

    phi_r = [np.pi/4, np.pi/27]
    beta  = [] 

    for phi in phi_r:
        t, hp, hc = utils.gen_wf(40., 20., 0.5, -0.1, 0, 0, return_zero=False, additional_pars={'coalescence_angle': phi,
                                                                       'inclination' : 0.,
                                                                       'use_mode_lm' : [1]})
        this_beta = np.unwrap(np.angle(hp - 1j*hc))
        beta.append(this_beta)

    q0 = np.mod(-2.*(phi_r[0]-phi_r[1]),2.*np.pi)
    q1 = np.mod(beta[0][0]-beta[1][0],2.*np.pi)
    assert q0 == q1
    

    


