from pyTEOBResumS import pyTEOBResumS
"""
    h = pyTEOBResumS(m1,
    m2,
    spin1x,
    spin1y,
    spin1z,
    spin2x,
    spin2y,
    spin2z,
    inclination,
    polarisation,
    f_min,
    sampling_rate,
    LambdaAl2,
    LambdaBl2,
    distance,
    NQC,
    tidal,
    speedy,
    RWZ,
    dynamics,
    lm,
    solver_scheme)
"""

def resize_time_series(inarr, N):
    # zero-pad to the required length
    outarr = np.pad(inarr, (0,N-len(inarr)), mode='constant', constant_values=0)
    return outarr


def fd_from_td(t_arr, strain_t, srate = 4096, N = 4096):
#    N = int(2**np.ceil(np.log2(len(t_arr)))) # round it up to the next power of two
    strain_t = resize_time_series(strain_t, N)
    strain_t[len(t_arr):] = 0.0
    strain_f = np.fft.rfft(strain_t)/srate #rfft means real fft, i.e. it's optimized when the input is real. It's divided by the srate to respect LAL conventions
    ff = np.fft.rfftfreq(len(strain_t), np.mean(np.diff(t_arr)))
    return (ff, strain_f)

def FDWaveform(params):
    m1 = params['m1']
    m2 = params['m2']
    spin1z = params['spin1z']
    spin2z = params['spin2z']
    inclination = params['iota']
    polarisation = params['psi']
    f_min = 20.0
    sampling_rate = 4096.
    segment_length = int(sampling_rate*32)
    t = np.linspace(0.0,32,segment_length)
    
    h = pyTEOBResumS(m1,
                     m2,
                     0.0,
                     0.0,
                     spin1z,
                     0.0,
                     0.0,
                     spin2z,
                     inclination,
                     polarisation,
                     f_min,
                     sampling_rate,
                     0.0,
                     0.0,
                     1.0,
                     1,
                     0,
                     1,
                     0,
                     0,
                     0,
                     0)
                     
#    return h
    f, hptilde = fd_from_td(t, h[:,0], srate = 4096, N = segment_length)
    return f, hptilde

if __name__ == "__main__":
    params = {'m1':40.0,
              'm2':4.0,
              'spin1z': 0.9,
              'spin2z': 0.0,
              'iota':0.0,
              'psi':0.0
              }
    from pylab import *
    f, h = FDWaveform(params)

    plot(f, h)

    show()
