# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import lal
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, tukey

def downsample(strain, old_sampling_rate, new_sampling_rate):
    factor = int(old_sampling_rate/new_sampling_rate)
    return strain[::factor]

def resize_time_series(inarr, N):
    # zero-pad to the required length
    outarr = np.pad(inarr, (0,N-len(inarr)), mode='constant', constant_values=0)
    return outarr

fname='H-H1_LOSC_4_V1-1126259446-32.txt'
tevent = 1126259462.423
srate=4096

def load_data(fname,
              chunk_size=4.0,
              trigtime=tevent,
              sampling_rate=4096,
              psd_file=None):

    # Extract some metadata from the file name
    ifo,fr_type,starttime,T=fname.strip('.txt').split('-')
    starttime=float(starttime)
    T=float(T)
    print('Loading {0} starting at {1} length {2}s'.format(fname,starttime,T))
    rawstrain = np.loadtxt(fname)
    N=len(rawstrain)
    # sampling timestep (s)
    dt = T/N
    # Sampling rate (Hz)
    srate=1./dt
    
    if sampling_rate is not None:
        strain = downsample(rawstrain, srate, sampling_rate)
        srate = sampling_rate
        dt = 1./srate
    else:
        strain = rawstrain

    time = np.linspace(starttime,starttime+T,len(strain))
    # find the index corresponding to the trigger time
    index_trigtime = int((trigtime-starttime)*srate)
    # Number of samples in the chunk
    chunksize=int(chunk_size*srate)
    # Starting time for the signal chunk
    # We want the trigger time 1s before the end of the segment
    index_chunk_start = index_trigtime - int(srate*(chunk_size-1))
    chunk_start = starttime+dt*index_chunk_start
    # signal chunk
    signal_chunk=np.zeros(chunksize,dtype=np.float64)
    for i in range(chunksize):
        signal_chunk[i] = strain[index_chunk_start+i]

    mask = np.ones(len(strain), dtype=bool)
    mask[range(index_chunk_start,index_chunk_start+chunksize,chunksize)] = False
    # window the data
    padding = 0.1
    window=tukey(chunksize,padding)
    signal_chunk*=window
    # zero-pad to the required length

    windowNorm = chunksize/np.sum(window**2)
    # Compute the frequency domain strain
    sf = np.fft.rfft(signal_chunk)*windowNorm*dt
    frequencies = np.fft.rfftfreq(chunksize, dt)
#    sf = sf[np.where(np.logical_and(freqs>=flow,freqs<fhigh))]

    if psd_file is None:
        # Compute the PSD
        psd, freqs = mlab.psd(strain[mask], Fs = srate, NFFT = np.int(srate), window=tukey(np.int(srate),padding))
        psd_int = interp1d(freqs, psd, bounds_error=False, fill_value=np.inf)
    else:
        f, psd = np.loadtxt(psd_file,unpack=True)
        psd *= psd
        # generate an interpolant for the PSD
        psd_int = interp1d(f, psd, bounds_error=False, fill_value=np.inf)

    # compute times and frequencies for convenience
    
    times = chunk_start+np.linspace(0,chunk_size,chunksize)

    return times, signal_chunk, frequencies, sf, psd_int(frequencies)

def generate_data(psd_file,
                  T=4.0,
                  starttime=1126259446.,
                  sampling_rate=4096.,
                  zero_noise = False):

    f, psd = np.loadtxt(psd_file,unpack=True)
    psd *= psd
    # generate an interpolant for the PSD
    psd_int = interp1d(f, psd, bounds_error=False, fill_value=np.inf)
    df = 1./T
    N = int(sampling_rate*T)
    times = np.linspace(starttime,starttime+T,N)
    # generate the FD noise
    frequencies = np.arange(0,sampling_rate/2.,df)
    if zero_noise is False:
        frequency_series = np.array([np.random.normal(0.0,np.sqrt(psd_int(f)/df/2.))+1j*np.random.normal(0.0,np.sqrt(psd_int(f)/df/2.)) for f in frequencies])
    else: frequency_series = np.zeros(len(frequencies),dtype=np.complex64)
    return times, frequencies, frequency_series, psd_int(frequencies)

# function to writen data
def whiten(strain, interp_psd, dt):
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, dt)
    bb, ab = butter(4, [20/(0.5*srate), 300 / (0.5*srate) ], btype='band')
    strain = filtfilt(bb, ab, strain)
    # whitening: transform to freq domain, divide by asd, then transform back,
    # taking care to get normalization right.
    hf = np.fft.rfft(strain)
    white_hf = hf / (np.sqrt(interp_psd(freqs) /dt/2.))
    white_ht = np.fft.irfft(white_hf, n=Nt)
    return white_ht

if __name__ == "__main__":
    T, F, strainF, psd = generate_data('/Users/wdp/src/lalsuite/lalsimulation/src/LIGO-T0900288-v3-ZERO_DET_high_P.txt', sampling_rate = 4096, starttime=1126259446., T=4, zero_noise = False)
    from matplotlib import pyplot as plt
    plt.figure()
    plt.loglog(F, psd)

    plt.figure()
    plt.plot(F,strainF)
    plt.show()
