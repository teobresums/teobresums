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
from scipy.signal import butter, filtfilt, welch

def chunks(times,strain,chunksize,avoid=None):
    # Skip the 0th chunk which has filter ringing
    if avoid is None:
        avoid = times[0]-1e6 # dummy value
    return (strain[chunksize*j:chunksize*(j+1)]
                   for j in range(1,len(strain)//chunksize)
                       if not times[chunksize*j] < avoid < times[chunksize*(j+1)] )

def resize_time_series(inarr, N):
    # zero-pad to the required length
    outarr = np.pad(inarr, (0,N-len(inarr)), mode='constant', constant_values=0)
    # roll the array so that the peak of the time series is in the center of the frame
    outarr = np.roll(outarr,N//2-np.argmax(outarr))
    return outarr

def fd_from_td(t_arr, strain_t, srate = 4096, N = 4096):
    #    N = int(2**np.ceil(np.log2(len(t_arr)))) # round it up to the next power of two
    strain_t = resize_time_series(strain_t, N)
    strain_t[len(t_arr):] = 0.0
    strain_f = np.fft.rfft(strain_t)/srate #rfft means real fft, i.e. it's optimized when the input is real. It's divided by the srate to respect LAL conventions
    ff = np.fft.rfftfreq(len(strain_t), np.mean(np.diff(t_arr)))
    return (ff, strain_f)

fname='H-H1_LOSC_4_V1-1126259446-32.txt'
tevent = 1126259462.423

srate=4096

def get_bandpassed_strain(fname,low=20,high=2028):
    ifo,fr_type,starttime,T=fname.strip('.txt').split('-')
    starttime=float(starttime)
    T=float(T)
    print('Loading {0} starting at {1} length {2}s'.format(fname,starttime,T))
    rawstrain = np.loadtxt(fname)
    dt = T/len(rawstrain)
    # Sampling rate (Hz)
    srate=1/dt
    bb, ab = butter(4, [low/(0.5*srate), high / (0.5*srate) ], btype='band')
    strain = filtfilt(bb, ab, rawstrain)
    times = np.linspace(starttime,starttime+T,len(strain))
    return times,strain,srate


def load_data(fname, chunk_size=4.0, trigtime=tevent, injection=False):
    """
    Read the strain data and divide it into 
    Nsegments.
    Find the on-source chunk X of length chunk_size.
    Estimate the covariance matrix Sigma
    Returns (X,Sigma)
    If trigtime is not given X is None    
    """
    # Extract some metadata from the file name
    ifo,fr_type,starttime,T=fname.strip('.txt').split('-')
    starttime=float(starttime)
    T=float(T)
    print('Loading {0} starting at {1} length {2}s'.format(fname,starttime,T))
    rawstrain = np.loadtxt(fname)
    # sampling timestep (s)
    dt = T/len(rawstrain)
    # Sampling rate (Hz)
    srate=1/dt
    
    # bandpass between 20 hz and 10 hz below nyquist
    bb, ab = butter(4, [20/(0.5*srate), (0.5*srate-10) / (0.5*srate) ], btype='band')
    strain = filtfilt(bb, ab, rawstrain)
    # Round to integer number of chunks
    Nchunks = np.int(T//chunk_size)

    seglen = np.float(T)/Nchunks
    # Number of samples in each chunk
    chunksize=int(seglen*srate)

    # signal chunk
    signal_chunksize=np.int(srate*chunk_size)

    # Find the on-source chunk, centred on trigtime
    if trigtime is not None:
          i = int((trigtime-starttime)*srate)
          on_source = strain[i-signal_chunksize//2:i+signal_chunksize//2+1]
    else:
        on_source = None
        print('No trigtime given')

    # Compute the times for convenience
    i = int((trigtime-starttime)*srate)
    times = np.linspace(starttime,starttime+T-dt,len(strain))[i-signal_chunksize//2:i+signal_chunksize//2+1]

    # Compute the PSD
    psd, freqs = mlab.psd(strain, Fs = srate, NFFT = np.int(srate))

    # Compute the frequency domain strain
    frequency, on_source_fd = fd_from_td(times, on_source, srate = srate, N = signal_chunksize+1)
    psd_int = np.interp(frequency, freqs, psd)
    return times, on_source, frequency, on_source_fd, psd_int

# function to writen data
def whiten(strain, interp_psd, dt):
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, dt)
    # whitening: transform to freq domain, divide by asd, then transform back,
    # taking care to get normalization right.
    hf = np.fft.rfft(strain)
    white_hf = hf / (np.sqrt(interp_psd(freqs) /dt/2.))
    white_ht = np.fft.irfft(white_hf, n=Nt)
    return white_ht

if __name__ == "__main__":
    times, on_source, frequency, on_source_fd, psd = load_data('data/H-H1_LOSC_4_V1-1126259446-32.txt')
    from matplotlib import pyplot as plt
    plt.loglog(frequency, psd)
    plt.show()
