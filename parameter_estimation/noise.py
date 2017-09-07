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

def chunks(times,strain,chunksize,avoid=None):
    # Skip the 0th chunk which has filter ringing
    if avoid is None:
        avoid = times[0]-1e6 # dummy value
    return (strain[chunksize*j:chunksize*(j+1)]
                   for j in range(1,len(strain)//chunksize)
                       if not times[chunksize*j] < avoid < times[chunksize*(j+1)] )

def downsample(strain, old_sampling_rate, new_sampling_rate):
    factor = int(old_sampling_rate/new_sampling_rate)
    return strain.reshape(-1, factor).mean(axis=1)

def resize_time_series(inarr, N):
    # zero-pad to the required length
    outarr = np.pad(inarr, (0,N-len(inarr)), mode='constant', constant_values=0)
    # roll the array so that the peak of the time series is in the center of the frame
#    outarr = np.roll(outarr,N//2-np.argmax(np.abs(outarr)))
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
    srate=1/dt
#    bb, ab = butter(4, [20/(0.5*srate), 2028 / (0.5*srate) ], btype='band')
#    strain = filtfilt(bb, ab, rawstrain)

    strain = rawstrain
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
    for i in range(chunksize): signal_chunk[i] = strain[index_chunk_start+i]
    
    # window the data
    padding = 0.5
    window=tukey(chunksize,padding)
    signal_chunk*=window
    # zero-pad to the required length

    windowNorm = chunksize/np.sum(window**2)
    # Compute the frequency domain strain
    sf = np.fft.rfft(signal_chunk)*windowNorm
    # Compute the PSD
    psd, freqs = mlab.psd(strain, Fs = srate, NFFT = np.int(srate))
    psd_int = interp1d(freqs, psd)
    
    # compute times and frequencies for convenience
    
    times = chunk_start+np.linspace(0,chunk_size,chunksize)
    frequencies = np.linspace(0,srate/2.,chunksize/2 +1)

    return times, signal_chunk, frequencies, sf, psd_int(frequencies)

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
    T, strainT, F, strainF, psd = load_data('data/H-H1_LOSC_4_V1-1126259446-32.txt')
    from matplotlib import pyplot as plt
    fmin = 20
    fmax = 2048
    imin = int(fmin*4)
    imax = int(fmax*4)
    print imin, imax
    plt.figure()
    plt.plot(1126259459.42+np.linspace(0,4,int(4*4096)),strainT)
    plt.axvline(tevent)
    plt.figure()
    plt.plot(1126259459.42+np.linspace(0,4,int(4*4096)),whiten(strainT,psd,1./4096.))
    plt.axvline(tevent)
    plt.show()
