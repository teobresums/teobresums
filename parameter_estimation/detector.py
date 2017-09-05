from __future__ import division
import numpy as np
import noise
from pylal import inject
from lal import ComputeDetAMResponse, GreenwichMeanSiderealTime, TimeDelayFromEarthCenter
from scipy.signal import tukey

detMap = {'H1': 'LHO_4k', 'H2': 'LHO_2k', 'L1': 'LLO_4k', 'G1': 'GEO_600', 'V1': 'VIRGO', 'T1': 'TAMA_300'}

class GravitationalWaveDetector(object):
    
    def __init__(self, name, datafile, chunk_size=4.0, trigtime=1126259462.423, **kwargs):
        self.name = name
        self.trigtime = trigtime
        self.lal_detector = inject.cached_detector[detMap[self.name]]
        self.location = self.lal_detector.location
        self.Times, self.TimeSeries, self.Frequency, self.FrequencySeries, self.PowerSpectralDensity = noise.load_data(datafile, chunk_size=chunk_size, trigtime=trigtime, **kwargs)
        self.Epoch = self.Times[0]
        self.dt = np.diff(self.Times)[0]
        self.T = self.Times[-1]-self.Times[0]
        self.df = np.diff(self.Frequency)[0]
        self.sampling_rate = 1./self.dt
        self.segment_length = int(self.T*self.sampling_rate)
        self.SigmaSq = self.PowerSpectralDensity*self.dt**2
        self.Flow = 20.0
        self.Fhigh = self.sampling_rate/2.
        self.kmin = int(self.Flow/self.df)
        self.kmax = int(self.Fhigh/self.df)
    
    def Project(self, hptilde, hctilde, ra, dec, psi, tc):
        """
        projects and timeshifts the GW signal onto the detector
        """
        gmst = GreenwichMeanSiderealTime(tc)
        fp,fc = ComputeDetAMResponse(self.lal_detector.response, ra, dec, psi, gmst)
        
        timeShift = self.Epoch-tc + TimeDelayFromEarthCenter(self.location, ra, dec, tc)

        return np.exp(-1j*2.0*np.pi*timeShift*self.Frequency[self.kmin:self.kmax])*(fp*hptilde[self.kmin:self.kmax]+fc*hctilde[self.kmin:self.kmax])

    def logLikelihood(self, hptilde, hctilde, ra, dec, psi, tc):
        
        TwoDeltaToverN = 2.0 * self.dt / self.segment_length;
        residuals = self.FrequencySeries[self.kmin:self.kmax] - self.Project(hptilde, hctilde, ra, dec, psi, tc)
        
        numerator = residuals*np.conj(residuals)
#        for i,l in enumerate(self.Frequency):
#            print hptilde[i],hctilde[i],numerator[i],self.InversePowerSpectralDensity[i]
#        exit()
        logLseries = -TwoDeltaToverN*np.real(numerator/self.SigmaSq[self.kmin:self.kmax])
        
        for i,l in enumerate(logLseries):
            print self.Frequency[i],logLseries[i],hptilde[i],numerator[i],self.SigmaSq[i]
        print np.sum(logLseries)
        exit()
        return np.sum(logLseries)

if __name__ == "__main__":
    H = GravitationalWaveDetector('H1','data/H-H1_LOSC_4_V1-1126259446-32.txt', trigtime = 1126259462.43)
    L = GravitationalWaveDetector('L1','data/L-L1_LOSC_4_V1-1126259446-32.txt', trigtime = 1126259462.43)
    
    for attr, value in H.__dict__.iteritems():
        print attr, value

    from pyTEOBResumS import pyTEOBResumS
    
    f_min = 20.0
    sampling_rate = 4096.
    
    tc = 1126259462.43
    ra = 1.82161142311
    dec = -1.27626704919
    inc = 2.92033171962
    psi = 0.0
    
    h = pyTEOBResumS(40.0,
                     40.0,
                     0.0,
                     0.0,
                     0.01,
                     0.0,
                     0.0,
                     0.01,
                     inc,
                     psi,
                     H.Flow,
                     H.sampling_rate,
                     0.0,
                     0.0,
                     400,
                     1,
                     0,
                     1,
                     0,
                     0,
                     0,
                     0)

    from pylab import *
    fig = figure()
    ax = fig.add_subplot(111)
#    ax.plot(h[:,0])
    hp = noise.resize_time_series(h[:,0],H.segment_length)
    hc = noise.resize_time_series(h[:,1],H.segment_length)
#    ax.plot(H.Times,hp)
    # roll the array so that the peak of the waveform is 1s from the end of the frame
#    # find the index corresponding to the trigger time
    index_trigtime = np.argmax(np.abs(hp))
#    
#    # Starting time for the signal chunk
#    # We want the trigger time 1s before the end of the segment
    index_wf_start = -(index_trigtime - int(H.sampling_rate*(H.T-1)))

    hp = np.roll(hp,index_wf_start)
    hc = np.roll(hp,index_wf_start)
    # window the data
    padding = 0.4
    window=tukey(H.segment_length,2.0*H.sampling_rate*padding/H.T)
    hp*=window
    hc*=window
    hptilde = np.fft.rfft(hp)
    hctilde = np.fft.rfft(hc)
    
    tc = H.trigtime+np.linspace(-0.05,0.05,1001)

#    ax.plot(H.Frequency[H.kmin:H.kmax],hptilde[H.kmin:H.kmax])
#    ax.plot(H.Frequency[H.kmin:H.kmax],H.FrequencySeries[H.kmin:H.kmax])
#
#    plt.show()
#    exit()

#    fig = figure()
#    ax = fig.add_subplot(111)
    logL = np.array([H.logLikelihood(hptilde, hctilde, ra, dec, psi, t)+L.logLikelihood(hptilde, hctilde, ra, dec, psi, t) for t in tc])
    C = ax.plot(tc,logL)
    print "%.15f"%tc[logL.argmax()]
    show()
                                



