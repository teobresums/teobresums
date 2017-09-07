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
        self.InversePowerSpectralDensity = 1./(self.PowerSpectralDensity*0.5/self.dt)
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
        
        timeShift = (self.Epoch-1)-tc + TimeDelayFromEarthCenter(self.location, ra, dec, tc)

        return np.exp(-1j*2.0*np.pi*timeShift*self.Frequency[self.kmin:self.kmax])*(fp*hptilde[self.kmin:self.kmax]+fc*hctilde[self.kmin:self.kmax])

    def logLikelihood(self, hptilde, hctilde, ra, dec, psi, tc):
        
        residuals = self.FrequencySeries[self.kmin:self.kmax] - self.Project(hptilde, hctilde, ra, dec, psi, tc)
        
        numerator = residuals*np.conj(residuals)

        return -np.sum(np.real(numerator)*self.InversePowerSpectralDensity[self.kmin:self.kmax]*self.dt)

if __name__ == "__main__":
    H = GravitationalWaveDetector('H1','data/H-H1_LOSC_4_V1-1126259446-32.txt', trigtime = 1126259462.43)
    L = GravitationalWaveDetector('L1','data/L-L1_LOSC_4_V1-1126259446-32.txt', trigtime = 1126259462.43)
    
    for attr, value in H.__dict__.iteritems():
        print attr, value

    from pyTEOBResumS import pyTEOBResumS

    ra = 1.82161142311
    dec = -1.27626704919
    tc = 1126259462.43
    m1=40.0
    m2=40.0
    spin1x = 0.0
    spin1y = 0.0
    spin1z = 0.001
    spin2x = 0.0
    spin2y = 0.0
    spin2z = 0.001
    inclination = 2.92033171962
    polarisation = 2.31550047039
    f_min = 20.0
    sampling_rate = 4096.
    dt = 1./sampling_rate
    LambdaAl2 = 0.0
    LambdaBl2 = 0.0
    LambdaAl3 = 0.0
    LambdaBl3 = 0.0
    LambdaAl4 = 0.0
    LambdaBl4 = 0.0
    distance = 1000.0

    flags ={'NQC':'1',
            'tidal':0,
            'speedy':1,
            'dynamics':0,
            'solver_scheme':0,
            'RWZ':0,
            'Yagi_fits':1,
            'spin':1,
            'multipoles':0,
            'geometric_units':0,
            'set':0
        }

    lm = -1

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
                     dt,
                     LambdaAl2,
                     LambdaAl3,
                     LambdaAl4,
                     LambdaBl2,
                     LambdaBl3,
                     LambdaBl4,
                     distance,
                     lm,
                     flags)

    from pylab import *
    fig = figure()
    ax = fig.add_subplot(111)
#    ax.loglog(H.Frequency,np.sqrt(H.PowerSpectralDensity))
#    show()

    # window the waveform
    padding = 0.5
    window=tukey(H.segment_length,padding)
    windowNorm = H.segment_length/np.sum(window**2)
#    h[:,0]*=window
#    ax.plot(h[:,0])
#    show()
#    exit()
#    h[:,1]*=window
    hp = noise.resize_time_series(h[:,0],H.segment_length)
    hc = noise.resize_time_series(h[:,1],H.segment_length)
    #    # Starting time for the signal chunk
    #    # We want the trigger time 1s before the end of the segment
    # find the index corresponding to the trigger time
#    index_trigtime = np.argmax(np.abs(hp-1j*hc))
#    index_wf_start = -(index_trigtime - int(H.sampling_rate*(H.T-1)))
#    hp = np.roll(hp,index_wf_start)
#    hc = np.roll(hp,index_wf_start)
    hp*=window
    hc*=window
    
#    ax.plot(hp)
#    ax.plot(H.Times,H.TimeSeries,alpha=0.5)
#    ax.plot(H.Times,hp)
#    ax.axvline(H.trigtime)
    # roll the array so that the peak of the waveform is 1s from the end of the frame
#
    hptilde = np.fft.rfft(hp)*windowNorm
    hctilde = np.fft.rfft(hc)*windowNorm
#    show()
#    exit()
    tc = H.trigtime+np.linspace(-0.05,0.05,1001)
#*
#    print len(np.exp(-1j*2.0*np.pi*((H.Epoch-1)-tc+0.007)*H.Frequency[H.kmin:H.kmax]))
#    print H.Epoch
#    ax.plot(H.Frequency[H.kmin:H.kmax],H.FrequencySeries[H.kmin:H.kmax],alpha=0.5)
#    ax.plot(H.Frequency[H.kmin:H.kmax],hptilde[H.kmin:H.kmax]*np.exp(-1j*2.0*np.pi*((H.Epoch-1)-H.trigtime-0.003)*H.Frequency[H.kmin:H.kmax]))
#    plt.xlim(20,300)
#    plt.ylim(-1e-19,1e-19)
#    plt.show()
#    exit()

#    fig = figure()
#    ax = fig.add_subplot(111)
    logL = np.array([H.logLikelihood(hptilde, hctilde, ra, dec, polarisation, t)+L.logLikelihood(hptilde, hctilde, ra, dec, polarisation, t) for t in tc])
    C = ax.plot(tc,logL)
    print "%.15f"%tc[logL.argmax()]
    show()
                                



