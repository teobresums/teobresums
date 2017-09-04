from __future__ import division
import numpy as np
import noise
from pylal import inject
from lal import ComputeDetAMResponse, GreenwichMeanSiderealTime, TimeDelayFromEarthCenter

detMap = {'H1': 'LHO_4k', 'H2': 'LHO_2k', 'L1': 'LLO_4k', 'G1': 'GEO_600', 'V1': 'VIRGO', 'T1': 'TAMA_300'}

class GravitationalWaveDetector(object):
    
    def __init__(self, name, datafile, **kwargs):
        self.name = name
        self.lal_detector = inject.cached_detector[detMap[self.name]]
        self.location = self.lal_detector.location
        self.Times, self.TimeSeries, self.Frequency, self.FrequencySeries, self.PowerSpectralDensity = noise.load_data(datafile, **kwargs)
        self.Epoch = self.Times[0]
        self.dt = np.diff(self.Times)[0]
        self.T = self.Times[-1]-self.Times[0]
        self.df = np.diff(self.Frequency)[0]
        self.sampling_rate = 1./self.dt
        self.segment_length = int(self.T*self.sampling_rate)
        self.InversePowerSpectralDensity = 1./self.PowerSpectralDensity
        self.Flow = 20.0
        self.Fhigh = 300.0
        self.kmin = np.argmin(np.abs(self.Frequency-self.Flow))
        self.kmax = np.argmin(np.abs(self.Frequency-self.Fhigh))
    
    def Project(self, hptilde, hctilde, ra, dec, psi, tc):
        """
        projects and timeshifts the GW signal onto the detector
        """
        gmst = GreenwichMeanSiderealTime(tc)
        fp,fc = ComputeDetAMResponse(self.lal_detector.response, ra, dec, psi, gmst)
        timeShift = self.Epoch-tc + TimeDelayFromEarthCenter(self.location, ra, dec, tc)
        return np.exp(-1j*2.0*np.pi*timeShift)*(fp*hptilde+fc*hctilde)

    def logLikelihood(self, hp, hc, ra, dec, psi, tc, domain = 'time'):

        if domain == 'time':
            f, hp = noise.fd_from_td(t, hp, srate = self.sampling_rate, N = self.segment_length)
            f, hc = noise.fd_from_td(t, hc, srate = self.sampling_rate, N = self.segment_length)
        
        TwoDeltaTOverN = 2.0*self.dt/self.segment_length
        residuals = self.FrequencySeries - self.Project(hp, hc, ra, dec, psi, tc)
        numerator = residuals[self.kmin:self.kmax]*np.conj(residuals[self.kmin:self.kmax])
        return -TwoDeltaTOverN*np.sum(np.real(numerator)*(self.InversePowerSpectralDensity[self.kmin:self.kmax]/(self.dt*self.dt)))

if __name__ == "__main__":
    H = GravitationalWaveDetector('H1','data/H-H1_LOSC_4_V1-1126259446-32.txt')
    L = GravitationalWaveDetector('L1','data/L-L1_LOSC_4_V1-1126259446-32.txt')
    
    for attr, value in H.__dict__.iteritems():
        print attr, value

    from pyTEOBResumS import pyTEOBResumS
    
    f_min = 20.0
    sampling_rate = 4096.
    segment_length = int(sampling_rate*32)
    t = np.linspace(0.0,32,segment_length)
    
    tc = 1126259462.423
    ra = 0.0
    dec = 0.0
    inc = 0.0
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
    ra = np.linspace(0.0,2.0*np.pi,32)
    dec = np.linspace(-np.pi/2.0,np.pi/2.0,32)
    
    fig = figure()
    ax = fig.add_subplot(111,projection='aitoff')
    C = ax.contourf([[H.logLikelihood(h[:,0], h[:,1], r, d, psi, tc)+L.logLikelihood(h[:,0], h[:,1], r, d, psi, tc) for r in ra] for d in dec],10)
    colorbar(C)
    show()
                                



