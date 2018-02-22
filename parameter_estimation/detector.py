from __future__ import division
import numpy as np
import noise
from pylal import inject
from lal import ComputeDetAMResponse, GreenwichMeanSiderealTime, TimeDelayFromEarthCenter, LIGOTimeGPS
from scipy.signal import tukey

detMap = {'H1': 'LHO_4k', 'H2': 'LHO_2k', 'L1': 'LLO_4k', 'G1': 'GEO_600', 'V1': 'VIRGO', 'T1': 'TAMA_300'}

class GravitationalWaveDetector(object):

    def __init__(self, name,
                 datafile   = None,
                 psd_file   = None,
                 T          = 4.0,
                 starttime  = 1126259446.,
                 trigtime   = 1126259462.423,
                 sampling_rate = 4096.,
                 flow       = 20,
                 fhigh      = 1024,
                 zero_noise = False,
                 **kwargs):
        
        self.name       = name
        self.datafile   = datafile
        self.psd_file   = psd_file
        self.Epoch      = starttime
        self.sampling_rate = sampling_rate
        self.flow       = flow
        self.fhigh      = fhigh
        self.trigtime   = trigtime
        self.zero_noise = zero_noise
        self.T          = T
        
        self.lal_detector = inject.cached_detector[detMap[self.name]]
        
        if datafile is not None:
            self.Times, self.TimeSeries, self.Frequency, self.FrequencySeries, self.PowerSpectralDensity = noise.load_data(self.datafile, chunk_size=self.T, trigtime=self.trigtime, psd_file=self.psd_file, **kwargs)
        elif psd_file is not None:
            self.Times, self.Frequency, self.FrequencySeries, self.PowerSpectralDensity = noise.generate_data(self.psd_file, T=self.T, starttime=self.Epoch, sampling_rate=self.sampling_rate, zero_noise = self.zero_noise, **kwargs)
        else:
            print("User must specify either a datafile or a PSD file")
            exit(-1)
        
        self.df     = np.diff(self.Frequency)[0]
        self.dt     = 1.0/self.sampling_rate
        self.segment_length = int(self.T*self.sampling_rate)
        self.kmin   = int(self.flow/self.df)
        self.kmax   = int(self.fhigh/self.df)
        self.sigmasq =  self.PowerSpectralDensity * self.dt * self.dt
        self.TwoDeltaTOverN = 2.0*self.dt/float(self.segment_length)
    
    def Project(self, hptilde, hctilde, ra, dec, psi, tc):
        """
        projects and timeshifts the GW signal onto the detector
        """
        tgps  = LIGOTimeGPS(tc)
        gmst  = GreenwichMeanSiderealTime(tgps)
        fp,fc = ComputeDetAMResponse(self.lal_detector.response, ra, dec, psi, gmst)
        timeShift = (tc - self.Epoch) + TimeDelayFromEarthCenter(self.lal_detector.location, ra, dec, tgps)

        return np.exp(-1j*2.0*np.pi*timeShift*self.Frequency)*(fp*hptilde+fc*hctilde)

    def logLikelihood(self, hptilde, hctilde, ra, dec, psi, tc):
        
        template    = self.Project(hptilde, hctilde, ra, dec, psi, tc)
        data        = self.FrequencySeries[self.kmin:self.kmax]
        residuals   = (data - template)
        overlap     = np.conj(residuals)*residuals/self.sigmasq[self.kmin:self.kmax]

        return -self.TwoDeltaTOverN*np.sum(overlap).real

    def inject(self, hptilde, hctilde, ra, dec, psi, tc):
        template = self.Project(hptilde, hctilde, ra, dec, psi, tc)
        self.FrequencySeries[self.kmin:self.kmax] += template
        self.SNR = np.sqrt(4.0*self.df*np.sum(np.conj(template)*template/self.PowerSpectralDensity[self.kmin:self.kmax]).real)
        print("Injected SNR in %s = %.2f\n"%(self.name,self.SNR))
        return template

if __name__ == "__main__":
    H = GravitationalWaveDetector('H1',
                                  psd_file = '/Users/wdp/src/lalsuite/lalsimulation/src/LIGO-T0900288-v3-ZERO_DET_high_P.txt',
                                  trigtime = 1126259462.43,
                                  starttime= 1126259459.43,
                                  T = 8.,
                                  flow = 20.,
                                  fhigh = 1024.,
                                  sampling_rate = 2048.,
                                  zero_noise = False)

    from pylab import *

    fig = figure()
    ax = fig.add_subplot(111)
    C = ax.plot(H.Frequency,H.FrequencySeries,alpha=0.5)

    from pyTEOBResumS import pyTEOBResumS

    ra = 1.82161142311
    dec = -1.27626704919
    tc = 1126259462.43
    m1 = 40.0
    m2 = 40.0
    spin1x = 0.0
    spin1y = 0.0
    spin1z = 0.1
    spin2x = 0.0
    spin2y = 0.0
    spin2z = 0.7
    inclination = 0.0
    polarisation = 0.0
    f_min = H.flow
    sampling_rate = H.sampling_rate
    dt = H.dt
    LambdaAl2 = 0.0
    LambdaBl2 = 0.0
    LambdaAl3 = 0.0
    LambdaBl3 = 0.0
    LambdaAl4 = 0.0
    LambdaBl4 = 0.0
    distance = 100.0

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

    lm = 1

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

    # window the waveform
    padding = 0.1
    window  = tukey(H.segment_length,padding)
    windowNorm = H.segment_length/np.sum(window**2)

    hp = noise.resize_time_series(h[:,0], H.segment_length)
    hc = noise.resize_time_series(h[:,1], H.segment_length)

    hp*=window
    hc*=window

    hptilde = (np.fft.rfft(hp)*windowNorm*H.dt)[1:]
    hctilde = (np.fft.rfft(hc)*windowNorm*H.dt)[1:]

    teob = H.Project(hptilde, hctilde, ra, dec, polarisation, tc)
    import lalsimulation as lalsim
    amp_order = 0
    phase_order = -1
    wave_flags = None
    non_GR_params = None
    approx = lalsim.IMRPhenomPv2

    hptilde1, hctilde1 = lalsim.SimInspiralChooseFDWaveform(3.14/2.,
                       H.df,
                       m1*lalsim.lal.MSUN_SI,
                       m2*lalsim.lal.MSUN_SI,
                       spin1x, spin1y, spin1z,
                       spin2x, spin2y, spin2z,
                       H.flow, H.sampling_rate/2, 100.0,
                       distance*1e6*lalsim.lal.PC_SI,
                       0.0,
                       0.0, 0.0,
                       wave_flags, non_GR_params, amp_order, phase_order, approx)

    imr = H.Project(hptilde1.data.data[1:], hctilde1.data.data[1:], ra, dec, polarisation, tc)
    C = ax.plot(H.Frequency,hptilde,alpha=0.5,color='k')
    C = ax.plot(H.Frequency,hptilde1.data.data[1:],alpha=0.5,color='r')

    show()
                                



