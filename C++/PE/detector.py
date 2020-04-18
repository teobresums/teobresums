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
                 fhigh      = None,
                 zero_noise = False,
                 **kwargs):
        
        self.name       = name
        self.datafile   = datafile
        self.psd_file   = psd_file
        self.Epoch      = starttime
        self.sampling_rate = sampling_rate
        self.flow       = flow
        if fhigh is None:
            self.fhigh  = self.sampling_rate/2.
        else:
            self.fhigh  = fhigh
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
        
        self.df     = 1./self.T
        self.dt     = 1.0/self.sampling_rate
        self.segment_length = int(self.T*self.sampling_rate)
        self.kmin   = int(self.flow/self.df)
        self.kmax   = int(self.fhigh/self.df)+1
        self.sigmasq =  self.PowerSpectralDensity * self.dt * self.dt
        self.TwoDeltaTOverN = 2.0*self.dt/float(self.segment_length)
    
    def Project(self, hptilde, hctilde, ra, dec, psi, tc, domain = "F"):
        """
        projects and timeshifts the GW signal onto the detector
        """
        timeShift = 0.0
        tgps  = LIGOTimeGPS(tc)
        gmst  = GreenwichMeanSiderealTime(tgps)
        fp,fc = ComputeDetAMResponse(self.lal_detector.response, ra, dec, psi, gmst)
        timedelay = TimeDelayFromEarthCenter(self.lal_detector.location, ra, dec, tgps)
        if domain == "F":
            timeShift = self.Epoch - tc
        timeShift += timedelay
        #np.exp(1j*2.0*np.pi*timeShift*self.Frequency[self.kmin:self.kmax])*
        return (fp*hptilde[self.kmin:self.kmax]+fc*hctilde[self.kmin:self.kmax])*np.exp(1j*2.0*np.pi*timeShift*self.Frequency[self.kmin:self.kmax])

    def logLikelihood(self, hptilde, hctilde, ra, dec, psi, tc, domain = "F"):
        
        template    = self.Project(hptilde, hctilde, ra, dec, psi, tc, domain =  domain)
        data        = self.FrequencySeries[self.kmin:self.kmax]
        residuals   = (data - template)
        overlap     = np.conj(residuals)*residuals/self.sigmasq[self.kmin:self.kmax]

        return -self.TwoDeltaTOverN*np.sum(overlap).real

    def inject(self, hptilde, hctilde, ra, dec, psi, tc, domain = "F"):
        template = self.Project(hptilde, hctilde, ra, dec, psi, tc, domain = domain)
        self.FrequencySeries[self.kmin:self.kmax] += template
        self.SNR = np.sqrt(4.0*self.df*np.sum(np.conj(template)*template/self.PowerSpectralDensity[self.kmin:self.kmax]).real)
        print("Injected SNR in %s = %.2f\n"%(self.name,self.SNR))
        return template

if __name__ == "__main__":
    
#    starttime = 1126259460.0
#    trigtime  = 1126259462.43

    starttime = 0.0
    trigtime  = 2.0

    H = GravitationalWaveDetector('H1',
                                  psd_file = '/Users/wdp/src/lalsuite/lalsimulation/src/LIGO-T0900288-v3-ZERO_DET_high_P.txt',
                                  trigtime = trigtime,
                                  starttime= starttime,
                                  T = 4.0,
                                  flow = 30.,
                                  sampling_rate = 4096.,
                                  zero_noise = True)

    import lalsimulation as lalsim
    from pylab import *
    from pyTEOBResumS import pyTEOBResumS

    amp_order = 0
    phase_order = -1
    wave_flags = None
    non_GR_params = None
    tcs = np.linspace(trigtime-0.1,trigtime+0.1,101)
    logl = []
    logl1 = []
    # waveform parameters
    
    ra = 0.0
    dec = 0.0
    m1 = 32.0
    m2 = 27.0
    spin1x = 0.0
    spin1y = 0.0
    spin1z = 0.0
    spin2x = 0.0
    spin2y = 0.0
    spin2z = 0.0
    inclination = 0.0
    polarisation = 0.0
    phase = 0.0
    f_min = H.flow-5.0
    sampling_rate = H.sampling_rate
    dt = H.dt
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
            'Yagi_fits':0,
            'spin':1,
            'multipoles':0,
            'geometric_units':0,
            'set':0
        }

    lm = -1

    # window the waveform
    padding = 0.1
    window  = tukey(H.segment_length,padding)
    windowNorm = np.sqrt(H.segment_length/np.sum(window**2))

    approx = lalsim.IMRPhenomPv2

    # generate the injection
    if 0:
#        f = plt.figure()
#        ax = f.add_subplot(111)
        hp_imr_inj, hc_imr_inj = lalsim.SimInspiralChooseFDWaveform(phase,
                           H.df,
                           m1*lalsim.lal.MSUN_SI,
                           m2*lalsim.lal.MSUN_SI,
                           spin1x, spin1y, spin1z,
                           spin2x, spin2y, spin2z,
                           H.flow, H.sampling_rate/2., 0.0,
                           distance*1e6*lalsim.lal.PC_SI,
                           inclination,
                           0.0, 0.0,
                           wave_flags, non_GR_params, amp_order, phase_order, approx)

        injection_imr = H.inject(hp_imr_inj.data.data, hc_imr_inj.data.data, ra, dec, polarisation, trigtime)
#        ax.plot(H.Frequency[H.kmin:H.kmax],injection_imr, label='IMRPhenomPv2', lw=3)
    else:
        h_inj = pyTEOBResumS(m1,
                         m2,
                         spin1x,
                         spin1y,
                         spin1z,
                         spin2x,
                         spin2y,
                         spin2z,
                         inclination,
                         phase,
                         H.flow,
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
            
        hp_inj, hc_inj = noise.resize_time_series(h_inj, H.segment_length, dt, starttime, trigtime)

        hp_inj *= window
        hc_inj *= window

        hptilde_inj = (np.fft.rfft(hp_inj)*windowNorm*H.dt)
        hctilde_inj = (np.fft.rfft(hc_inj)*windowNorm*H.dt)
        H.FrequencySeries[:] = 0.0+1j*0.0
        injection = H.inject(hptilde_inj, hctilde_inj, ra, dec, polarisation, trigtime, domain = "T")
#        ax.plot(H.Frequency[H.kmin:H.kmax],injection,label='TEOB')
#    plt.legend()
#    plt.show()
#    exit()
    for k,tc in enumerate(tcs):

        h = pyTEOBResumS(m1,
                         m2,
                         spin1x,
                         spin1y,
                         spin1z,
                         spin2x,
                         spin2y,
                         spin2z,
                         inclination,
                         phase,
                         H.flow,
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

        hp, hc = noise.resize_time_series(h, H.segment_length, dt, starttime, tc)

        hp *= window
        hc *= window

        hptilde = (np.fft.rfft(hp)*windowNorm*H.dt)
        hctilde = (np.fft.rfft(hc)*windowNorm*H.dt)

        teob = H.Project(hptilde, hctilde, ra, dec, polarisation, tc)

        hptilde1, hctilde1 = lalsim.SimInspiralChooseFDWaveform(phase,
                           H.df,
                           m1*lalsim.lal.MSUN_SI,
                           m2*lalsim.lal.MSUN_SI,
                           spin1x, spin1y, spin1z,
                           spin2x, spin2y, spin2z,
                           H.flow, H.sampling_rate/2., 0.0,
                           distance*1e6*lalsim.lal.PC_SI,
                           inclination,
                           0.0, 0.0,
                           wave_flags, non_GR_params, amp_order, phase_order, approx)
                           
        logl.append(H.logLikelihood(hptilde, hctilde, ra, dec, polarisation, tc, domain = "T"))
        logl1.append(H.logLikelihood(hptilde1.data.data, hctilde1.data.data, ra, dec, polarisation, tc))
        print k,logl[-1],logl1[-1]
    plt.plot(tcs,logl)
    plt.plot(tcs,logl1)
    print tcs[np.argmax(logl)], tcs[np.argmax(logl1)], (tcs[np.argmax(logl)]-trigtime)/np.diff(tcs)[0]
    plt.axvline(trigtime,color='k',linestyle='dashed')
    plt.show()



