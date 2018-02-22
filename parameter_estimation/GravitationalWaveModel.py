#!/usr/bin/env python
import numpy as np
import cpnest.model
import sys
import os
from optparse import OptionParser
import itertools as it
import noise
from scipy.signal import tukey

from detector import GravitationalWaveDetector
from pyTEOBResumS import pyTEOBResumS

import lalsimulation as lalsim

detMap = {'H1': 'LHO_4k',
          'H2': 'LHO_2k',
          'L1': 'LLO_4k',
          'G1': 'GEO_600',
          'V1': 'VIRGO',
          'T1': 'TAMA_300'}

def McQ2Masses(mc, q):
    """
    Simple utility to convert between mc an d q to component masses
    """
    factor = mc * np.power(1. + q, 1.0/5.0);
    m1 = factor * np.power(q, -3.0/5.0);
    m2 = factor * np.power(q, +2.0/5.0);
    return m1, m2

def Masses2McQ(m1, m2):
    """
    Simple utility to convert between component masses to mc q
    """
    q = m2/m1
    eta = m1*m2/(m1+m2)
    mc = (m1*m2)**(3./5.)/(m1+m2)**(1./5.)
    return mc, q

def PolarToCartesian(a, th, ph):
    """
    Simple utility function to convert between polar and cartesian representations
    for the spin quantities
    """
    return a*np.sin(th)*np.cos(ph), a*np.sin(th)*np.sin(ph), a*np.cos(th)

class GravitationalWaveModel(cpnest.model.Model):
    """
    Gravitational Wave signal model
    """

    names = []
    bounds = []

    def __init__(self,
                 detector_names,
                 datafiles = None,
                 psd_files = None,
                 injection = False,
                 zero_noise=False,
                 sampling_rate=4096.,
                 T=4.0,
                 starttime = 1126259446.,
                 trigtime=1126259462.423,
                 template='LAL',
                 flow = 20,
                 fhigh = 500,
                 **kwargs):

        super(GravitationalWaveModel,self).__init__(**kwargs)
        
        self.detectors = detector_names
        self.T = T
        self.epoch = starttime
        self.zero_noise = zero_noise
        if datafiles is None: self.datafiles = len(self.detectors)*[None]
        else: self.datafiles = datafiles
        if psd_files is None: self.psd_files = len(self.detectors)*[None]
        else: self.psd_files = psd_files
        self.trigtime = trigtime
        self.template = template
        self.injection=injection
        self.sampling_rate = sampling_rate
        self.dt = 1./self.sampling_rate
        self.flow = flow
        self.fhigh = fhigh
        self.segment_length = int(self.sampling_rate*self.T)
        
        self.detectors = [GravitationalWaveDetector(name,
                                                    datafile = datum,
                                                    psd_file = psd_file,
                                                    T=self.T,
                                                    starttime = self.epoch,
                                                    trigtime = self.trigtime,
                                                    flow = self.flow,
                                                    fhigh = self.fhigh,
                                                    zero_noise = self.zero_noise,
                                                    **kwargs) for name,datum,psd_file in zip(self.detectors,self.datafiles,self.psd_files)]

        self.df = self.detectors[0].df
        self.padding = 0.1
        self.window=tukey(self.segment_length,self.padding)
        self.windowNorm = self.segment_length/np.sum(self.window**2)

        self.plain_template = None
        self.injected_template = None
        
        if self.injection:
            
            self.injection_parameters = {'mc':1.2,
                                        'q':1.0,
                                        'spin1':0.0,
                                        'theta_1l':0.0,
                                        'phi_1l':0.0,
                                        'spin2':0.0,
                                        'theta_2l':0.0,
                                        'phi_2l':0.0,
                                        'distance':100.0,
                                        'inclination':0.0,
                                        'ra':4.2,
                                        'dec':1.1,
                                        'psi':1.0,
                                        'phi0':0.0,
                                        'lambdaAl2':500.0,
                                        'lambdaBl2':500.0,
                                        'tc':1126259462.423}
            
            sys.stderr.write("Injection parameters:\n")
            for key, value in self.injection_parameters.iteritems():
                sys.stderr.write("%s --> %.3f\n"%(key,value))
            print "component masses are :", McQ2Masses(self.injection_parameters['mc'], self.injection_parameters['q'])
            if self.template == 'LAL':
                amp_order = 0
                phase_order = -1
                wave_flags = None
                non_GR_params = None
                approx = lalsim.IMRPhenomPv2
                
                spin1x, spin1y, spin1z = PolarToCartesian(self.injection_parameters['spin1'],self.injection_parameters['theta_1l'],self.injection_parameters['phi_1l'])
                spin2x, spin2y, spin2z = PolarToCartesian(self.injection_parameters['spin2'],self.injection_parameters['theta_2l'],self.injection_parameters['phi_2l'])
                m1,m2 = McQ2Masses(self.injection_parameters['mc'], self.injection_parameters['q'])
                hptilde, hctilde = lalsim.SimInspiralChooseFDWaveform(self.injection_parameters['phi0'],
                                   self.df,
                                   m1*lalsim.lal.MSUN_SI,
                                   m2*lalsim.lal.MSUN_SI,
                                   spin1x, spin1y, spin1z,
                                   spin2x, spin2y, spin2z,
                                   self.flow, self.fhigh, 100.0,
                                   self.injection_parameters['distance']*1e6*lalsim.lal.PC_SI,
                                   self.injection_parameters['inclination'],
                                   self.injection_parameters['lambdaAl2'], self.injection_parameters['lambdaBl2'],
                                   wave_flags, non_GR_params, amp_order, phase_order, approx)
                hptilde = hptilde.data.data
                hctilde = hctilde.data.data

            else:
                self.flags ={'NQC':'0',
                'tidal':1,
                'speedy':1,
                'dynamics':0,
                'solver_scheme':0,
                'RWZ':0,
                'Yagi_fits':1,
                'spin':0,
                'multipoles':0,
                'geometric_units':0,
                'set':0
                }
                self.injection_parameters['spin1z'] = self.injection_parameters['spin1']
                self.injection_parameters['spin2z'] = self.injection_parameters['spin2']
                self.injection_parameters['iota'] = self.injection_parameters['inclination']
                
                self.calculate_plain_template(self.injection_parameters)
                hptilde, hctilde = self.plain_template[0], self.plain_template[1]

            self.injected_template = []
            
            for d in self.detectors:
                self.injected_template.append(d.inject( hptilde, hctilde, self.injection_parameters['ra'], self.injection_parameters['dec'], self.injection_parameters['psi'], self.trigtime))

        self.logZnoise = self.log_nulllikelihood()
        
        #parameters
        if self.template == 'LAL':
            self.names=['phi0', 'ra', 'dec', 'tc', 'mc', 'q',
                        'iota', 'psi', 'distance',
                        'spin1','theta_1l','phi_1l',
                        'spin2','theta_2l','phi_2l']

            self.bounds=[[0,2.0*np.pi],
                         [0,2.0*np.pi],
                         [-np.pi/2.0,np.pi/2.0],
                         [self.trigtime-0.05,self.trigtime+0.05],
                         [10.0,50.0],
                         [0.1,1.0],
                         [0.0,np.pi],
                         [0.0,np.pi],
                         [1.0,1000.0],
                         [0.0,1.0],[-np.pi/2.0,np.pi/2.0],[0.0,2.0*np.pi],
                         [0.0,1.0],[-np.pi/2.0,np.pi/2.0],[0.0,2.0*np.pi]]
        else:
            self.names=['phi0', 'ra', 'dec', 'tc', 'mc', 'q',
                'iota', 'psi', 'distance','spin1z','spin2z','lambdaAl2','lambdaBl2']

            self.bounds=[[0,2.0*np.pi],
                         [0,2.0*np.pi],
                         [-np.pi/2.0,np.pi/2.0],
                         [self.trigtime-0.05,self.trigtime+0.05],
                         [1.0,2.0],
                         [0.5,1.0],
                         [-np.pi/2.,np.pi/2.],
                         [0.0,np.pi],
                         [1.0,500.0],
                         [-0.8,0.8],[-0.8,0.8],
                         [0.0,1000.0],[0.0,1000.0]]
#            self.bounds=[[0,2.0*np.pi],
#                         [4.1,4.3],
#                         [1.0,1.2],
#                         [self.trigtime-0.0005,self.trigtime+0.0005],
#                         [27.0,28.0],
#                         [0.89,0.91],
#                         [-0.0001,0.0001],
#                         [0.0,0.001],
#                         [999.0,1001.0],
#                         [-0.0001,0.0001],[-0.0001,0.0001]]
            self.flags ={'NQC':'1',
                'tidal':1,
                'speedy':1,
                'dynamics':0,
                'solver_scheme':0,
                'RWZ':0,
                'Yagi_fits':1,
                'spin':0,
                'multipoles':0,
                'geometric_units':0,
                'set':0
            }
#        self.padding = 0.1
#        self.window=tukey(self.segment_length,self.padding)
#        self.windowNorm = self.segment_length/np.sum(self.window**2)


    def calculate_plain_template(self,x):
        mc = x['mc']
        q = x['q']
        d = x['distance']
        m1, m2 = McQ2Masses(mc, q)
        
        if self.template == 'LAL':
            amp_order = 0
            phase_order = -1
            wave_flags = None
            non_GR_params = None
            approx = lalsim.IMRPhenomPv2
            spin1x, spin1y, spin1z = PolarToCartesian(x['spin1'], x['theta_1l'], x['phi_1l'])
            spin2x, spin2y, spin2z = PolarToCartesian(x['spin2'], x['theta_2l'], x['phi_2l'])
            hptilde, hctilde = lalsim.SimInspiralChooseFDWaveform(x['phi0'],
                                           self.df,
                                           m1*lalsim.lal.MSUN_SI,
                                           m2*lalsim.lal.MSUN_SI,
                                           spin1x, spin1y, spin1z,
                                           spin2x, spin2y, spin2z,
                                           self.flow, self.fhigh, 100.0,
                                           d*1e6*lalsim.lal.PC_SI,
                                           x['iota'],
                                           0.0, 0.0,
                                           wave_flags, non_GR_params, amp_order, phase_order, approx)

            hp = hptilde.data.data
            hc = hctilde.data.data
        
        else:
            h = pyTEOBResumS(m1,
                             m2,
                             0.0,
                             0.0,
                             x['spin1z'],
                             0.0,
                             0.0,
                             x['spin2z'],
                             x['iota'],
                             x['phi0'],
                             self.flow,
                             self.dt,
                             x['lambdaAl2'],
                             x['lambdaBl2'],
                             0.0,
                             0.0,
                             0.0,
                             0.0,
                             d,
                             -1,
                             self.flags)

            hp = noise.resize_time_series(h[:,0],self.segment_length)
            hc = noise.resize_time_series(h[:,1],self.segment_length)

            hp*=self.window
            hc*=self.window

            hp = np.fft.rfft(hp)*self.windowNorm*self.dt# we multiply by dt to get the dimensionfull FFT
            hc = np.fft.rfft(hc)*self.windowNorm*self.dt
        
        self.plain_template = (hp,hc)
    
    def log_nulllikelihood(self):
        
        Nfd = len(self.detectors[0].FrequencySeries)
        h = np.zeros(Nfd)
        
        return np.sum([d.logLikelihood(h, h, 0.0, 0.0, 0.0, 0.0) for d in self.detectors])
    
    def log_likelihood(self,x):
        
        self.calculate_plain_template(x)
        if np.any(np.isnan(self.plain_template)):
            return -np.inf
        logL = np.sum([d.logLikelihood(self.plain_template[0], self.plain_template[1], x['ra'], x['dec'], x['psi'], x['tc']) for d in self.detectors])

        return logL
    
    def log_prior(self, x):
        if np.isfinite(super(GravitationalWaveModel,self).log_prior(x)):
            logP = 2.0*np.log(x['distance'])
            logP += np.log(np.abs(np.cos(x['dec'])))
            logP += np.log(np.abs(np.cos(x['iota'])))
            mc = x['mc']
            q = x['q']
            m1, m2 = McQ2Masses(mc, q)
            logP += np.log(m1*m1/mc)
            if self.template == 'LAL': logP += np.log(np.abs(np.cos(x['theta_1l'])))+np.log(np.abs(np.cos(x['theta_2l'])))
            return logP
        else:
            return -np.inf

    def potential_energy(self, x):
        return -self.log_prior(x)

    def force(self, x):
        self.gradient = {}
        self.gradient['distance'] = 2.0/x['distance']
        self.gradient['dec'] = -np.tan(x['dec'])
        self.gradient['iota'] = -np.tan(x['iota'])
        mc = x['mc']
        q = x['q']
        m1, m2 = McQ2Masses(mc, q)
        m1_derivative = np.power(1. + q, 1.0/5.0)*np.power(q, -3.0/5.0)
        self.gradient['mc'] = -m1_derivative/mc
        self.gradient['q'] = -mc * (2.0*q+3)/(5.0*np.power(q,8.0/5.0)*np.power(1. + q, 5.0/5.0))
        if self.template == 'LAL':
            self.gradient['theta_1l'] = -np.tan(x['theta_1l'])
            self.gradient['theta_2l'] = -np.tan(x['theta_2l'])
        return self.gradient

if __name__=='__main__':
    parser=OptionParser()
    parser.add_option('-o','--out-dir',default=None,type='string',metavar='DIR',help='Directory for output: defaults to gw150914/')
    parser.add_option('-t','--threads',default=None,type='int',metavar='N',help='Number of threads (default = 1/core)')
    parser.add_option('-f','--full-run',default=0,type='int',metavar='full_run',help='perform a full PE run')
    parser.add_option('--inject',default=False,action='store_true',help='Inject signal')
    parser.add_option('--zero-noise',default=False,action='store_true',help='Generate a 0 noise realisation')
    parser.add_option('--template',default='LAL',type='str',metavar='template',help='template to use for the run')
    parser.add_option('--seglen',default=4,type='float',metavar='seglen',help='length of the data stretch to analyse')
    parser.add_option('--flow',default=20,type='float',metavar='flow',help='low frequency cutoff')
    parser.add_option('--fhigh',default=500,type='float',metavar='fhigh',help='high frequency cutoff')
    parser.add_option('--nlive',default=1024,type='int',metavar='n',help='Live points')
    parser.add_option('--maxmcmc',default=1024,type='int',metavar='m',help='max MCMC points')
    parser.add_option('--poolsize',default=1000,type='int',metavar='k',help='number of points in the ensemble sampler pool')
    (opts,args)=parser.parse_args()

    trigtime = 1126259462.423
    starttime = trigtime - opts.seglen + 2
    
    if opts.out_dir is None:
        opts.out_dir='./gw150914/'

    if opts.full_run:
        signal_model = GravitationalWaveModel(['H1','L1','V1'],
                                              T         = opts.seglen,
                                              template  = opts.template,
                                              sampling_rate = 4096.,
                                              injection = opts.inject,
                                              zero_noise    = opts.zero_noise,
                                              starttime = starttime,
                                              trigtime  = trigtime,
#                                              psd_files = ['H-GW15-asd.txt',
#                                                           'L-GW15-asd.txt'],
                                              psd_files = ['/Users/wdp/src/lalsuite/lalsimulation/src/LIGO-P1200087-v18-AdV_DESIGN.txt',
                                                           '/Users/wdp/src/lalsuite/lalsimulation/src/LIGO-P1200087-v18-AdV_DESIGN.txt',
                                                           '/Users/wdp/src/lalsuite/lalsimulation/src/LIGO-P1200087-v18-AdV_DESIGN.txt'],
#                                                               '/Users/wdp/src/lalsuite/lalsimulation/src/LIGO-P1200087-v18-AdV_DESIGN.txt'],
#                                              datafiles = ['data/H-H1_LOSC_4_V1-1126259446-32.txt','data/L-L1_LOSC_4_V1-1126259446-32.txt'],
                                              flow=opts.flow,
                                              fhigh=opts.fhigh)
        print('Noise evidence {0}'.format(signal_model.logZnoise))

        work=cpnest.CPNest(signal_model,
                           verbose=3,
                           Poolsize=opts.poolsize,
                           Nthreads=opts.threads,
                           Nlive=opts.nlive,
                           maxmcmc=opts.maxmcmc,
                           output=opts.out_dir,
                           balance_samplers=False)
        work.run()
        print('Signal evidence {0}'.format(work.NS.logZ))
        logB = work.NS.logZ-signal_model.logZnoise
        print('log B {0}'.format(logB))
    else:
        signal_model = GravitationalWaveModel(['H1','L1','V1'],
                                              T=opts.seglen,
                                              template = opts.template,
                                              sampling_rate = 2048.,
                                              injection = opts.inject,
                                              zero_noise = opts.zero_noise,
                                              starttime = 1126259459.423,
                                              trigtime = 1126259462.423,
                                              psd_files = ['/Users/wdp/src/lalsuite/lalsimulation/src/LIGO-P1200087-v18-AdV_DESIGN.txt',
                                                           '/Users/wdp/src/lalsuite/lalsimulation/src/LIGO-P1200087-v18-AdV_DESIGN.txt',
                                                           '/Users/wdp/src/lalsuite/lalsimulation/src/LIGO-P1200087-v18-AdV_DESIGN.txt'],
#                                              datafiles = ['data/H-H1_LOSC_4_V1-1126259446-32.txt','data/L-L1_LOSC_4_V1-1126259446-32.txt'],
                                              flow=opts.flow,
                                              fhigh=opts.fhigh)
        x = np.genfromtxt(os.path.join(opts.out_dir,'posterior.dat'),names=True)
        logZ = np.loadtxt(os.path.join(opts.out_dir,'chain_{0}_1234.txt_evidence.txt'.format(opts.nlive)))[0]
        print('Signal evidence {0}'.format(logZ))
        logB = logZ-signal_model.logZnoise
        print('log B {0}'.format(logB))
    import corner
    import matplotlib.pyplot as plt
    intrinsic = ['mc','q','spin1','spin2']

    figure = corner.corner(np.array([x[n] for n in intrinsic]).T, labels=[r"$\mathcal{M}$", r"$q$", r"$s_{1}$", r"$s_{2}$"],
                       quantiles=[0.16, 0.5, 0.84], truths = [signal_model.injection_parameters[n] for n in intrinsic],
                       show_titles=True, title_kwargs={"fontsize": 12}, smooth1d=0.5)
    plt.savefig(os.path.join(opts.out_dir,'intrinsic.pdf'),bbbox_inches='tight')
    extrinsic = ['ra','dec','distance','tc']

    figure = corner.corner(np.array([x[n] for n in extrinsic]).T, labels=[r"$\mathrm{RA}$", r"$\mathrm{dec}$", r"$D_L/Mpc$", r"$t_c/s$"],
               quantiles=[0.16, 0.5, 0.84], truths = [signal_model.injection_parameters[n] for n in extrinsic],
               show_titles=True, title_kwargs={"fontsize": 12}, smooth1d=0.5)
    plt.savefig(os.path.join(opts.out_dir,'extrinsic.pdf'),bbbox_inches='tight')
