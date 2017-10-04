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
    factor = mc * np.power(1. + q, 1.0/5.0);
    m1 = factor * np.power(q, -3.0/5.0);
    m2 = factor * np.power(q, +2.0/5.0);
    return m1, m2

class GravitationalWaveModel(cpnest.model.Model):
    
    names = []
    bounds = []

    def __init__(self,
                 inject=False,
                 chunk_size=8.0,
                 trigtime=1126259462.423,
                 template='LAL',
                 flow = 20,
                 fhigh = 500,
                 **kwargs):
        
        super(GravitationalWaveModel,self).__init__(**kwargs)
        # this is the merger time in H1
        self.tevent = trigtime
        self.template = template
        self.inject=inject
        self.flow = flow
        self.fhigh = fhigh
        
        self.detectors = [GravitationalWaveDetector('H1','data/H-H1_LOSC_4_V1-1126259446-32.txt', chunk_size=chunk_size, trigtime=trigtime, flow = self.flow, fhigh = self.fhigh, **kwargs),
                          GravitationalWaveDetector('L1','data/L-L1_LOSC_4_V1-1126259446-32.txt', chunk_size=chunk_size, trigtime=trigtime, flow = self.flow, fhigh = self.fhigh, **kwargs)]

        self.sampling_rate = self.detectors[0].sampling_rate
        self.segment_length = self.detectors[0].segment_length
        self.dt = 1./self.sampling_rate
        self.df = self.detectors[0].df
        
        #parameters
        if self.template == 'LAL':
            self.names=['phi0', 'ra', 'dec', 'tc', 'mc', 'q',
                        'iota', 'psi', 'distance',
                        'spin1x','spin1y','spin1z',
                        'spin2x','spin2y','spin2z']

            self.bounds=[[0,2.0*np.pi],
                         [0,2.0*np.pi],
                         [-np.pi/2.0,np.pi/2.0],
                         [self.tevent-0.05,self.tevent+0.05],
                         [25.0,35.0],
                         [0.5,1.0],
                         [0.0,np.pi],
                         [0.0,np.pi],
                         [1.0,2000.0],
                         [-0.5,0.5],[-0.5,0.5],[-0.5,0.5],
                         [-0.5,0.5],[-0.5,0.5],[-0.5,0.5]]
        else:
            self.names=['phi0', 'ra', 'dec', 'tc', 'mc', 'q',
                'iota', 'psi', 'distance','spin1z','spin2z']

            self.bounds=[[0,2.0*np.pi],
                         [0,2.0*np.pi],
                         [-np.pi/2.0,np.pi/2.0],
                         [self.tevent-0.05,self.tevent+0.05],
                         [25.0,35.0],
                         [0.5,1.0],
                         [0.0,np.pi],
                         [0.0,np.pi],
                         [1.0,2000.0],
                         [-0.5,0.5],[-0.5,0.5]]
                         
            self.flags ={'NQC':'1',
                'tidal':0,
                'speedy':1,
                'dynamics':0,
                'solver_scheme':2,
                'RWZ':0,
                'Yagi_fits':0,
                'spin':0,
                'multipoles':0,
                'geometric_units':0,
                'set':0
            }
        self.window=tukey(self.segment_length,0.5)
        self.windowNorm = self.segment_length/np.sum(self.window**2)
            
    def log_likelihood(self,x, template = 'TEOB'):
        
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
            
#            hptilde, hctilde = lalsim.SimInspiralChooseFDWaveform(x['phi0'],
#                                           self.df,
#                                           m1*lalsim.lal.MSUN_SI,
#                                           m2*lalsim.lal.MSUN_SI,
#                                           0.0, 0.0, x['spin1z'],
#                                           0.0, 0.0, x['spin2z'],
#                                           self.flow, self.fhigh, 100.0,
#                                           d*1e6*lalsim.lal.PC_SI,
#                                           x['iota'],
#                                           0.0, 0.0,
#                                           wave_flags, non_GR_params, amp_order, phase_order, approx)
            hptilde, hctilde = lalsim.SimInspiralChooseFDWaveform(x['phi0'],
                                           self.df,
                                           m1*lalsim.lal.MSUN_SI,
                                           m2*lalsim.lal.MSUN_SI,
                                           x['spin1x'], x['spin1y'], x['spin1z'],
                                           x['spin2x'], x['spin2y'], x['spin2z'],
                                           self.flow, self.fhigh, 100.0,
                                           d*1e6*lalsim.lal.PC_SI,
                                           x['iota'],
                                           0.0, 0.0,
                                           wave_flags, non_GR_params, amp_order, phase_order, approx)
            hp = noise.resize_time_series(hptilde.data.data,self.segment_length/2+1)
            hc = noise.resize_time_series(hctilde.data.data,self.segment_length/2+1)
                
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
                             x['psi'],
                             self.flow,
                             self.dt,
                             0.0,
                             0.0,
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

        return np.sum([d.logLikelihood(hp, hc, x['ra'], x['dec'], x['psi'], x['tc']) for d in self.detectors])
    
    def log_prior(self, x):
        if np.isfinite(super(GravitationalWaveModel,self).log_prior(x)):
            logP = 2.0*np.log(x['distance'])
            logP += np.log(np.abs(np.cos(x['dec'])))
            logP += np.log(np.abs(np.cos(x['iota'])))
            mc = x['mc']
            q = x['q']
            m1, m2 = McQ2Masses(mc, q)
            logP += np.log(m1*m1/mc)
            return logP
        else:
            return -np.inf

class NoiseModel(cpnest.model.Model):
    
    names = []
    bounds = []

    def __init__(self, chunk_size=8.0, trigtime=1126259462.43, flow = 20, fhigh = 500, **kwargs):
        
        super(NoiseModel,self).__init__(**kwargs)
        # this is the merger time in H1
        self.tevent = 1126259462.423
        
        self.detectors = [GravitationalWaveDetector('H1','data/H-H1_LOSC_4_V1-1126259446-32.txt', chunk_size=chunk_size, trigtime=trigtime, flow = flow, fhigh = fhigh, **kwargs),
                          GravitationalWaveDetector('L1','data/L-L1_LOSC_4_V1-1126259446-32.txt', chunk_size=chunk_size, trigtime=trigtime, flow = flow, fhigh = fhigh, **kwargs)]

    def log_likelihood(self,x):
        
        Nfd = self.detectors[0].segment_length/2+1
        h = np.zeros(Nfd)
        
        return np.sum([d.logLikelihood(h, h, 0.0, 0.0, 0.0, 0.0) for d in self.detectors])

if __name__=='__main__':
    parser=OptionParser()
    parser.add_option('-o','--out-dir',default=None,type='string',metavar='DIR',help='Directory for output: defaults to gw150914/')
    parser.add_option('-t','--threads',default=None,type='int',metavar='N',help='Number of threads (default = 1/core)')
    parser.add_option('-f','--full-run',default=0,type='int',metavar='full_run',help='perform a full PE run')
    parser.add_option('--inject',default=False,action='store_true',help='Inject NR Signal')
    parser.add_option('--template',default='LAL',type='str',metavar='template',help='template to use for the run')
    parser.add_option('--seglen',default=4,type='float',metavar='seglen',help='length of the data stretch to analyse')
    parser.add_option('--flow',default=20,type='float',metavar='flow',help='low frequency cutoff')
    parser.add_option('--fhigh',default=500,type='float',metavar='fhigh',help='high frequency cutoff')
    parser.add_option('--nlive',default=1024,type='int',metavar='n',help='Live points')
    parser.add_option('--maxmcmc',default=1024,type='int',metavar='m',help='max MCMC points')
    (opts,args)=parser.parse_args()

    if opts.out_dir is None:
        opts.out_dir='./gw150914/'
    noise_model  = NoiseModel(chunk_size=opts.seglen,
                              flow=opts.flow,
                              fhigh=opts.fhigh)
    logZnoise=noise_model.log_likelihood(noise_model.new_point())
    print('Noise evidence {0}'.format(logZnoise))
    if opts.full_run:
        signal_model = GravitationalWaveModel(chunk_size=opts.seglen,
                                              template = opts.template,
                                              flow=opts.flow,
                                              fhigh=opts.fhigh)
        work=cpnest.CPNest(signal_model,
                           verbose=3,
                           Poolsize=1024,
                           Nthreads=opts.threads,
                           Nlive=opts.nlive,
                           maxmcmc=opts.maxmcmc,
                           output=opts.out_dir)
        work.run()
        print('Signal evidence {0}'.format(work.NS.logZ))
        logB = work.NS.logZ-logZnoise
        print('log B {0}'.format(logB))
        x = work.posterior_samples.ravel()
    else:
        x = np.genfromtxt(os.path.join(opts.out_dir,'posterior.dat'),names=True)

        logZ = np.loadtxt(os.path.join(opts.out_dir,'chain_{0}_1234.txt_evidence.txt'.format(opts.nlive)))[0]
        print('Signal evidence {0}'.format(logZ))
        logB = logZ-logZnoise
        print('log B {0}'.format(logB))
