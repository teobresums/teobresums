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

    def __init__(self, inject=False, chunk_size=4.0, trigtime=1126259462.423, **kwargs):
        
        super(GravitationalWaveModel,self).__init__(**kwargs)
        # this is the merger time in H1
        self.tevent = trigtime

        self.inject=inject

        self.detectors = [GravitationalWaveDetector('H1','data/H-H1_LOSC_4_V1-1126259446-32.txt', chunk_size=chunk_size, trigtime=trigtime, **kwargs),GravitationalWaveDetector('L1','data/L-L1_LOSC_4_V1-1126259446-32.txt', chunk_size=chunk_size, trigtime=trigtime, **kwargs)]

        self.sampling_rate = self.detectors[0].sampling_rate
        self.segment_length = self.detectors[0].segment_length
        self.dt = 1./self.sampling_rate
        self.Flow = self.detectors[0].Flow
        
        #parameters
        self.names=['ra',
           'dec',
           'tc',
           'mc',
           'q',
           'spin1z',
           'spin2z',
           'iota',
           'psi',
           'logdistance']

        self.bounds=[[0,2.0*np.pi],
                [-np.pi/2.0,np.pi/2.0],
                [self.tevent-0.05,self.tevent+0.05],
                [20,40],
                [0.3,1.0],
                [-0.3,0.3],
                [-0.3,0.3],
                [0.0,np.pi],
                [0.0,np.pi],
                [np.log(1),np.log(2000)]]
                
        self.flags ={'NQC':'1',
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
        self.window=tukey(self.segment_length,0.5)
        self.windowNorm = self.segment_length/np.sum(self.window**2)
            
    def log_likelihood(self,x):
        
        mc = x['mc']
        q = x['q']
        d = np.exp(x['logdistance'])
        m1, m2 = McQ2Masses(mc, q)

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
                         self.Flow,
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
        
        hptilde = np.fft.rfft(hp)*self.windowNorm
        hctilde = np.fft.rfft(hc)*self.windowNorm

        return np.sum([d.logLikelihood(hptilde, hctilde, x['ra'], x['dec'], x['psi'], x['tc']) for d in self.detectors])
    
    def log_prior(self, x):
        if np.isfinite(super(GravitationalWaveModel,self).log_prior(x)):
            logP = 2.0*x['logdistance']
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

    def __init__(self, chunk_size=4.0, trigtime=1126259462.43, **kwargs):
        
        super(NoiseModel,self).__init__(**kwargs)
        # this is the merger time in H1
        self.tevent = 1126259462.423

        self.detectors = [GravitationalWaveDetector('H1','data/H-H1_LOSC_4_V1-1126259446-32.txt', chunk_size=chunk_size, trigtime=trigtime,**kwargs),
                          GravitationalWaveDetector('L1','data/L-L1_LOSC_4_V1-1126259446-32.txt', chunk_size=chunk_size, trigtime=trigtime,**kwargs)]

    def log_likelihood(self,x):
        
        Nfd = self.detectors[0].segment_length/2+1
        h = np.zeros(Nfd)
        
        return np.sum([d.logLikelihood(h, h, 0.0, 0.0, 0.0, 0.0) for d in self.detectors])

if __name__=='__main__':
    parser=OptionParser()
    parser.add_option('-o','--out-dir',default=None,type='string',metavar='DIR',help='Directory for output: defaults to gw150914/')
    parser.add_option('-t','--threads',default=None,type='int',metavar='N',help='Number of threads (default = 1/core)')
    parser.add_option('-f','--full-run',default=1,type='int',metavar='full_run',help='perform a full PE run')
    parser.add_option('--inject',default=False,action='store_true',help='Inject NR Signal')
    parser.add_option('--nlive',default=1024,type='int',metavar='n',help='Live points')
    parser.add_option('--maxmcmc',default=1024,type='int',metavar='m',help='max MCMC points')
    (opts,args)=parser.parse_args()

    if opts.out_dir is None:
        opts.out_dir='./gw150914/'
    noise_model  = NoiseModel()
    logZnoise=noise_model.log_likelihood(noise_model.new_point())
    print('Noise evidence {0}'.format(logZnoise))
    if opts.full_run:
        signal_model = GravitationalWaveModel()
        work=cpnest.CPNest(signal_model,
                           verbose=3,
                           Poolsize=128,
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
#        print('Signal evidence {0}'.format(logZ))
#        logB = logZ-logZnoise
#        print('log B {0}'.format(logB))
