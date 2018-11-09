from pyTEOBResumS import pyTEOBResumS
import numpy as np
#import lalsimulation as lalsim

m1=61.27
m2=59.77
spin1x = 0.0
spin1y = 0.0
spin1z = -0.58
spin2x = 0.0
spin2y = 0.0
spin2z = -0.57
inclination = 1.6604149
phase = 0.6868080
f_min = 10.0
sampling_rate = 2048.
dt = 1./sampling_rate

LambdaAl2 = 0.0
LambdaBl2 = 0.0
LambdaAl3 = 0.0
LambdaBl3 = 0.0
LambdaAl4 = 0.0
LambdaBl4 = 0.0
distance = 500.0

flags ={'NQC':1,
        'tidal':0,
        'speedy':1,
        'dynamics':0,
        'solver_scheme':0,
        'RWZ':0,
        'Yagi_fits':0,
        'spin':1,
        'multipoles':-1,
        'geometric_units':0,
        'set':0
        }

lm = -1

wave_flags = None
non_GR_params = None
amp_order=0
phase_order=7
#approx = lalsim.SEOBNRv4

import matplotlib.pyplot as plt
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

for f_min in [5,10,20]:#,np.pi/3.,np.pi/2.]:
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
                    f_min,
                    dt,
                    LambdaAl2,
                    LambdaBl2,
                    LambdaAl3,
                    LambdaBl3,
                    LambdaAl4,
                    LambdaBl4,
                    distance,
                    lm,
                    flags)
                    
#    hp,hc = lalsim.SimInspiralChooseTDWaveform(phase,
#                                 dt,
#                                 m1*lalsim.lal.MSUN_SI, m2*lalsim.lal.MSUN_SI,
#                                 spin1x, spin1y, spin1z,
#                                 spin2x, spin2y, spin2z,
#                                 f_min, f_min,
#                                 distance*lalsim.lal.PC_SI*1e6,
#                                 inclination,
#                                 LambdaAl2, LambdaBl2,
#                                 wave_flags, non_GR_params,
#                                 amp_order, phase_order,
#                                 approx)

    t = dt*np.arange(h.shape[0])
    ax1.plot(t, h[:,0],label=str(f_min))
#    ax1.plot(hp.data.data,linestyle='dashed')
    ax2.plot(t, h[:,1])
#    ax2.plot(hc.data.data,linestyle='dashed')
ax1.legend()
plt.show()
