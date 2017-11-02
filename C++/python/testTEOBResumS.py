from pyTEOBResumS import pyTEOBResumS
import numpy as np

m1=40.0
m2=40.0
spin1x = 0.0
spin1y = 0.0
spin1z = 0.1
spin2x = 0.0
spin2y = 0.0
spin2z = 0.2
inclination = 0.0
polarisation = 0.0
f_min = 20.0
sampling_rate = 2048.
dt = 1./sampling_rate
LambdaAl2 = 0.0
LambdaBl2 = 0.0
LambdaAl3 = 0.0
LambdaBl3 = 0.0
LambdaAl4 = 0.0
LambdaBl4 = 0.0
distance = 400.0

flags ={'NQC':'0',
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
                LambdaBl2,
                 LambdaAl3,
                 LambdaBl3,
                 LambdaAl4,
                 LambdaBl4,
                 distance,
                lm,
                flags)

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(h[:,0])
ax.plot(h[:,1])

hf = np.fft.rfft(h[:,0])
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(hf.real)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(h[:,0]**2+h[:,1]**2)
plt.show()
