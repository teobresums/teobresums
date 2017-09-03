from pyTEOBResumS import pyTEOBResumS
import numpy as np

m1=40.0
m2=40.0
spin1x = 0.0
spin1y = 0.0
spin1z = 0.2
spin2x = 0.0
spin2y = 0.0
spin2z = 0.2
inclination = 0.0
polarisation = 0.0
f_min = 20.0
sampling_rate = 4096.
LambdaAl2 = 0.0
LambdaBl2 = 0.0
distance = 40.0
NQC = 0
tidal = 0
speedy = 1
RWZ = 0
dynamics = 0
waveform = 0
lm = 1
solver_scheme = 0

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
                sampling_rate,
                LambdaAl2,
                LambdaBl2,
                distance,
                NQC,
                tidal,
                speedy,
                RWZ,
                dynamics,
                waveform,
                lm,
                solver_scheme)

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
