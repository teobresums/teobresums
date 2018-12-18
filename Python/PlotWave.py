#!/usr/bin/python3

"""
Script to plot TEOBResumS waveform

Example

$ python Run.py sim1/waveform.txt sim2/waveform.txt ...
$ python Run.py sim?/waveform.txt 

SB 12/2018
"""

import sys, glob
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__": 
    
    args = sys.argv[1:]
    if not len(args):
        print("Usage: python "+sys.argv[0]+" <files>")
    else:
        for a in args:
            for fname in glob.glob(a):
              
              # t:0 real:1 imag:2 Ampli:3 Phase:4
              #t, Reh, Imh, A, phi = np.loadtxt(fname, unpack=True)
              t, Reh, Imh = np.loadtxt(fname, unpack=True)

              A = np.sqrt(Reh**2 +Imh**2)
              phi = np.unwrap(-np.angle(Reh+1j*Imh))
              omg = np.diff(phi)/np.diff(t)

              plt.plot(t, A, label=r"$A/M$")
              plt.plot(t, Reh, label=r"$\Re{h}$")
              plt.plot(t[1:], omg, label=r"$M\omega$")
              plt.xlabel('t/M')
              #plt.ylabel('...')
              plt.title(fname)
              #plt.grid(True)
              #plt.xlim(0, t[-1])
              plt.ylim(-.75, .75)
              plt.legend(loc='upper left')
              plt.savefig(".png")
              plt.show()

        
