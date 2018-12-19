#!/usr/bin/python3

"""
Simple script to plot TEOBResumS waveform
SB 12/2018
"""

import sys, os
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import numpy as np
from ast import literal_eval

from EOBUtils import * # id_data_file

def build_output_fname(fname):
    """
    Build a file name from the 
    """
    x = os.path.abspath(fname).split("/")
    return x[-2]+"_"+os.path.splitext(x[-1])[0]    
    
if __name__ == "__main__": 

    # Parse input
    parser = ArgumentParser(description="Simple plot of waveform files from TEOBResumS")
    parser.add_argument("-i", dest="filenames", 
                        nargs='+', required=True,
                        help="Input files", metavar="FILE")
    parser.add_argument("-o", dest="outputdir", 
                        nargs='?', default=os.getcwd(),
                        help="Output directory")    
    parser.add_argument("-m", dest="mode", 
                        nargs=1, default="aor", 
                        help="Type of plot. A string containing one or more of the chars 'apori' for amplitude, phase, frequency, real, and imaginary")
    parser.add_argument("-s", dest="show", 
                        nargs=1, default="yes",
                        help="Show plot")
    parser.add_argument("-a", dest="yrange",
                        type=float, nargs='+',
                        default=[-0.75,0.75], 
                        help="Ranges of yaxis")    

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    
    # Info/checks
    args = parser.parse_args()
    files_found = []
    for f in args.filenames:
        if os.path.isfile(f):
            files_found.append(f)
            #print(f)
    if not files_found:
        print("No file was found")
        sys.exit()        

    if not os.path.exists(args.outputdir):
        os.mkdir( args.outputdir, 0O755 )
    print("Process %d files of %d" % (len(files_found),len(args.filenames)))

    mode = str(args.mode) # "apori"
    
    for f in files_found:

        # Load data
        t = id_data_file(f)
        if t is "triap":
            # t:0 real:1 imag:2 Ampli:3 Phase:4
            t, Reh, Imh, A, phi = np.loadtxt(f, unpack=True)
        if t is "tri":
            # t:0 real:1 imag:2
            t, Reh, Imh = np.loadtxt(f, unpack=True)
            A = np.sqrt(Reh**2 +Imh**2)
            phi = np.unwrap(-np.angle(Reh+1j*Imh))
        if t is "tap":
            # t:0 Ampli:1 Phase:2
            t, A, phi = np.loadtxt(f, unpack=True)
            h = A* np.exp(-1j*phi)
            Reh = np.real(h)
            Imh = np.imag(h)
        if t is None:
            continue
        
        # Compute Frequency
        omg = np.diff(phi)/np.diff(t)

        # Plots
        if "a" in mode:
            plt.plot(t, A, color='blue', label=r"$A/(M\nu)$")
        if "p" in mode:
            plt.plot(t, phi, color='gray', label=r"$\phi$")
        if "o" in mode:
            plt.plot(t[1:], omg, color='red', label=r"$M\omega$")
        if "r" in mode:
            plt.plot(t, Reh, color='cyan', label=r"$\Re{h/(M\nu)}$", alpha=0.6)
        if "i" in mode:
            plt.plot(t, Imh, color='orange', label=r"$\Im{h/(M\nu)}$", alpha=0.6)

        name = build_output_fname(f)
        
        plt.xlabel('$t/M$')
        plt.title(name)
        #plt.grid(True)
        #plt.xlim(0, t[-1])
        plt.ylim(args.yrange)
        plt.legend(loc='upper left')

        plt.savefig(args.outputdir+"/"+name+".png", dpi=400)

        if str(args.show) is "yes": plt.show()

        
