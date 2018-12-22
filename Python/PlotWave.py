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
    Build a file name from fname
    """
    x = os.path.abspath(fname).split("/")
    return x[-2]+"_"+os.path.splitext(x[-1])[0]    
    
if __name__ == "__main__": 

    #
    # Parse input
    #
    parser = ArgumentParser(description="Simple plot of waveform files from TEOBResumS")
    parser.add_argument('-i', '--input-files', dest="filenames", nargs='+', required=True,
                        help="Input files", metavar="FILE")
    parser.add_argument('-o-', '--output-dir', dest="outputdir", nargs='?', default=os.getcwd(),
                        help="Output directory")    
    parser.add_argument('-m', '--mode', dest='mode', nargs=1, default="aor", 
                        help="Type of plot. A string containing one or more of the chars 'apori' for amplitude, phase, frequency, real, and imaginary")
    parser.add_argument('-x', '--xrange', dest="xrange", type=float, nargs='+', default=[],
                        help="Range of xaxis") # Note argparse might have problems with scientific notation.   
    parser.add_argument('-y', '--yrange', dest="yrange", type=float, nargs='+', default=[],
                        help="Range of yaxis")    
    parser.add_argument('--no-title', dest='showtitle', action='store_false', default=True,
                        help="Do not show title")
    parser.add_argument('--no-legend', dest='showlegend', action='store_false', default=True,
                        help="Do not show legend")
    parser.add_argument('--no-show', dest='showplot', action='store_false', default=True,
                        help="Do not show plot")
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    #
    # Checks
    #
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

        #
        # Load data
        #
        t = id_data_file(f)
        if t is "triap":
            # t:0 real:1 imag:2 Ampli:3 Phase:4
            t, Reh, Imh, A, phi = np.loadtxt(f, unpack=True)
            if not A.any():
                A = np.sqrt(Reh**2 +Imh**2)
            if not phi.any():
                phi = np.unwrap(-np.angle(Reh+1j*Imh))
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

        #
        # Compute Frequency
        #
        omg = np.diff(phi)/np.diff(t)

        #
        # Plots
        #
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

        #plt.xlabel('$t/M$')            
        plt.xlabel('time')
        plt.ylabel('h')
        #plt.grid(True)
        if args.xrange: plt.xlim(args.xrange)
        if args.yrange: plt.xlim(args.yrange)
        if args.showlegend: plt.legend(loc='upper left')
            
        #
        # Build a name for output and title
        #
        name = build_output_fname(f)
        
        if args.showtitle: plt.title(name)        

        #
        # Show/save
        #
        if args.showplot: plt.show()
        plt.savefig(args.outputdir+"/"+name+".png", dpi=400)

                
