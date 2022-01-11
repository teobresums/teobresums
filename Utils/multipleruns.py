#!/usr/bin/python3

"""
Script to auto generate parfiles and run them, if required

Given a template parfile with basic setup, 
generates parfiles for all combination of a given subset of parameters
(e.g. to vary binary masses and spins)

SB 10/2018
RG 12/2021
"""

import optparse as op; import eobutils as utils; import matplotlib.pyplot as plt; import numpy as np
import os 

ELLMAX = 8

parser=op.OptionParser()
parser.add_option('--fname', dest='basep', type='str', help='Name of the baseline parfile')
parser.add_option('--path',  dest='based', type='str', help='Name of the path to the parfile')
parser.add_option('--gen_pars', action='store_true', default=True, help='Generate parfiles from a baseline')
parser.add_option('--run', action='store_true', default=False, help='Run and remove the parfiles')
parser.add_option('--multiprocess', action='store_true', default=False, help='Run with multiprocess')
parser.add_option('--nproc', dest='nproc', type='int', help='Number of processes to run with multiprocess')
parser.add_option('--plot_wf', action='store_true', default=False, help='Plot the polarization')
parser.add_option('--plot_hlm', action='store_true', default=False, help='Plot the multipoles')

(opts,args) = parser.parse_args()

# Example usage:
# python3 multipleruns.py --fname GW150914.par --path ./ --gen_pars --run --multiprocess --nproc 2 --plot_wf

if __name__ == '__main__':

    # Base dir & parfile
    based = opts.based
    basep = opts.basep

    # Set new values/ranges for parameters (Use lists)
    q    = [1., 1.2, 1.4]
    chi1 = [0.]
    chi2 = [0.]

    # Pack them into a dictionary
    # NOTE: keys must match those in parfile otherwise ignored
    n = {'q': q,
         'chi1': chi1,
         'chi2': chi2
        }

    ###############################################################
    # Do not modify below here unless you know what you are doing!
    ###############################################################

    # Generate the parfiles
    if opts.gen_pars == True:
        parlist = utils.generate_parfiles(n,based,basep)

    # Run
    if opts.run == True:
        if opts.multiprocess == True:
            utils.run_batch_parfiles_MP(parlist, opts.nproc)
        else:
            utils.run_batch_parfiles(parlist)

    # Visualize
    for d in parlist:
        # plot the waveform
        d = d.strip(".par")

        if (opts.plot_wf):
            t, hp, hc, _, _ = np.loadtxt(based+d+'/waveform.txt', unpack=True)
            A = np.sqrt(hp**2+hc**2)
            plt.plot(t, A)
            plt.plot(t, hp)
            plt.plot(t, hc)
            plt.show()

        # plot the multipoles
        if(opts.plot_hlm):

            fig, axs = plt.subplots(3, 1)

            for ell in range(2, ELLMAX+1):
                for emm in range(1, ell+1):
                    
                    fn1 = "/hlm_l"+str(ell)+"_m"+str(emm)+".txt"
                    fn2 = "/hlm_interp_l"+str(ell)+"_m"+str(emm)+".txt"
                    full_path = None

                    if os.path.exists(based+d+fn1):
                        full_path = based+d+fn1
                    elif os.path.exists(based+d+fn2):
                        full_path = based+d+fn2

                    if isinstance(full_path, str):
                        t, A, p = np.loadtxt(full_path, unpack=True)
                        omg = np.zeros_like(p)
                        omg[1:] = np.diff(p)/np.diff(t)

                        axs[0].plot(t,abs(A),label=r'$'+str(ell)+str(emm)+'$')
                        axs[1].plot(t,p,  label=r'$\phi_{'+str(ell)+','+str(emm)+'}$')
                        axs[2].plot(t,omg,label=r'$\omega_{'+str(ell)+','+str(emm)+'}$')

            axs[0].legend()
            axs[0].set_ylabel('$ A $')
            axs[1].set_ylabel('$ \phi $')
            axs[2].set_ylabel('$ \omega $')
            axs[2].set_xlabel('$ t $')
            axs[1].text(0, 0.75*p[-1], d, fontsize = 12)

            plt.show()
