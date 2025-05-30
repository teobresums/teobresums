import numpy as np
from numpy.lib import recfunctions as rfn
import matplotlib.pyplot as plt
import matplotlib
import multiprocessing as mp
from parspace_and_plots import write_dict_to_txt, gen_random_pars
from utilities import CreateDict
import EOBRun_module as EOB
import time
import argparse
import tqdm

try:
    from PyART.models import teob
    from PyART import waveform as wf
    from PyART.analysis import match
    from PyART.utils import wf_utils as wfu
    pyart = True
except:
    print("Warning! PyART not installed; not computing mismatch or phase difference.")
    pyart = False

MSUN = 4.925490947641267e-06

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.size']   = 15

match_settings = {'kind'                : 'single-mode',
                  'final_frequency_mm'  : 2048.,
                  'dt'                  : 1./8192.,
                  'iota'                : 0.,
                  'coa_phase'           : 0.,
                  'eff_pols'            : 0.,
                  'debug'               : False,
                  'modes-or-pol'        : 'modes',
                  'modes'               : [(2,2)],
                  'resize_factor'       : 4,
                  'taper_start'         : 0.1,
                  'psd'                 : 'flat',
                  'asd_file'            : None,
                  'cut_longer'          : False,
                  'cut_second_waveform' : True,
                  'taper_alpha'         : 0.1,
                  'pad_end_frac'        : 0.5,
                  'pre_align_shift'     : 0.,
                  'pre_align'           : True,
                  'taper'               : 'sigmoid',
                }

def run_tols(delta_log_tol, pardic, f0):
    """
    Generate EOB waveform with default and modified tolerances.
    Compute mismatch, walltime, phase difference.
    Return everything essentially.
    """
    newdic = pardic.copy()
    newdic = CreateDict(M=pardic['M'],
                        q=pardic['q'],
                        f0=f0,
                        chi1=[pardic['chi1x'], pardic['chi1y'], pardic['chi1z']],
                        chi2=[pardic['chi2x'], pardic['chi2y'], pardic['chi2z']],
                        use_geom="no",
                        ecc=pardic['ecc'],
                        anomaly=pardic['anomaly'],)

    try:
        clock0 = time.time()
        t0, _, _, hlm0, _ = EOB.EOBRunPy(newdic)
        clock0 = time.time() - clock0
        t0 = t0 - t0[0]

        # Build waveform object by hand for match
        if pyart:
            eob0 = wf.Waveform()
            eob0._u      = t0
            eob0._hlm    = teob.convert_hlm(hlm0)
            eob0._kind   = 'EOB'
            eob0.domain  = 'Time'
            # Alignment time interval
            t_align = np.linspace(eob0.u[0] + 0.1*(eob0.u[-1] - eob0.u[0]), eob0.u[0] + 0.75*(eob0.u[-1] - eob0.u[0]), 10000)
        
        change_tols = [{'ode_abstol': 10**(-13 + delta_log_tol), 'ode_reltol': 10**(-11)}, 
                       {'ode_abstol': 10**(-13),                 'ode_reltol': 10**(-11 + delta_log_tol)}, 
                       {'ode_abstol': 10**(-13 + delta_log_tol), 'ode_reltol': 10**(-11 + delta_log_tol)}]
        times   = []
        dphis   = []
        matches = []
        for tol_dict in change_tols:
            newdic = {**newdic, **tol_dict}
            clock = time.time()
            t, _, _, hlm, _ = EOB.EOBRunPy(newdic)
            times.append((time.time() - clock)/clock0)
            t = t - t[0]

            # Build waveform object by hand for match
            if pyart:
                eob = wf.Waveform()
                eob._u     = t
                eob._hlm   = teob.convert_hlm(hlm)
                eob._kind  = 'EOB'
                eob.domain = 'Time'

                # Align with base case
                tau, dphi, _ = wfu.Align(t_align, t_align[-1], 0.1*t_align[-1], eob.u, eob.hlm[(2,2)]['p'], eob0.u, eob0.hlm[(2,2)]['p'])
                eob._u = eob.u + tau
                eob._hlm[(2,2)]['p'] = eob.hlm[(2,2)]['p'] - dphi
                eob._hlm[(2,2)]['z'] = eob.hlm[(2,2)]['z'] * np.exp(-1j*dphi)
                eob._hlm[(2,2)]['real'] = eob.hlm[(2,2)]['z'].real
                eob._hlm[(2,2)]['imag'] = eob.hlm[(2,2)]['z'].imag

                # Get phase difference
                dphis.append((np.interp(eob0.u, eob.u, eob.hlm[(2,2)]['p']) - eob0.hlm[(2,2)]['p'])[np.argmax(eob0.hlm[(2,2)]['A'])])

                # Compute mismatch
                mismatch = match.Matcher(eob0, eob, settings={**match_settings, **{'M': newdic['M'], 'initial_frequency_mm': 1.35*newdic['initial_frequency']}})
                matches.append(mismatch.mismatch)
            else:
                dphis   = [1, 1, 1]
                matches = [1, 1, 1]
    except:
        times   = [1, 1, 1]
        dphis   = [1, 1, 1]
        matches = [1, 1, 1]

    return {**pardic, **{'delta_log_tol': delta_log_tol, 'time_abs': times[0],  'dphi_abs': dphis[0],  'match_abs': matches[0],
                                                         'time_rel': times[1],  'dphi_rel': dphis[1],  'match_rel': matches[1],
                                                         'time_both': times[2], 'dphi_both': dphis[2], 'match_both': matches[2]}}

def dict_to_vec(dic):
    """
    Transpose a dictionary of vectors into a vector of dictionaries
    """

    out = []
    for j in range(len(dic[list(dic.keys())[0]])):
        data = {key: dic[key][j] for key in dic}
        out.append(data)
    
    return out

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Explore the effect of different tolerances on the EOB waveform generation.')
    parser.add_argument('--N',          type=int,            default=10000,          help='Number of random parameter sets to generate.')
    parser.add_argument('--M',          type=float, nargs=2, default=[10., 200.],    help='Mass range (in solar masses) for the binary.')
    parser.add_argument('--q',          type=float, nargs=2, default=[1., 10.],      help='Mass ratio range for the binary.')
    parser.add_argument('--f0',         type=float,          default=10.,            help='Initial frequency.')
    parser.add_argument('--chi',        type=float, nargs=2, default=[-0.99, 0.99],  help='Spin range for the binary.')
    parser.add_argument('--eccentric',  action='store_true',                         help='Eccentric?')
    parser.add_argument('--ecc',        type=float, nargs=2, default=[0., 0.6],      help='Eccentricity range for the binary.')
    parser.add_argument('--anomaly',    type=float, nargs=2, default=[0., 2.*np.pi], help='Anomaly range for the binary.')
    parser.add_argument('--precessing', action='store_true',                         help='Precessing?')
    parser.add_argument('--dlogtol',    type=int,            default=4,              help='Log tolerance increase w.r.t. default.')
    parser.add_argument('--nprocs',     type=int,            default=8,              help='Number of processes if running parallel.')
    parser.add_argument('--plot',       action='store_true',                         help='Plot results?')
    parser.add_argument('--file',       type=str,            default=None,           help='Input file')
    parser.add_argument('--px',         type=str,            default='q',            help='x variable on plots')
    parser.add_argument('--py',         type=str,            default='chieff',       help='y variable on plots')
    parser.add_argument('--pc',         type=str,            default='match',        help='Color variable on plots')
    parser.add_argument('--worst',      action='store_true',                         help='Plot the worst matching wfs according to the given pc variable')

    args = parser.parse_args()

    if args.eccentric:
        ecc_int     = args.ecc
        anomaly_int = args.anomaly
    else:
        ecc_int     = None
        anomaly_int = None
    outfile = f"tolerance_N{args.N}_f0{args.f0}_dlogtol{args.dlogtol}_M{args.M[1]}_q{args.q[1]}{f'_ecc{ecc_int[1]}' if ecc_int is not None else ''}{'_prec' if args.precessing else ''}.txt"
    
    if args.file is None:
        pars = gen_random_pars(Mint=args.M, qint=args.q, chi_int=args.chi, ecc_int=ecc_int, anomaly_int=anomaly_int, precessing=args.precessing, N=args.N)
        pars = dict_to_vec(pars)
    
        # Create output file and write header
        dummy = {**pars[0], **{'delta_log_tol': 0, 'time_abs':  0,  'dphi_abs': 0,  'match_abs': 0,
                                                   'time_rel':  0,  'dphi_rel': 0,  'match_rel': 0,
                                                   'time_both': 0, 'dphi_both': 0, 'match_both': 0}}
        write_dict_to_txt(outfile, dummy, True)
    
        with mp.Pool(args.nprocs) as p:
            # for jj in range(args.N):
            #     thing = p.apply_async(run_tols, (args.dlogtol, pars[jj], args.f0))
            #     outdict = thing.get()
            #     write_dict_to_txt(outfile, outdict, False)
            for jj in tqdm.tqdm([[ii, p.apply_async(run_tols, (args.dlogtol, pars[ii], args.f0))] for ii in range(args.N)]):
                outdict = jj[1].get()
                write_dict_to_txt(outfile, outdict, False)
    else:
        outfile = args.file
    
    data = np.genfromtxt(outfile, names=True)
    temp = (data['q']*data['chi1z'] + data['chi2z'])/(1. + data['q'])
    data = rfn.append_fields(data, 'chieff', temp)
    temp = np.maximum((np.sqrt(data['chi1x']**2 + data['chi1y']**2)), (4. + 3.*data['q'])/((4.*data['q'] + 3.)*data['q'])*(np.sqrt(data['chi2x']**2 + data['chi2y']**2)))
    data = rfn.append_fields(data, 'chip',   temp)
    temp = args.f0*data['M']*MSUN
    data = rfn.append_fields(data, 'f0',     temp)

    varlabels = {'match_abs' : r'$\mathcal{F}_{\rm abs}$',
                 'match_rel' : r'$\mathcal{F}_{\rm rel}$',
                 'match_both': r'$\mathcal{F}_{\rm both}$',
                 'time_abs'  : r'$\Delta t_{\rm abs}/\Delta t_0$',
                 'time_rel'  : r'$\Delta t_{\rm rel}/\Delta t_0$',
                 'time_both' : r'$\Delta t_{\rm both}/\Delta t_0$',
                 'dphi_abs'  : r'$\delta \phi_{22}^{\rm abs}$',
                 'dphi_rel'  : r'$\delta \phi_{22}^{\rm rel}$',
                 'dphi_both' : r'$\delta \phi_{22}^{\rm both}$',
                 'q'         : r'$q$',
                 'chieff'    : r'$\chi_{\rm eff}$',
                 'chip'      : r'$\chi_{\rm p}$',
                 'ecc'       : r'$e_0$',
                 'anomaly'   : r'$\zeta_0$',
                 'M'         : r'$M$',
                 'f0'        : r'$M f_0$'}

    if args.plot:
        for tolt in ['abs', 'rel', 'both']:
            fig, ax = plt.subplots(layout='constrained', figsize=(8, 6))
            var = f'{args.pc}_{tolt}'
            data[var] = abs(data[var])
            norm_var = matplotlib.colors.LogNorm(vmin=max(min(data[var]), 1.e-12), vmax=max(data[var]) if 'time' not in var else 1.)
            tp = ax.scatter(data[args.px], data[args.py], c=data[var], norm=norm_var, s=10)
            cbar = plt.colorbar(tp)
            ax.set_xlabel(varlabels[args.px])
            ax.set_ylabel(varlabels[args.py])
            cbar.ax.set_yscale('log')
            cbar.ax.set_ylabel(varlabels[var])
        
        if args.worst:
            jw = np.argmax(data[f'{args.pc}_both'])
            print("Worst matching waveform:")
            print(f"M          = {data['M'][jw]}")
            print(f"q          = {data['q'][jw]}")
            print(f"chi1       = ({data['chi1x'][jw]}, {data['chi1y'][jw]}, {data['chi1z'][jw]})")
            print(f"chi2       = ({data['chi2x'][jw]}, {data['chi2y'][jw]}, {data['chi2z'][jw]})")
            print(f"ecc        = {data['ecc'][jw]}")
            print(f"anomaly    = {data['anomaly'][jw]}")
            print(f"time_both  = {data[f'time_both'][jw]:.6e}")
            print(f"dphi_both  = {data[f'dphi_both'][jw]:.6e}")
            print(f"match_both = {data[f'match_both'][jw]:.6e}")

            pardic = CreateDict(M=data['M'][jw],
                                q=data['q'][jw],
                                chi1=[data['chi1x'][jw], data['chi1y'][jw], data['chi1z'][jw]],
                                chi2=[data['chi2x'][jw], data['chi2y'][jw], data['chi2z'][jw]],
                                f0=args.f0,
                                use_geom="no",
                                ecc=data['ecc'][jw],
                                anomaly=data['anomaly'][jw])
            clock = time.time()
            t0, hp0, hc0, hlm0, dyn0 = EOB.EOBRunPy(pardic)
            t0 = t0 - t0[0]
            pardic.update({'ode_abstol': 10**(-13 + args.dlogtol),
                           'ode_reltol': 10**(-11 + args.dlogtol)})
            t,  hp,  hc,  hlm,  dyn  = EOB.EOBRunPy(pardic)
            t = t - t[0]

            if pyart:
                t_align = np.linspace(t0[0] + 0.1*(t0[-1] - t0[0]), t0[0] + 0.75*(t0[-1] - t0[0]), 10000)
                tau, dphi, _ = wfu.Align(t_align, t_align[-1], 0.1*t_align[-1], t, hlm['1'][1], t0, hlm0['1'][1])
                t = t + tau
                deltaphi = np.interp(t0, t, hlm['1'][1]) - dphi - hlm0['1'][1]

                fig, ax = plt.subplots(2, 1, layout='constrained', figsize=(9,8), sharex=True)
                ax[0].plot(t0, hlm0['1'][0]*np.cos(hlm0['1'][1]), color='k')
                ax[0].plot(t,  hlm['1'][0]*np.cos(hlm['1'][1] - dphi))
                ax[1].plot(t0, deltaphi)
                ax[0].set_ylabel(r'$\Re h_{22}$')
                ax[1].set_xlabel(r'$t/M$')
                ax[1].set_ylabel(r'$\delta \phi_{22}$')
            else:
                fig, ax = plt.subplots(layout='constrained', figsize=(9,6), sharex=True)
                ax.plot(t0, hlm0['1'][0]*np.cos(hlm0['1'][1]), color='k')
                ax.plot(t,  hlm['1'][0]*np.cos(hlm['1'][1]))
                ax.set_ylabel(r'$\Re h_{22}$')
                ax.set_xlabel(r'$t/M$')
        plt.show()