import os; import re; import json; import collections; import argparse

VERBOSE = 1

# name : [description, type, default value]
PARKEYS = {
    'M': ['Total mass of the system', 'float', 1. ],
    'q': ['Mass ratio of the system','float',  1. ],
    'LambdaAl2': ['Quadrupolar tidal parameter of body 1','float', 0. ],
    'LambdaBl2': ['Quadrupolar tidal parameter of body 2', 'float', 0. ],
    'chi1': ['Spin of body 1', 'float', 0. ],
    'chi2': ['Spin of body 2', 'float', 0. ],
    'chi1x': ['x component of the spin of body 1', 'float', 0. ],
    'chi1y': ['y component of the spin of body 1', 'float', 0. ],
    'chi1z': ['z component of the spin of body 1', 'float', 0. ],
    'chi2x': ['x component of the spin of body 2', 'float', 0. ],
    'chi2y': ['y component of the spin of body 2', 'float', 0. ],
    'chi2z': ['z component of the spin of body 2', 'float', 0. ],
    #
    'r0': ['Initial separation in mass-rescaled units', 'float', 20. ],
    'initial_frequency' : ['Initial frequency in Hz or mass-rescaled units, depending on use_geometric_units', 'float', 0.004 ],
    'use_geometric_units': ['Use mass rescaled units or physical units', 'str',  "yes"],
    #
    'use_spins': ['0 = nonspinning (deprecated), 1 = spin-aligned, 2 =precessing spins', 'int', 1 ],
    'tides': ['Tidal model (NNLO, TEOBRESUM, TEOBRESUM3)', 'str', "TEOBRESUM" ],
    'tides_gravitomagnetic': ['Use gravitomagnetic tidal coefficients', 'str', "no" ],
    'use_a6c_fits': ['Use fits for a6c', 'str', "HM"],
    'use_cN3LO_fits' : ['Use fits for cN3LO', "HM"],
    'use_fmode': ['Use the fmode resonance model', 'str', "no" ],
    'use_lambda234_fits': ['Use quasi-universal relations to fit Lambda4,5,6,7,8 (YAGI13, GODZIEBA20)', 'str', "YAGI13" ],
    'use_speedytail': ['Speed up the computation of the tail factor in h_lm and F_lm', 'str', "no" ],
    'pGSF_tidal': ['GSF p-exponent (default=4)', 'float', 4. ],
    'dt_merger_interp': ['dt used for the interpolation of the merger dynamics, after the peak of Omega', 'float', 0.1 ],
    'dt_interp': ['dt used to interpolate the waveform', 'float', 0.1 ],
    'dt': ['dt used if solving the ODEs with a uniform timestep', 'float', 0.5 ],
    'srate_interp': ['Sampling rate, 1/dt at which to interpolate the waveform', 'float', 0.5 ],
    'srate': ['Sampling rate, 1/dt at which the waveform is computed if using a uniform timestep', 'float', 0.5 ],
    'size': ['Initial guess for the size of the waveform', 'int', 100 ],
    'ringdown_extend_array': ['Extend the ringdown array', 'int', 100 ],
    'distance': ['Distance of the source from Earth [Mpc,  ]', 'float', 1 ],
    'inclination': ['Angle between \hat{N} and the orbital angular momentum at the initial time', 'float', 0 ],
    'polarization': ['Polarization of the waves (unused)', 'float', 0 ],
    'coalescence_angle': ['Reference coalescence phase', 'float', 0 ],
    'use_mode_lm': ['List of modes to use', 'int list', "[1]" ],
    'output_lm': ['List of modes to output', 'int list', "[1]" ],
    'kpostpeak': ['List of modes for which NQCs are extracted from the pospeak template', 'int list', "[0,3,6,7,8]"] ,
    #
    'postadiabatic_dynamics': ['Turn on/off the postadiabatic approximation', 'str', "yes" ],
    'postadiabatic_dynamics_N': ['Use the N-th order PA (default=8)', 'int', 8 ],
    'postadiabatic_dynamics_size': ['Guess for the size of the postadiabatic dynamics', 'int', 500 ],
    'postadiabatic_dynamics_rmin': ['Stop the postadiabatic dynamics at this value of separation (default=14)', 'float', 20 ],
    'postadiabatic_dynamics_dr': ['Step dr of the postadiabatic dynamics', 'float', 0.1 ],
    'postadiabatic_dynamics_stop':['Stop after the postadiabatic dynamics', 'str', "no" ],
    #
    'centrifugal_radius': ['Model for the centrifugal radius', 'str', "CENTRAD_NLO" ],
    'use_flm': ['Model for the radiation reaction', 'str', "USEFLM_HM" ],
    'nqc': ['Choose how to set the NQCs (manual, auto, no)', 'str', "auto" ],
    'nqc_coefs_flx': ['Option for the NQC model used in the flux', 'str', "nrfit_spin202002"],
    'nqc_coefs_hlm': ['Option for the NQC model used in the waveform', 'str', "nrfit_spin202002"  ],
    'nqc_coefs_flx_file' : ['File containing the NQC coefficients for the flux', 'str',  None],
    'compute_LR': ['Compute the light ring (LR)', 'str', "no" ],
    'compute_LR_guess': ['Initial guess for the light ring (LR)', 'str', 4 ],
    'compute_LSO': ['Compute the last stable orbit (LSO)', 'str',  "no"],
    'compute_LSO_guess': ['Guess for the last stable orbit (LSO)', 'str', 4 ],
    'project_spins': ['Project the spins along the dynamics', 'str', "yes" ],
    'spin_interp_domain': ['Evaluate spins(time) or spins(frequency)', 'int',  1],
    'spin_flx': ['Model for the spin-evolution \dot\omega', 'str', "SPIN_FLX_PN" ],
    'ringdown_eulerangles': ['Option for the extension of the Euler angles beyond merger', 'str', "RD_EULERANGLES_QNMs" ],
    #
    'output_hpc': ['Output the plus and cross polarizations', 'str', "no" ],
    'output_multipoles': ['Output the waveform multipoles h_lm', 'int list', "no" ],
    'output_dynamics': ['Output the EOB dynamics', 'str', "no" ],
    'output_nqc': ['Output the NQCs', 'str', "no" ],
    'output_nqc_coefs': ['Output the NQC coefficients', 'str', "no" ],
    'output_ringdown': ['Output the ringdown portion of the waveform', 'str', "no" ],
    #
    'ode_timestep': ['Choose the ODE timestep between adaptive, uniform, adaptive+uniform after LSO', 'str', "adaptive" ],
    'ode_abstol': ['ODE absolute numerical tolerance', 'float', 1e-11 ],
    'ode_reltol': ['ODE relative numerical tolerance', 'float', 1e-13 ],
    'ode_tmax': ['Time stopping condition for the ODE solver', 'float', 1e7 ],
    'ode_stop_at_radius': ['Stop the ODE evolution at r', 'float', 2 ],
    'df': ['Frequency spacing for TEOBResumSPA', 'float', 1 ],
    'time_shift_FD': ['Time shift the FD waveform (merger at zero)', 'str', "yes" ],
    'interp_freqs': ['Use a user-input list of frequency to interpolate TEOBResumSPA', 'str', "no" ],
    'freqs': ['List of frequencies required by interp_freqs', 'float list', None ],
    'output_dir': ['directory where to perform the output', 'str', "data" ],
    "interp_uniform_grid": ['Interpolate on a uniform grid?', 'str', "yes" ],
    }

HEAD = """
## Parameters\n

Table automaticaly generated with  
`python3 parfile.py --docs 1`  

| Name      | Meaning                                                                             | Type    | Default | Reference |
|-----------|-------------------------------------------------------------------------------------|---------|---------|-----------|
"""

PARKEYS = collections.OrderedDict(PARKEYS)

class Parfile():
    """
    Class for writing, generating and handling Parfiles
    """
    def __init__(self, path='./', par=None):

        self.path = path
        self.data = self.init_parfile()
        if isinstance(par, str):
            if os.path.isfile(os.path.join(path,par)):
                self.update_fromfile(os.path.join(path,par))
            else:
                if VERBOSE:
                    print('WARNING: File {} not found'.format(par))
        elif isinstance(par, dict):
            self.update_fromdict(par)
        else:
            if VERBOSE:
                print("WARNING: Parfile is initialized empty")

    def info(self):
        """
        Print info on the parfile keys
        """
        print('\n *** Parfile ***\n')
        for key, val in PARKEYS.items():
            print('{} : \t{}.'.format(key,val[0]))

    def write_params_table(self):
        with open('parameters.md', 'w') as f:
            f.write(HEAD)
            for key, val in PARKEYS.items():
                f.write("| `{}`\t| {}\t | {}\t | {}\t | |\n".format(key,val[0],val[1],val[2]))

    def init_parfile(self):
        """
        Initialize a Parfile with all the keys
        """
        return dict.fromkeys(PARKEYS.keys())

    def update_fromfile(self,fname):
        """
        Update Par from a file
        """
        self.data.update(self.read_parfile_par(fname))

    def update_fromdict(self,dat):
        """
        Update a Parfile from a dict
        """
        self.data.update(dat)

    def read_parfile_par(self, fname):
        """
        Read a parfile.par into a py dictionary
        """
        d = {}
        with open(fname) as f:
            lines = [l for l in f.readlines() if l.strip()] # rm empty
            for line in lines:
                if line[0]=="#": continue
                line = line.rstrip("\n")
                line = line.split("#", 1)[0]
                key, val = line.split("=")
                if key.strip() not in PARKEYS.keys():
                    if VERBOSE:
                        print("WARNING: input parfile contains non-standard option {}".format(key))
                else:
                    d[key.strip()] = val.strip().replace('"','')
        return d

    def write_parfile_par(self, fname='parfile.par'):
        """
        Write a parfile to file
        """
        d = self.remove_keys_None_vals()
        with open(fname, 'w') as f:
            for key, val in d.items():
                if key not in PARKEYS.keys():
                    if VERBOSE:
                        print("WARNING: skip writing option {}".format(key))
                    continue
                if re.match('^[a-z]+', str(val), re.IGNORECASE):
                    f.write('%s = "%s"\n' % (key, str(val)))
                else:
                    f.write('%s = %s\n' % (key, str(val)))
        return

    def remove_keys_None_vals(self):
        return {k: v for k, v in self.data.items() if v is not None}

if __name__ == "__main__":

    # Parse input
    parser = argparse.ArgumentParser(description="Handle C parfiles and parameters")
    parser.add_argument("--examples", dest="examples", type=int, default=0, 
                        help="Generate all example parfiles in parfiles.json (default: 0)")
    parser.add_argument("--tests", dest="tests", type=int, default=0, 
                        help="Generate all test parfiles in parfiles.json (default: 0)")
    parser.add_argument("--exbns", dest="ex_bns", type=int, default=0, 
                        help="Load a general default BNS parfile, then use another dictionary to generate simulations and write the parfiles (default: 0)")
    parser.add_argument("--docs", dest="docs", type=int, default=0, 
                        help="Generate a parameters.md documenting the various parameters (default: 0)")
    args = parser.parse_args()

    if (args.examples):
        # Generate all the example parfiles also present under C/par
        with open("./parfiles.json", "r") as f:
            data = json.load(f)
            
        for s in data['examples']:
            print("Writing {}.par ".format(s['name']))
            par = Parfile(par=s)
            par.write_parfile_par(fname=s['name']+".par")

    if (args.tests):
        # Generate all the test parfiles
        with open("./parfiles.json", "r") as f:
            data = json.load(f)

        for s in data['tests']:
            print("Writing {}.par ".format(s['name']))
            par = Parfile(par=s)
            par.write_parfile_par(fname=s['name']+".par")

    if (args.ex_bns):
        # Load a general default BNS parfile, then use another dictionary to
        # generate simulations and write the parfiles
        with open("./parfiles.json", "r") as f:
            data = json.load(f)

        bns = data['tests'][1]
        par = Parfile(par=bns)
        d = {'use_fmode': "yes",
             'tides': "TEOBRESUM3"}
        par.update_fromdict(d)
        par.write_parfile_par(fname=bns['name']+"_updated.par")

    if(args.docs):
        # Generate the documentation for the parfile and parameters. 
        # Do: python3 parfile.py --docs 1
        VERBOSE = 0
        par = Parfile()
        par.write_params_table()