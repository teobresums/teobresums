import os; import re; import json; import collections

VERBOSE = 1

PARKEYS = {
    'M': 'Total mass of the system',
    'q': 'Mass ratio of the system',
    'LambdaAl2': 'Quadrupolar tidal parameter of body 1',
    'LambdaBl2': 'Quadrupolar tidal parameter of body 2',
    'chi1': 'Spin of body 1',
    'chi2': 'Spin of body 2',
    'chi1x': 'x component of the spin of body 1',
    'chi1y': 'y component of the spin of body 1',
    'chi1z': 'z component of the spin of body 1',
    'chi2x': 'x component of the spin of body 2',
    'chi2y': 'y component of the spin of body 2',
    'chi2z': 'z component of the spin of body 2',
    #
    'r0': 'Initial separation in mass-rescaled units',
    'initial_frequency' : 'Initial frequency in Hz or mass-rescaled units, depending on use_geometric_units',
    'use_geometric_units': 'Use mass rescaled units or physical units',
    #
    'use_spins': '0 = nonspinning (deprecated), 1 = spin-aligned, 2 =precessing spins',
    'tides': 'Tidal model (NNLO, TEOBRESUM, TEOBRESUM3) ',
    'tides_gravitomagnetic': 'Use gravitomagnetic tidal coefficients',
    'use_fmode': 'Use the fmode resonance model',
    'use_lambda234_fits': 'Use quasi-universal relations to fit Lambda4,5,6,7,8 (YAGI13, GODZIEBA20)',
    'use_speedytail': 'Speed up the computation of the tail factor in h_lm and F_lm',
    'pGSF_tidal': 'GSF p-exponent (default=4)',
    'dt_merger_interp': 'dt used for the interpolation of the merger dynamics, after the peak of Omega',
    'dt_interp': 'dt used to interpolate the waveform',
    'dt': 'dt used if solving the ODEs with a uniform timestep',
    'srate_interp': 'Sampling rate, 1/dt at which to interpolate the waveform',
    'srate': 'Sampling rate, 1/dt at which the waveform is computed if using a uniform timestep',
    'size': 'Initial guess for the size of the waveform',
    'ringdown_extend_array': 'Extend the ringdown array',
    'distance': 'Distance of the source from Earth [Mpc]',
    'inclination': 'Angle between \hat{N} and the orbital angular momentum at the initial time',
    'polarization': 'Polarization of the waves (unused)',
    'coalescence_angle': 'Reference coalescence phase',
    'use_mode_lm': 'List of modes to use',
    'output_lm': ' List of modes to output',
    #
    'postadiabatic_dynamics': 'Turn on/off the postadiabatic approximation',
    'postadiabatic_dynamics_N': 'Use the N-th order PA (default=8)',
    'postadiabatic_dynamics_size': 'Guess for the size of the postadiabatic dynamics',
    'postadiabatic_dynamics_rmin': 'Stop the postadiabatic dynamics at this value of separation (default=14)',
    'postadiabatic_dynamics_dr': 'Step dr of the postadiabatic dynamics',
    'postadiabatic_dynamics_stop':'Stop after the postadiabatic dynamics',
    #
    'centrifugal_radius': 'Model for the centrifugal radius',
    'use_flm': 'Model for the radiation reaction',
    'nqc': 'Choose how to set the NQCs (manual, auto, no)',
    'nqc_coefs_flx': 'Option for the NQC model used in the flux',
    'nqc_coefs_hlm': 'Option for the NQC model used in the waveform',
    'nqc_coefs_flx_file' : 'File containing the NQC coefficients for the flux',
    'compute_LR': 'Compute the light ring (LR)',
    'compute_LR_guess': 'Initial guess for the light ring (LR)',
    'compute_LSO': 'Compute the last stable orbit (LSO)',
    'compute_LSO_guess': 'Guess for the last stable orbit (LSO)',
    'project_spins': 'Project the spins along the dynamics',
    'spin_interp_domain': 'Evaluate spins(time) or spins(frequency)',
    'spin_flx': 'Model for the spin-evolution \dot\omega',
    'ringdown_eulerangles': 'Option for the extension of the Euler angles beyond merger',
    #
    'output_hpc': 'Output the plus and cross polarizations',
    'output_multipoles': 'Output the waveform multipoles h_lm',
    'output_dynamics': 'Output the EOB dynamics',
    'output_nqc': 'Output the NQCs',
    'output_nqc_coefs': 'Output the NQC coefficients',
    'output_ringdown': ' Output the ringdown portion of the waveform',
    #
    'ode_timestep': 'Choose the ODE timestep between adaptive, uniform, adaptive+uniform after LSO',
    'ode_abstol': 'ODE absolute numerical tolerance',
    'ode_reltol': 'ODE relative numerical tolerance',
    'ode_tmax': 'Time stopping condition for the ODE solver',
    'ode_stop_at_radius':'Stop the ODE evolution at r',
    'df': 'Frequency spacing for TEOBResumSPA',
    'time_shift_FD': 'Time shift the FD waveform (merger at zero)',
    'interp_freqs': 'Use a user-input list of frequency to interpolate TEOBResumSPA',
    'freqs': 'List of frequencies required by interp_freqs',
    'output_dir': 'directory where to perform the output',
    "interp_uniform_grid": 'Interpolate on a uniform grid?',
    }


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
            print('{} :\t{}.'.format(key,val))

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

    if (1):
        # Generate all the example parfiles also present under C/par
        with open("./parfiles.json", "r") as f:
            data = json.load(f)
            
        for s in data['examples']:
            print("Writing {}.par ".format(s['name']))
            par = Parfile(par=s)
            par.write_parfile_par(fname=s['name']+".par")

    if (1):
        # Generate all the test parfiles
        with open("./parfiles.json", "r") as f:
            data = json.load(f)

        for s in data['tests']:
            print("Writing {}.par ".format(s['name']))
            par = Parfile(par=s)
            par.write_parfile_par(fname=s['name']+".par")

    if (0):
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
