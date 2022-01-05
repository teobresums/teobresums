TEOBResumS
##########

Complete waveforms from quasicircular compact binaries mergers.
Features:

 * Time-domain inspiral-merger-ringdown waveform for BBH
 * Time-domain inspiral-merger waveform for BNS
 * Frequency-domain inspiral-merger waveform with EOB-SPA 
 * Higher modes up to (8,8)
 * Non-precessing spins dynamics and waveform 
 * Tidal interaction in dynamics and waveform
 * Spin Precession (from v3)
 * ODE or post-adiabatic solution of Hamiltonian dynamics
 
Main links:
 
 * `TEOBResumS repo <https://bitbucket.org/eob_ihes/teobresums/src/master/>`_

 * `TEOBResumS wiki <https://bitbucket.org/eob_ihes/teobresums/wiki/Home>`_

Quick start
###########

Clone the ``TEOBResumS`` repo::

   git clone https://bitbucket.org/eob_ihes/teobresums/src/master/

Compile the code and setup the ``python`` module::

  cd teobresums/Python/ # check options in setup.py
  python setup.py build_ext --inplace # see also the Makefile

``PyCBC`` should now automatically detect the package and make the waveform available for use.

May the force be with you!

Usage
#####

An example of usage in given in ``test.py``.

Please consult `TEOBResumS wiki
<https://bitbucket.org/eob_ihes/teobresums/wiki/Home>`_ and the
literature for the binary parameter space covered by the model and its accuracy.


Parameters
##########

There is a set of required parameters:

:param mass1: mass of compact object 1, in Solar masses
:type mass1: float
:param mass2: mass of compact object 2, in Solar masses
:type mass2: float
:param spin1z: z-component of dimensionless spin of compact object 1
:type spin1z: float
:param spin2z: z-component of dimensionless spin of compact object 2
:type spin2z: float
:param lambda1: quadrupolar tidal polarizability parameter of compact object 1
:type lambda1: float
:param lambda2: quadrupolar tidal polarizability parameter of compact object 2
:type lambda2: float
:param distance: source distance, in Mpc
:type distance: float
:param inclination: source inclination angle, in rad
:type inclination: float
:param mode_array: list of modes to be used, 
:type mode_array: list ot tuples (l,m), optional. Example: [(2,1),(2,2),(3,2)]

By default ``TEOBResumS`` employs:

 * Post-adiabatic dynamics + ODE solution below a minimum radius (14)
 * 2,2 mode
 * Nonprecessing spins
 * Tides (if input tidal parameters are positive)
 * ...

The default behaviour can be changed via an extended set of input
parameters. Here, some knowledge of the model and of the code might be
important. For a starting point, see the full list of input parameters
in ``teobresums_pars_update()``.

	      



