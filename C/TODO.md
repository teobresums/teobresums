# TODO list

## Spin dynamics

Urgent:

 - [ ] PARS: fix EOB pars where needed (see //FIXME comments)
 - [ ] RHS: check PN coefs
 - [ ] RHS: check dot gamma angle (enforce >0?)
 - [ ] RHS: check q convention in  https://arxiv.org/abs/1307.4418, > 1 ? < 1. Now assumed ours, q>1
 - [ ] RHS: check definition of v  = omg^{1/3}, with or without M ?
 - [ ] RHS: add N4LO pieces (now NLO), Current N4LO equations in paper do not look good
 - [ ] RHS: Compute PN coefs at first call (as in other part of the code)
 - [ ] ODE stopping criterion (Momega = ?, should be sufficiently large to allow interpolation)
 - [ ] Twist: check mapping real/imag <---> phase/ampli
 - [ ] interface with main: add routines to update the spin component along L in the appropriate places of the dynamics and aligned spin waveforms computation.

Later:

 - [ ] RHS: Write spin equations in mass-rescaled variables ...
 - [ ] Twist: add FD routine
 - [ ] RHS: add rhs in the angles only
 - [ ] RHS: add EOB rhs 
 - [ ] add merger ...
 - [ ] stand alone test code

Note for the logic:

 - use " usespin == MODE_SPIN_GENERIC " to switch/call/activate generic spin options
 - spin dynamics should be called/computed at the beginning of main if required.
 - during aligned-spin dynamics and waveform, add a call to " eob_spin_dyn_Sp_interp " to compute the projection chi1 and chi2 at the given time
 - before h+,x computation, allocate the mem for the twisted modes hTlm and twist the aligned hlm with twist_hlm_*

## Other

 - [x] remove parameters from Dynamics structure to void duplication, use EOBPars
 - [x] move compute_hpc in *Waveform.c 
