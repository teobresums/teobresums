# Summary of code validation 

S.Akcay, S.Bernuzzi

--- WORK IN PROGRESS ---

The new C implementation on branch `development`

    commit ???

is tested against the previous C++ implementation on branch `master`  

    commit bc8c1a7e49b6ecbb098eebff51678dfd3acf2112

Tests

  * BBH runs of mass ratios q = 1, 3, 10 nospinning
  * BNS runs of mass ratios q = 1, 2 nospinning
  * BBH runs of mass ratios q = 1, 10 and equal spins Sz^A = Sz^B = +0.7 
  * BBH runs of mass ratios q = 1, 10 and unequal spins Sz^A = -0.7 , Sz^B = +0.7
  * BBH runs of mass ratios q = 1, 10 and unequal spins Sz^A = 0 Sz^B = 0.7
  * BBH runs of mass ratios q = 1, 10 and unequal spins Sz^A = -0.1 Sz^B = 0.7
  * BNS runs of mass ratios q = 1, 2 and equal spins Sz^A = Sz^B = +0.4
  * BNS runs of mass ratios q = 1, 2 and unequal spins Sz^A = -0.3, Sz^B = +0.4

The parfiles for these runs can be found in the directory `../benchmark/`.

The validation process shows agreement between the two codes with relative differences at or below the level <~ 1e-6 in 

  * Dynamical variables (file names dyn*.txt )
  * Multipolar waveforms (all multipoles up to l=8, file name Waveform_ ...)

Various plots showing relative differences can be found in this directory.

TODO: Add statements for NQC and ringdown

During the validation process we found minor bugs in the C++ implementation

   1. Mistake in the computation of the second derivative of the metric A function.
      The mistake affects only the initial data for tidal and nonspinning runs.
      For comparative tests we imposed manually the values of A'' in the C++ code.
     
   2. Memory bug in the routine computing the finite difference.
      The bug does not affect any physical result.
     
   3. Multipolar spinning waveform in (2,1),(3,1) and (3,3) modes.
      Abuse of `const double` code, resulting in wrong values of `p2`
      and `p3` in Newtonian waveform for spin case. 
      The bug does not affect results using just the (2,2) mode

   4. Waveform was not computed in first step, resulting in small time
      shift after the interpolation (if performed). 
      
   5. Waveform for spin and tides employs nonspinning
      Hamiltonian. This result in differences of order ~ 1e-5 in 
      dynamical variables that remain of that order due to
      cancellations, but it larger in variable ddotr. 

All the above bugs except 1. have been fixed in several commits up to

    commit XXX


