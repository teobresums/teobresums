# Summary of code validation 

S.Ackay, S.Bernuzzi

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

  * Dynamical variables
  * Multipolar waveforms (all multipoles up to l=8)

Various plots showing relative differences can be found in this directory.

TODO: Add statements for NQC and ringdown

During the validation process we found minor bugs in the C++ implementation

   * Mistake in the second derivative of the metric A function.
     The mistake affect only the initial data for tidal and nonspinning runs.     
   * Memory bug in the routine computing the finite difference
     The bug does not affect any physical result
   * ...
   

