# Summary of code validation 

## S.Akcay, S.Bernuzzi

UPDATE as of July:

The C code now contains new terms that are not in the C++ as we have ceased to support it.

The last time the C/C++ code were checked to be FULLY consistent was on 18 June 2018 by comparing the dynamics and the inspiral+merger+ringdown waveforms 9 BBH cases and dynamics+inspiral waveform for 6 BNS cases.

The checked versions are:

- C on branch `development`     commit  2ac8bb48a5804010c8951a6e23904d89117cf62d    on 18 June 2018

- C++ on branch `master`        commit e97b08cb43d25b23cc9524c46ed74e7b9a1e2c64     on 18 June 2018



--------------------------------------------------------



The new C implementation on branch `development`

    commit 8e64500adef87e1f1b969af8c3c9c6d39bccbb0b

is tested against the previous C++ implementation on branch `master`

    commit 1fde4f901af43e724c43959faed89509b7d1435c

We validate by comparing the output of the C code with that of the C++ code.
The results are in the directory `../validation` under the branch `development`

More specifically, we computed the relative disagreement between the two codes " \Delta_{rel}"
for: 

 - 7 dynamical quantities: {r, phi, p_phi, M*Omega, \ddot{r}, p_r*, M*Omega_orb} during the inspiral
 - Waveforms h_lm in both amplitude and phase for all non-zero {l,m} modes during the inspiral
 - Waveform for the (2,2) mode in both amplitude and phase including NQC and RINGDOWN.

We performed our comparison for a total of 11 different binary black hole (BBH/bbh) and 6 different binary neutron star (BNS/bns) runs each of which is parametrized in terms of

	 (q, Sz^A, Sz^B, f_min)

where q >= 1 is the mass ratio, s1 and s2 are the dimensionless spins of the compact objects, and f_min is the dimensionless initial Keplerian frequency.

We chose the following values for the above parameters in our comparisons.

  * BBH runs of mass ratios q = 1, 3, 12 non-spinning
  * BBH runs of mass ratios q = 1, 2, 3, 5 with various spins aligned or anti-aligned with the orbital angular momentum
	- q = 1, 3 with equal spins Sz^A = Sz^B = +0.7 
	- q = 1 with Sz^A = Sz^B = +0.98 with f_min = {1e-3, 2e-3, 4e-3}
	- q = 1 with Sz^A = Sz^B = -0.99 with f_min = 4e-3
	- q = 2 with Sz^A = -0.5, Sz^B = -0.8, f_min = 3e-3
	- q = 5 with Sz^A = 0.6, Sz^B = -0.9, f_min = 3e-3

  * BNS runs of mass ratios q = 1, 2 non-spinning
  * BNS runs of mass ratios q = 1, 1.5, 2
	- q = 1 with Sz^A = Sz^B = 0.4 with f_min = 3e-3
	- q = 2 with Sz^A = Sz^B = 0.2 with f_min = 2e-3
	- q = 1.5 with Sz^A = -0.2, Sz^B = 0.3, f_min = 2.5e-3
	- q = 1.5 with Sz^A = 0.2, Sz^B = 0.4, f_min = 5e-4

The parfiles for these runs can be found in the directory `../benchmark/`.

We repeated each comparison three times because both codes have "uniform" or "adaptive" time stepping routines in their'respective ODE integrators. Therefore, the three cases we considered are

 - Adaptive(C++) vs. Adaptive(C) labelled as "Adaptive" in the subfolders.
 - Adaptive(C++) vs. Uniform(C) labelled as "Adaptive_Uniform" in the subfolders.
 - Uniform(C++) vs. Uniform(C) labelled as "Uniform" in the subfolders.

As we considered the Adaptive_Uniform case to be representative of the Uniform_Adaptive case, we did not repeat the comparisons for the latter. We had an additional motivation for this as C++ code only seems to perform the ringdown computations for the "adaptive" case.

We should remark that for the two "low"-frequency cases: 
	1) BBH (1,0.98,0.98,0.001) and 
	2) BNS (1.5, 0.2, 0.4, 5e-4) 
we used ONLY the adaptive ODE routine since uniform timestepping takes too long and generates gigabytes of data.

Because of this "caveat" the "adaptive" folders contain 17 (= 11 BBH+6 BNS) cases whereas the rest contain 10 BBH + 5 BNS = 15 cases.

Each case is labelled by 
		`bbh_q<x>_s<y>_s<z>_f<w>`
		`bns_q<x>_s<y>_s<z>_f<w>`
where {x, y, z, w} are the values for {q, Sz^A, SZ^B, f_min}. 

The validation process shows agreement between the two codes with relative differences at or below the level <~ 1e-6 or much smaller in almost every single comparison we present.

There are cases in the multipolar waveform comparisons in which the relative disagreement exceeds the 1e-3 threshold. These happen only for high modes and are due to the fact that the absolute magnitudes of these modes are becoming comparable to the 1e-11 ODE integrator accuracy. These relative differences can be decreased by increasing the accuracy of the ODE solvers in both C and C++ at the expense of runtime. Part of this disagreement is due to the external Mathematica interpolator we used to interpolate the C++ data in the comparisons of Adaptive and Adaptive_Uniform cases where the ODE integrators use different time steps. In the case of Uniform timestepping, there is no need for such interpolation hence the disagreements are always <~ 1e-5 except maybe at the last few time steps.

For the NQC checks, we used a smaller set as ringdown checks serve as a proxy for NQC.

For the ringdown, we compared only the (2,2) mode as this is what is currently available from the C++ code. As an extra check, we once again compared the waveform from t=0 until the end of the ringdown. We additionally zoomed in on the attachment point of the ringdown and showed the amplitude, phase comparisons in subfolders called "Attachment_22". As ringdown only occurs for BBHs, our comparison set contains 11 cases for "Adaptive" routines and 10 for "Adaptive_Uniform".

There is a clear disagreement between the two codes during ringdown, but in all our comparisons we found this to be less than 1e-3.
We discovered that this disagreement is due to a numerical cancellation in the ringdown routine that amplifies a <~ 1e-6 relative disagreement to nearly 1e-3 in some cases. This <~1e-6 disagreement comes out of the function that finds the time of the peak frequency. We believe the routine in the C code is more accurate and confirmed this to be mostly true using Mathematica. We further confirmed that this is the cause of the disagreement by forcing the C code to use C++'s result for this peak time which resulted in the agreement of the ringdown waveform amplitudes to better than 1e-12.

During the validation process we found minor bugs in the C++ implementation:

   1. Mistake in the computation of the second derivative of the metric A function.
      The mistake affects only the initial data for tidal and nonspinning runs.
      For comparative tests we imposed manually the values of A'' in the C++ code.
      We did not fix the bug in C++.
     
   2. Memory bug in the routine computing the finite difference.
      The bug does not affect any physical result.
     
   3. Multipolar spinning waveform in (2,1),(3,1) and (3,3) modes.
      Abuse of `const double` code, resulting in wrong values of `p2`
      and `p3` in Newtonian waveform for spin case. 
      The bug does not affect results using just the (2,2) mode

   4. Waveform was not computed in first step, resulting in small time
      shift after the interpolation (if performed).
      Bug has been fixed.
      
   5. Amplitudes rholm (3,1) and (3,3) did not contain the nu-dependent corrections.
      These corrections have been added.

   6. (8,8) mode had been labelled as the (8,1) mode and the remaining l=8 modes were missing both in the Newtonian flux and the Newtonian waveform functions. These have been added.


