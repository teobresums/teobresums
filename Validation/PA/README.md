# Validation of PA ("rush")

Used version in commit c1883829eeb6cf0b9af807c466e8c7b402b914e7  (2 November 2018) for these checks


This report is in regards to minor bug fix in the waveform output interpolation routine for when using PHYSICAL units whether for RUSH (post-adiabatic, henceforth PA) or no-rush (regular ODE solving) runs. 

## Sanity/consistency check

First, since the bug was related to physical units, I prepared REFERENCE waveforms from 8 GEOMETRIC unit runs, 4 BNS and 4BBH,  BEFORE the bug fix, and compared these with the corresponding 8 waveforms POSTbug fix to indeed confirm the machine precision agreement that we expected. This was my "sanity check".

## Geometric vs physical units check 

Second, I wanted to make sure that the code produced consistent results between a given geometric run and its physical-unit counterpart.
To this end, I prepared 10 geometric-unit par files, 5 BNS and 5 BBH, and 10 in physical units corresponding to these. 
For example, the par file bbh_q1_s0s0_f0004's counterpart in physical units is bbh_q1_s0s0_f40 where the initial frequency of 40Hz is for a system with total mass M  = 20Msun given in the par file. 

The 10 geometric par files and the 10 corresponding physical files are:

- 'bbh_q1_s0s0_f0005' 		bbh_q1_s0s0_f100

- 'bbh_q1_s0s0_f0004' 		bbh_q1_s0s0_f40

- 'bbh_q1_s07s07_f0004' 	bbh_q1_s07s07_f54

- 'bbh_q1_s098s098_f0004' 	bbh_q1_s098s098_f23

- 'bbh_q3_s07s07_f0004' 	bbh_q3_s07s07_f20

- 'bns_q1_s0s0_f0004' 		bns_q1_s0s0_f300

- 'bns_q1_s0s0_f00027' 		bns_q1_s0s0_f180

- 'bns_q1_s04s04_f0003' 	bns_q1_s04s04_f217

- 'bns_q1_s014s014_f0005' 	bns_q1_s014s014_f375

- 'bns_q1.5_s02s04_f0004'	bns_q1.5_s02s04_f270

For each geom-phys unit pair, I made 4 comparisons where I used the following time stepping routines for the ODE solvers:

	- uniform (for geometric) - uniform (for physical)

	- adaptive - uniform

	- uniform  adaptive

	- adaptive - adaptive

This resulted in 40 comparisons of the (2,2) mode. The agreement in each comparison between the geometric run and the physical run was very good with the only issue coming from the fact that the conversion factor G M_sun/c^3 is known to 10 significant digits. Therefore, I systemically observed that each run would start with an agreement of ~10 digits at the initial frequency then usually only agree to 3 or 4 digits (or better) at the end.

NOTE: for physical unit runs to interpolate properly at the end, 

Set	
```
	interp_uniform_grid = 1
```
either in default.par or in the individual par files (we will be standardizing these better).

## Rush vs. No rush check

Third, having built confidence that the physical unit runs work, I created 20 more par files, namely the PA versions of the 20 listed above. These simply have the suffix _postadiab. With these, I made 20 more comparisons. In each case, I compared the hpc (plus,cross) file between  a no-rush run with uniform time step and a Rush run with adaptive time step and Rmin usually set to 12 (sometimes 10 or 15). These comparisons have been added to the repo under `Figures/` as `.png` figures and have been labelled by the par file names. Each figure contains the physical paramaters and the value of Rmin as well.

	BLUE curves are the no-rush waveforms

	RED are the rush run waveforms

In each plot, the top panel shows the waveforms for the entire run, at the beginning and at the end. The bottom panel shows the relative difference in the waveform, the relative difference in the amplitude and the phase difference where I computed the phase from ArcTan[hx/h+].

The agreement is usually very good. However, there are some weird behaviours that the PA experts might wish to check. 
For example, in the bns_q1_s014s014_f375 case, the rush run ends a few cycles before the no-rush run, but its geometric-unit counterpart, i.e., bns_q1_s014s014_f0005, exhibits no such behaviour.
There also seems to be a visible jump in certain runs at Rmin, where the PA run ends. The waveform for r(t) < Rmin gives a worse agreement with the norush case than the PA for the r(t) > Rmin part of the evolution. See, e.g., bns_q1_s0s0 runs in both geometric and physical units.

