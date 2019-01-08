# TEOBResumS: Effective-one-body model with spin and tidal interactions

To start, checkout the [WIKI](https://bitbucket.org/eob_ihes/teobresums/wiki/browse/).

## Code versions

Code versions are tagged on master

 * [![C version](https://img.shields.io/badge/C_release-v1.0-green.svg)](https://bitbucket.org/eob_ihes/teobresums/src/v1.0/C/) C implementation currently maintained and developed.
 * [![C++ version](https://img.shields.io/badge/C++_release-v0.1-orange.svg)](https://bitbucket.org/eob_ihes/teobresums/src/v1.0/C++/) C++ implementation abandoned in June 2018 and code development is frozen. `v0.0` was reviewed by LVC. Review was concluded on 11.04.2018. The review's final statement can be found [here](https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/TEOBResumS/FinalReviewStatement). The reviewed version is tagged as `LVC-review-end `.

## Authors

See the [`AUTHORS`](https://bitbucket.org/eob_ihes/teobresums/src/master/AUTHORS) file.

## License

See the [`LICENSE`](https://bitbucket.org/eob_ihes/teobresums/src/master/LICENSE) file.

```
TEOBResumS is a free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
 
TEOBResumS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
```

## Pre-requisites

 * C standard library
 * GSL
 * libconfig

## Compilation

```
$ export TEOBRESUMS=/PATH/TO/teobresums/C/
$ cd $TEOBRESUMS
$ make -f Makefile.TEOBResumS
```

## Usage

```
$ cd $TEOBRESUMS
$ ./TEOBResumS.x <parfile>
```

For example parfiles see `par/`

## References

Main reference

 * Nagar et al (2018) [Time-domain effective-one-body gravitational waveforms for coalescing compact binaries with nonprecessing spins, tides and self-spin effects](https://inspirehep.net/record/1676430)

**please cite the above reference** when using the code.

Other key references

 * Nagar, Messina, Rettegno, Bini, Damour, Geralico, Akcay, Bernuzzi (2018) [Nonlinear-in-spin effects in effective-one-body waveform models of spin-aligned, inspiralling, neutron star binaries](https://inspirehep.net/record/1710050)
 * Akcay, Bernuzzi, Messina, Nagar, Ortiz, Rettegno [Effective-one-body multipolar waveform for tidally interacting binary neutron stars up to merger](http://inspirehep.net/record/1707624) 
 * Nagar, Riemenschneider, Pratten (2017) [Impact of Numerical Relativity information on effective-one-body waveform models](http://inspirehep.net/record/1518392)
 * Del Pozzo, Nagar (2016) [Analytic family of post-merger template waveforms](http://inspirehep.net/record/1469053)
 * Nagar, Damour, Reisswig, Pollney (2015) [Energetics and phasing of nonprecessing spinning coalescing black hole binaries](http://inspirehep.net/record/1380155) NRAR comparison, BBH spin
 * Bernuzzi, Nagar, Dietrich, Damour (2014) [Modeling the Dynamics of Tidally Interacting Binary Neutron Stars up to the Merger](http://inspirehep.net/record/1334334) TEOBResum
 * Bini, Damour (2014) [Gravitational self-force corrections to two-body tidal interactions and the effective one-body formalism](http://inspirehep.net/record/1318823)
 * Damour, Nagar (2014) [New effective-one-body description of coalescing nonprecessing spinning black-hole binaries](http://inspirehep.net/record/1303216) New formalism for spin with centrifugal radius
 * Damour, Nagar (2014) [A new analytic representation of the ringdown waveform of coalescing spinning black hole binaries](http://inspirehep.net/record/1298802)
 * Bini, Damour, Faye (2012) [Effective action approach to higher-order relativistic tidal interactions in binary systems and their effective one body description](http://inspirehep.net/record/1089309)
 * Damour, Nagar, Bernuzzi (2012) [Improved effective-one-body description of coalescing nonspinning black-hole binaries and its numerical-relativity completion](http://inspirehep.net/record/1207869)
 * Bernuzzi, Nagar, Thierfelder, Bruegmann (2012) [Tidal effects in binary neutron star coalescence](http://inspirehep.net/record/1114723) NRAR comparison, NNLO tidal model
 * Damour, Nagar (2009) [Effective One Body description of tidal effects in inspiralling compact binaries](http://inspirehep.net/record/838179) EOB formalism for tides
 * Damour, Nagar (2009) [Relativistic tidal properties of neutron stars](http://inspirehep.net/record/821786) Formalism relativistic Love numbers
 * Damour, Nagar (2009) [An Improved analytical description of inspiralling and coalescing black-hole binaries](http://inspirehep.net/record/812296)
 * Damour, Iyer, Nagar (2008) [Improved resummation of post-Newtonian multipolar waveforms from circularized compact binaries](http://inspirehep.net/record/802497) Resummed EOB waveform

## Developers

Please contact Sebastiano Bernuzzi and Alessandro Nagar if you wish to be included in the list of developers.

A short description of the code can be found in [here](C/src/README.md). Please, take a moment to understand the data structures and how the parameters are managed.

Since June 2018 the repository is managed (roughly) following a `master`/`development` workflow, see e.g. the description [here](https://bitbucket.org/eob_ihes/teobresums/wiki/EOB-repo-usage).

If you want to start a new branch or development make sure to

 * Talk to us 
 * Run the benchmark *before* applying any change 
 * Run the benchmark *after* you applied your changes and before pushing them; compare and make sure you did not break anything.

