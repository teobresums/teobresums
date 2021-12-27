# TEOBResumS: Effective-one-body model for compact binary waveforms

To start, checkout the [Wiki](https://bitbucket.org/eob_ihes/teobresums/wiki/browse/).

## Code versions

Code versions are tagged on master. Code versions used in scientific
papers are tagged with the arXiv number.

For detailed release notes of the most recent version see the
[`CHANGES file`](https://bitbucket.org/eob_ihes/teobresums/src/master/CHANGES). 

 * [![C version](https://img.shields.io/badge/C_release-v3.0-green.svg)]() `v3.0-GIOTTO` quasi-circular precessing models for BBH and BNS with higher modes, including EOB-SPA.
 * [![C version](https://img.shields.io/badge/C_release-v2.0-orange.svg)](https://bitbucket.org/eob_ihes/teobresums/src/v2.0/C/) `v2.0` quasi-circular non-precessing models for BBH and BNS with higher modes.
 * [![C version](https://img.shields.io/badge/C_release-v1.1-orange.svg)](https://bitbucket.org/eob_ihes/teobresums/src/v1.1/C/) `v1.0` quasi-circular non-precessing model for BBH and BNS, including post-adiabatic method and higher modes for BNS. This version was ported in [LAL](https://git.ligo.org/lscsoft/lalsuite) and the LAL implementation is reviewed by LVC.
 * [![C++ version](https://img.shields.io/badge/C++_release-v0.1-red.svg)](https://bitbucket.org/eob_ihes/teobresums/src/v1.1/C++/) C++ implementation abandoned in June 2018 and code development is frozen. `v0.0` was reviewed by LVC. Review was concluded on 11.04.2018. The review's final statement can be found [here](https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/TEOBResumS/FinalReviewStatement). The reviewed version is tagged as `LVC-review-end `.

## Authors 

See the [`AUTHORS file`](https://bitbucket.org/eob_ihes/teobresums/src/master/AUTHORS).

## License

TEOBResumS is a free software distributed under the terms of the GNU General Public License as published by the Free Software Foundation.

See the [`LICENSE file`](https://bitbucket.org/eob_ihes/teobresums/src/master/LICENSE).

## Pre-requisites

 * C compiler (e.g. GCC)
 * [GSL library](https://www.gnu.org/software/gsl/)

## Compilation

`TEOBResumS` can be built as an executable or a python module using [disutil](https://docs.python.org/3/library/distutils.html).

To build the executable, check the compiler and compilation options in `C/Makefile` and then:

```
$ cd /PATH/TO/teobresums/C/
$ make
```

Alternatively, to build the python module under the `teobresums/Python` directory do:

```
$ cd /PATH/TO/teobresums/Python/
$ python TEOBResumSWrap_setup.py build_ext --inplace
```

To build the module in your current environment do:

```
$ cd /PATH/TO/teobresums/Python/
$ python TEOBResumSWrap_setup.py install
```

We strongly advise users to use [anaconda](https://www.anaconda.com/products/individual-d).

## Usage

Run the executable:

```
$ cd /PATH/TO/teobresums/C/
$ ./TEOBResumS.x -p <parfile>
```

Example parfiles can be found in the directory `C/par/`.

Example scripts for using the python module can be instead found in `Python/Examples/`.

If you have a [PyCBC](https://pycbc.org/) installation, a plugin for `TEOBResumS` python module is available in `PyCBC/` and should be automatically detected by `PyCBC`; 
see `PyCBC/test.py`.

## References

Main reference

 * Nagar et al (2018) [Time-domain effective-one-body gravitational waveforms for coalescing compact binaries with nonprecessing spins, tides and self-spin effects](https://inspirehep.net/record/1676430)

**please cite the above reference** when using the code. Please, see
  also the complete list of [references on the Wiki](https://bitbucket.org/eob_ihes/teobresums/wiki/References) and cite papers according to what you use.

## Developers

Please contact Sebastiano Bernuzzi and Alessandro Nagar if you wish to be included in the list of developers.

A short description of the code can be found in [here](C/src/README.md). Please, take a moment to understand the data structures and how the parameters are managed.

Since June 2018 the repository is managed (roughly) following a `master`/`development` workflow, see e.g. the description [here](https://bitbucket.org/eob_ihes/teobresums/wiki/EOB-repo-usage).

If you want to start a new branch or development make sure to

 * Talk to us 
 * Run the benchmark *before* applying any change 
 * Run the benchmark *after* you applied your changes and before pushing them; compare and make sure you did not break anything.
