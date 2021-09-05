# TEOBResumS: Effective-one-body model with spin and tidal interactions

To start, checkout the [WIKI](https://bitbucket.org/eob_ihes/teobresums/wiki/browse/).

## Code versions

Code versions are tagged on master. 
For detailed release notes of the most recent version see the [`CHANGES`](https://bitbucket.org/eob_ihes/teobresums/src/master/CHANGES) file.

 * [![C version](https://img.shields.io/badge/C_release-v2.0-green.svg)](https://bitbucket.org/eob_ihes/teobresums/src/v2.0/C/) C implementation currently maintained and developed.
 * [![C version](https://img.shields.io/badge/C_release-v1.1-orange.svg)](https://bitbucket.org/eob_ihes/teobresums/src/v1.1/C/) Previous C implementation. `v1.0` was ported to [LAL](https://git.ligo.org/lscsoft/lalsuite) and that implementation is reviewed by LVC.
 * [![C++ version](https://img.shields.io/badge/C++_release-v0.1-red.svg)](https://bitbucket.org/eob_ihes/teobresums/src/v1.1/C++/) C++ implementation abandoned in June 2018 and code development is frozen. `v0.0` was reviewed by LVC. Review was concluded on 11.04.2018. The review's final statement can be found [here](https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/TEOBResumS/FinalReviewStatement). The reviewed version is tagged as `LVC-review-end `.

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

 * C compiler (e.g. GCC)
 * [GSL library](https://www.gnu.org/software/gsl/)
 * [libconfig](http://hyperrealm.github.io/libconfig/)

## Compilation

To build an executable, check compiler and compilation options in `Makefile.TEOBResumS` and then:

```
$ cd /PATH/TO/teobresums/C/
$ make -f Makefile.TEOBResumS
```

Alternatively, you can build a python module using [disutil](https://docs.python.org/3/library/distutils.html) and the following script  

```
$ cd /PATH/TO/teobresums/Python/
$ python TEOBResumSWrap_setup.py build_ext --inplace
```

to build the module under the teobresums/Python directory, or  

```
$ cd /PATH/TO/teobresums/Python/
$ python TEOBResumSWrap_setup.py install
```

to build it in your environment.  

We strongly advise users to use [anaconda](https://www.anaconda.com/products/individual-d) when building the Python module.

If you run in compilation errors, you may need to export the following variables:

```
$ export C_INCLUDE_PATH=PATH/TO/libconfig/include
$ export LIBRARY_PATH=$LIBRARY_PATH:PATH/TO/libconfig/lib
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:PATH/TO/libconfig/lib
```

## Usage

```
$ cd $TEOBRESUMS
$ ./TEOBResumS.x <parfile>
```

Example parfiles can be found in the directory `$TEOBRESUMS/C/par/`.

Example scripts for using the python module can be instead found in '$TEOBRESUMS/Python/Examples/`.

If you have a [PyCBC](https://pycbc.org/) installation, a plugin for `TEOBResumS` python module is available in `$TEOBRESUMS/pycbc/` and should be automatically detected by `PyCBC`; 
see `$TEOBRESUMS/pycbc/test.py`.

## References

Main reference

 * Nagar et al (2018) [Time-domain effective-one-body gravitational waveforms for coalescing compact binaries with nonprecessing spins, tides and self-spin effects](https://inspirehep.net/record/1676430)

**please cite the above reference** when using the code. Please see also the complete list of [`REFERENCES`](https://bitbucket.org/eob_ihes/teobresums/src/master/REFERENCES.md) and cite papers according to what you use.

## Developers

Please contact Sebastiano Bernuzzi and Alessandro Nagar if you wish to be included in the list of developers.

A short description of the code can be found in [here](C/src/README.md). Please, take a moment to understand the data structures and how the parameters are managed.

Since June 2018 the repository is managed (roughly) following a `master`/`development` workflow, see e.g. the description [here](https://bitbucket.org/eob_ihes/teobresums/wiki/EOB-repo-usage).

If you want to start a new branch or development make sure to

 * Talk to us 
 * Run the benchmark *before* applying any change 
 * Run the benchmark *after* you applied your changes and before pushing them; compare and make sure you did not break anything.

