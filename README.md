# TEOBResumS: Effective-one-body model for compact binary waveforms

To start, checkout the [Wiki](https://bitbucket.org/eob_ihes/teobresums/wiki/browse/).

## Code versions

Code versions are tagged on master. Code versions used in scientific
papers are tagged with the arXiv number. 

For detailed release notes of the most recent version see the
[`CHANGES` file](https://bitbucket.org/eob_ihes/teobresums/src/master/CHANGES). 

### Latest version  
[![C version](https://img.shields.io/badge/TEOBResumS-v4.1.4--GIOTTO-green)](https://bitbucket.org/eob_ihes/teobresums/src/v4.1.4-GIOTTO/)   
Quasi-circular precessing and non-precessing models for BBH, BNS and BHNS with higher modes, including EOB-SPA.

## Authors 

See the [`AUTHORS` file](https://bitbucket.org/eob_ihes/teobresums/src/master/AUTHORS).

## License

TEOBResumS is a free software distributed under the terms of the GNU General Public License as published by the Free Software Foundation.

See the [`LICENSE` file](https://bitbucket.org/eob_ihes/teobresums/src/master/LICENSE).

## Pre-requisites

 * C compiler (e.g. GCC)
 * [GSL library](https://www.gnu.org/software/gsl/)

## Installation

`TEOBResumS` can be built as an executable or a python module using [disutil](https://docs.python.org/3/library/distutils.html).

### PyPI
```
$ pip install teobresums
```

### From source
```
$ git clone https://bitbucket.org/eob_ihes/teobresums.git
$ cd teobresums/Python/
$ python setup.py install
```

We strongly advise users to use [anaconda](https://www.anaconda.com/products/individual-d).

### C executable
To build the executable, check the compiler and compilation options in [`C/Makefile` file](https://bitbucket.org/eob_ihes/teobresums/src/master/C/Makefile) and then:

```
$ git clone https://bitbucket.org/eob_ihes/teobresums.git
$ cd teobresums/C/
$ make
```

## Usage

### Python extension

```
$ python
>>> import EOBRun_module
```

Example scripts for using the python module can be found in `Python/Examples/`.

If you have a [PyCBC](https://pycbc.org/) installation, a plugin for `TEOBResumS` python module is available in `PyCBC/` and should be automatically detected by `PyCBC`; 
see `PyCBC/test.py`.

### C executable
Run the executable:

```
$ cd /PATH/TO/teobresums/C/
$ ./TEOBResumS.x -p <parfile>
```

Example parfiles can be found in the directory `C/par/`.

## References
Main references:

 * Riemenschneider et al (2021) [Assessment of consistent next-to-quasicircular corrections and postadiabatic approximation in effective-one-body multipolar waveforms for binary black hole coalescences](https://inspirehep.net/literature/1858441)
 * Nagar et al (2020) [Multipolar effective one body waveform model for spin-aligned black hole binaries](https://inspirehep.net/literature/1777194)
 * Nagar et al (2019) [Multipolar effective one body model for nonspinning black hole binaries](https://inspirehep.net/literature/1730540)
 * Nagar et al (2018) [Time-domain effective-one-body gravitational waveforms for coalescing compact binaries with nonprecessing spins, tides and self-spin effects](https://inspirehep.net/record/1676430)
 * Nagar et al (2015) [Energetics and phasing of nonprecessing spinning coalescing black hole binaries](https://inspirehep.net/literature/1380155) 
 * Damour et al (2014) [New effective-one-body description of coalescing nonprecessing spinning black-hole binaries](https://inspirehep.net/literature/1303216)
 
**please cite the above references** when using the code. Please, see
  also the complete list of [references on the Wiki](https://bitbucket.org/eob_ihes/teobresums/wiki/References) and cite papers according to what you use.

## Developers

Please contact Sebastiano Bernuzzi and Alessandro Nagar if you wish to be included in the list of developers.

A short description of the code can be found in [here](C/README.md). Please, take a moment to understand the data structures and how the parameters are managed.

See [here](https://bitbucket.org/eob_ihes/teobresums/wiki/EOB-repo-usage) for a description of the workflow of this repo.

If you want to start a new branch or development make sure to

 * Talk to us 
 * Run the benchmark *before* applying any change 
 * Run the benchmark *after* you applied your changes and before pushing them; compare and make sure you did not break anything.

## Contact

You can get in contact with TEOBResumS developers and fellow users through our mailing list. To subscribe please click on the badge below.

[![Mailing List](https://img.shields.io/badge/mailing%20list-join-informational)](https://lserv.uni-jena.de/mailman/listinfo/teobresums)
