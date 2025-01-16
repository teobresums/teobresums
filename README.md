# TEOBResumS: Effective-one-body model for compact binary waveforms

To start, checkout the [Wiki](https://bitbucket.org/teobresums/teobresums/wiki/browse/).

## Code versions

Code versions are tagged on master. Code versions used in scientific
papers are tagged with the arXiv number. 

For detailed release notes of the most recent version see the
[`CHANGES` file](https://bitbucket.org/teobresums/teobresums/src/GIOTTO/CHANGES). 

## Authors 

See the [`AUTHORS` file](https://bitbucket.org/teobresums/teobresums/src/GIOTTO/AUTHORS).

## License

TEOBResumS is a free software distributed under the terms of the GNU General Public License as published by the Free Software Foundation.

See the [`LICENSE` file](https://bitbucket.org/teobresums/teobresums/src/GIOTTO/LICENSE).

## Pre-requisites

 * C compiler (e.g. GCC)
 * [GSL library](https://www.gnu.org/software/gsl/)

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

## Developers

Please contact Sebastiano Bernuzzi and Alessandro Nagar if you wish to be included in the list of developers.

A short description of the code can be found in [here](C/README.md). Please, take a moment to understand the data structures and how the parameters are managed.

See [here](https://bitbucket.org/teobresums/teobresums/wiki/EOB-repo-usage) for a description of the workflow of this repo.

If you want to start a new branch or development make sure to

 * Talk to us 
 * Run the benchmark *before* applying any change 
 * Run the benchmark *after* you applied your changes and before pushing them; compare and make sure you did not break anything.

## Contact

You can get in contact with TEOBResumS developers and fellow users through our mailing list. To subscribe please click on the badge below.

[![Mailing List](https://img.shields.io/badge/mailing%20list-join-informational)](https://lserv.uni-jena.de/mailman/listinfo/teobresums)
