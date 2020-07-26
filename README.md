# TEOBResumS: Effective-one-body model with spin and tidal interactions

To start, checkout the [WIKI](https://bitbucket.org/eob_ihes/teobresums/wiki/EOB-repo-usage).

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

 * Standard C libs, including libconfig
 * GSL lib

## Compilation

```
$ cd $TEOBRESUMS
$ make -f Makefile.TEOBResumS
```

Alternatively, you can build a python module using [disutil](https://docs.python.org/3/library/distutils.html) and the following script
```
$ cd /PATH/TO/teobresums/Python/
$ python TEOBResumSWrap_setup.py build_ext --inplace
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

**please cite the above reference** when using the code. Please see also the complete list of [`REFERENCES`](https://bitbucket.org/eob_ihes/teobresums/src/master/REFERENCES.md) and cite papers according to what you use.
## Developers

Please contact Sebastiano Bernuzzi and Alessandro Nagar if you wish to be included in the list of developers.

A short description of the code can be found in [here](C/src/README.md). Please, take a moment to understand the data structures and how the parameters are managed.

Since June 2018 the repository is managed (roughly) following a `master`/`development` workflow, see e.g. the description [here](https://bitbucket.org/eob_ihes/teobresums/wiki/EOB-repo-usage).

If you want to start a new branch or development make sure to

 * Talk to us 
 * Run the benchmark *before* applying any change 
 * Run the benchmark *after* you applied your changes and before pushing them; compare and make sure you did not break anything.







