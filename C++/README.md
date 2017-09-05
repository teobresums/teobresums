# Quick start 

## Requirements

- C++ compiler (GNU, Intel)
- GSL/BLAS 

## How to compile

There are 2 options

- simple Makefile
- build that includes Python wrapper

### Makefile

A simple `Makefile` is provided in the main directory. Just run 

```
$ make -f Makefile.TEOBResumS
```

### Build and Python wrapper

Requirements

- autoreconf 2.69
- libtool
- automake 1.15.1

To compile 

```
$ ./bootstrap
$ ./configure --prefix=</abs/path/to/your/build/> CXXFLAGS="-std=c++11 -O3 -Wall -fPIC"
$ make
$ make install
```

To build and install the python wrapper

```
$ cd python 
$ python setup.py build_ext --inplace 
```

Once the installation is finished, add 

`/abs/path/to/your/build/`

to the path

## How to run 

Usage/help 

```
$ ./TEOBResumS.x 

USAGE:
        ./TEOBResumS.x -p <parfile>
        ./TEOBResumS.x [OPTIONS]

OPTIONS:
        -p                   <parfile>  reads input parameters from parfile. Overrides all other arguments. []
        -m1                  <double>   mass of the primary [Msun]. [40]
        -m2                  <double>   mass of the secondary [Msun].  [40]
        -chi1                <double>   dimensionless spin component along the orbital angular momentum of the primary. [0]
        -chi2                <double>   dimensionless spin component along the orbital angular momentum of the secondary. [0]
        -distance            <double>   source distance [Mpc]. [100]
        -inclination         <double>   (IOTA) inclination angle [rad]. [0]
        -polarisation        <double>   (PSI) polarisation angle [rad]. [0]
        -f_min               <double>   starting frequency [Hz]. [20]
        -srate               <double>   sampling rate [Hz]. [4096]
        -lambda1             <double>   tidal deformability for body 1 (Lambda/M^5). Only if tidal corrections are enabled. [0]
        -lambda2             <double>   tidal deformability for body 2 (Lambda/M^5). Only if tidal corrections are enabled. [0]
        -tidal               <int>      enable tidal corrections. [0 (false)]
        -NQC                 <int>      enable NQC corrections. [1 (true)]
        -speedy              <int>      faster tails calculations. [1 (true)]
        -dynamics            <int>      output dynamics evolution. [0 (false)]
        -RW                  <int>      Regge-Wheeler-Zerilli potential. [0 (false)]
        -multipoles          <int>      enable single multipole output, in geometrical units. [0 (false)]
        -mult_index          <int>      index for the output multipole. Requires multipoles output format. [-1]
        -output              <filename> output file. If multipoles is enable will contain t/M amplitude phase. Otherwise t(s) h+ hx. ['waveform.dat']
```

If you want to specify a parfile

```
$ ./TEOBResum -p <parfile>
```

please see

`examples/example.par`

for an example for a BBH q=1 chi1=chi2=0.7.



