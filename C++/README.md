# Quick start 

## Requirements

- C++ compiler (GNU, Intel)
- GSL/BLAS 

## How to compile

There are 2 options

- Simple Makefile
- Build that includes Python wrapper. NOTE THE BUILD WITH PYTHON IS POORLY TESTED AND MIGHT NOT WORK ON CERTAIN PLATFORMS

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
$ ./bootstrap.sh
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

Usage/help:

```
$ ./TEOBResumS.x 
```

The code can be run by specifying a parfile

```
$ ./TEOBResum -p <parfile>
```

or by specifying inputs of key quantities from command line. Please see `examples/`.

