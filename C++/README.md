# Quick start 

## Requirements

- C++
- GSL/BLAS 

## How to compile

### Makefile

A simple Makefile is provided in the main directory. Just run 

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
$ ./configure --prefix=</abs/path/to/your/build/> CXXFLAGS="-std=c++11 -Ofast -Wall -fPIC"
$ make
$ make install
```

To build and install the python wrapper

```
$ cd python 
$ python setup.py build_ext --inplace 
```

Once the installation is finished, add 
/abs/path/to/your/build/
to the path

## How to run 

```
$ ./TEOBResum <options>
```
